"""
_ascension_

This module contains functions for translating Ct 434Bu Uniprot protein
ascension numbers to CT numbers.
"""
from copy import deepcopy
import csv

import BeautifulSoup
import json
import requests

EB_PROTEOMICS_FILE = '../data/EB_Proteomics.csv'

def import_eb_proteomics_data(filename):
    # Map column name to header index
    headers = {}
    data = []

    # Open the input csv and pull out the fields we want
    # Make a header and column for our translated CT number
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile)

        # Get the first row, the headers, and build the output headers
        oheaders = reader.next()

        # Order the same as the data construction
        headers = [
            'CTL##',
            'Peptide Matches',
            'Fold Change (RifR/GspE)',
            'P-value'
        ]

        row_lookup = {
            'CTL##': oheaders.index('CTL##'),
            'Peptide Matches': oheaders.index('Peptide Matches'),
            'Fold Change (RifR/GspE)': oheaders.index('Fold Change (RifR/GspE)'),
            'P-value': oheaders.index('P-value')
        }

        # Grab the existing data we want 
        for row in reader:
            ctl = row[row_lookup['CTL##']]
            peptide_matches = row[row_lookup['Peptide Matches']]
            p_value = row[row_lookup['P-value']]
            fold_change = row[row_lookup['Fold Change (RifR/GspE)']]
            data.append([ctl, peptide_matches, fold_change, p_value])

    return headers, data


def process_eb_proteomics_data(data):
    thresholds = {
        'Peptide Matches': lambda x: float(x) >= 2,
        'Fold Change (RifR/GspE)': lambda x: abs(float(x)) >= 1.5,
        'P-value': lambda x: float(x) <= 0.05
    }
    out = []

    for d in data:
        ctl, peptide_matches, fold_change, p_value = d
        check_sum = sum([
            thresholds['Peptide Matches'](peptide_matches),
            thresholds['Fold Change (RifR/GspE)'](fold_change),
            thresholds['P-value'](p_value)
        ])
        if check_sum >= 2:
            out.append(d)

    return out


def write_reduced_csv(filename, headers, data):
    with open(filename, 'wb') as f:
        f.write('{0}\n'.format(', '.join(headers)))
        for d in data:
            f.write('{0}\n'.format(', '.join(d)))


def klee_dump(password, start_ctl, max_ctl):
    ctls = range(start_ctl, max_ctl + 1)
    bad_ctls = []

    for c in ctls:
        # Attempt to fetch full sequence
        seq = None
        try:
            seq = search_klee(c, first_seventy=False)
            print('Grabbed {ctl} {seq}...'.format(ctl=c, seq=seq))

            doc = {
                '_id': 'CTL00%03d' % c,
                'ctl': 'CTL00%03d' % c,
                'seq': '{seq}'.format(seq=seq),
            }

            response = requests.post(
                'https://ksnavely.cloudant.com/ctl_sequence',
                auth=('ksnavely', password),
                data=json.dumps(doc),
                headers={'Content-Type': 'application/json'}
            )
            response.raise_for_status()
            print('Saved {ctl}, resp: {resp}'.format(ctl=c, resp=response.text))

        except Exception as ex:
            bad_ctls.append((c, repr(ex)))
            print('Exception caught: {ex}'.format(ex=ex))
            continue

    with open('bad_ctls.txt', 'a') as f:
        for ctl, repr_ex in bad_ctls:
            f.write('{0}, {1}\n'.format(ctl, repr_ex))


def search_klee(ctl, first_seventy=True):
    KLEE_URL = 'http://www.genome.jp/dbget-bin/www_bget?ctb:CTL0%03d' % ctl

    print('Fetching {url}'.format(url=KLEE_URL))
    response = requests.get(
        KLEE_URL
    )
    response.raise_for_status()

    webpage = BeautifulSoup.BeautifulSoup(
        response.text
    )

    results = webpage.findAll(name='nobr', text='AA seq')
    if not results:
        raise Exception('No AA seq nobr found for {0}, resp: {1}'.format(ctl, response.text))

    if len(results) > 1:
        raise Exception('Results too long! len: {0}, resp: {1}'.format(
                len(results),
                response.text
            )
        )

    aa_seq_nobr = results[0]
    try:
        seq = aa_seq_nobr.parent.parent.parent.td.text.split(' aa')[1]
        print('Parsed seq {0}...'.format(seq[:10]))
    except:
        raise Exception('Unable to parse seq')

    if first_seventy:
        seq = seq[:70]

    return seq

def import_data(filename):
    # Output
    data = []

    # Open the input csv and pull out the fields we want
    # Make a header and column for our translated CT number
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile)

        # Get the first row, the headers, and build the output headers
        oheaders = reader.next()
        headers = [
            oheaders[1], # Original protein name
            oheaders[16], # Fold change
            oheaders[17], # P Value
        ]

        # Grab the existing data we want 
        for row in reader:
            protein_name = row[1]
            fold_change = row[16]
            p_value = row[17]
            data.append([protein_name, fold_change, p_value])

    return headers, data

def translate_data(headers, data):
    """
    :returns: headers and data for the original data plus additional protein
    prefix, protein name, ctl number, and ct number information.

    :param data: 2D list of rows. By row index:
      - 0 Protein description -- parse the name->CTL->CT number from this
      - 1 Fold Change (RifR/GspE)
      - 2 P-Value
    """
    # Extend the original headers
    theaders = deepcopy(headers)
    theaders.extend(
        [
            'Protein Prefix',
            'Protein Name',
            'CTL Number',
            'CT Number'
        ]
    )

    # Parse and convert the description to the desired data
    translated = []
    for row in data:
        # Parse out the protein name and prefix
        protein_desc = row[0]
        try:
            prefix, name = parse_protein_name(protein_desc)
        except ParseError:
            # Leave out mis-parses -- controls
            continue

        # Use the name to find CTL, CT numbers
        ctl = None #get_ctl_number(name)
        ct = None #get_ct_number(ctl)

        translated.append(
            row + [prefix, name, ctl, ct]
        )

    return theaders, translated

class ParseError(Exception):
    pass

def parse_protein_name(protein_desc):
    """
    Get the prefix and name of the protein that is listed at the start of the
    protein_desc.

    :param protein_desc: A protein description is a ~human readable string that
        begins with a descriptor followed by other text.
        
        The descriptor consists of a prefix, pipe ('|'), and the protein name.
        Ex: gb|AGJ65242.2

    :returns: The prefix and 'AGJ...' name of the protein as a tuple:
        (prefix, name)
    """
    protein_name = protein_desc.split(' ')[0]
    # Some controls won't parse, no worries
    if '|' not in protein_name:
        msg = 'Misparse: {protein_desc}'.format(protein_desc=protein_desc)
        print(msg)
        raise ParseError(msg)

    prefix = protein_name.split('|')[0]
    name = protein_name.split('|')[1]
    return prefix, name

def get_ctl_number(name):
    """
    Get the CTL number from the protein name.
    """
    pass

def get_ct_number(ctl):
    """
    Given the ctl number, return the ct number.
    """
    pass

def print_name_info(data):
    """
    Print name info for AGJ proteins
    """
    names = []
    for row in data:
        name = row[4]
        if name.startswith('AGJ'):
            names.append(row[4])

    names.sort()
    print("First AGJ protein: {first}".format(first=names[0]))
    print("Last AGJ protein: {last}".format(last=names[-1]))
    print("Non-AGJ proteins not considered.")


if __name__ == "__main__":
#   # Extract the data of interest from the Excel->CSV dump
#   headers, data = import_data('../data/ComparativeProteomics.csv')

#   # Get new headers and 2D list of orig data + protein prefix, protein name,
#   # ctl number, and ct number
#   headers, data = translate_data(headers, data)
#   print_name_info(data)
#   import ipdb; ipdb.set_trace()
    eb_headers, eb_data = import_eb_proteomics_data(EB_PROTEOMICS_FILE)
    reduced_eb_data = process_eb_proteomics_data(eb_data)
    write_reduced_csv('../data/Reduced_EB_Proteomics.csv', eb_headers, reduced_eb_data)
    import ipdb; ipdb.set_trace()
