"""
_ascension_

This module contains functions for translating Ct 434Bu Uniprot protein
ascension numbers to CT numbers.
"""
import csv

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
            'CT Number',
            oheaders[1], # Original protein name
            oheaders[16], # Fold change
            oheaders[17], # P Value
        ]

        # Grab the existing data we want 
        for row in reader:
            protein_name = row[1]
            fold_change = row[16]
            p_value = row[17]
            data.append([None, protein_name, fold_change, p_value])

    return headers, data

def translate_data(data):
    """
    Fill in the CT number data based on a parsed CTL number from a protein
    description.

    :param data: 2D list of rows. By row index:
      - 0 CT Number (is None before translation)
      - 1 Protein description -- parse the pretranslation CT number from this
      - 2 Fold Change (RifR/GspE)
      - 3 P-Value
    """
    for row in data:
        # Parse out the ctl number
        protein_desc = row[1]
        ctl = parse_ctl(protein_desc)

        # Translate to ct number
        row[0] = translate_ct_number(ctl)

def parse_ctl(protein_desc):
    pass

def translate_ct_number(ctl):
    pass


if __name__ == "__main__":
    headers, data = import_data('../data/ComparativeProteomics.csv')
    data = translate_data(data)
    import ipdb; ipdb.set_trace()
