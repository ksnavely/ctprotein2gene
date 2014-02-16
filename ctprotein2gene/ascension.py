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
            oheaders[1], # Original protein name
            oheaders[16], # Fold change
            oheaders[17], # P Value
        ]

        # Grab the existing data we want 
        for row in reader:
            protein_name = row[1]
            fold_change = row[16]
            p_value = row[17]
            data.extend([protein_name, fold_change, p_value])

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
    theaders = [].extend(headers).extend(
        [
            'Protein Prefix',
            'Protein Name',
            'CTL Number',
            'CT Number'
        ]
    )

    # Parse and convert the description to the desired data
    translated = [].extend(data)
    for row in translated:
        # Parse out the protein name and prefix
        protein_desc = row[1]
        prefix, name = parse_protein_name(protein_desc)

        # Use the name to find CTL, CT numbers
        ctl = get_ctl_number(name)
        ct = get_ct_number(ctl)
        row.extend([prefix, name, ctl, ct])

    return theaders, translated

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


if __name__ == "__main__":
    # Extract the data of interest from the Excel->CSV dump
    headers, data = import_data('../data/ComparativeProteomics.csv')

    # Get new headers and 2D list of orig data + protein prefix, protein name,
    # ctl number, and ct number
    headers, data = translate_data(headers, data)
    import ipdb; ipdb.set_trace()
