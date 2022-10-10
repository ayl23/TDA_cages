# Created by Aurelia Li, Adsorption and Advanced Materials Group (aam.ceb.cam.ac.uk),
# led by David Fairen-Jimenez from the Department of Chemical Engineering and
# Biotechnology, University of Cambridge.
# Thank you to Seth Wiggin from the Cambridge Crystallographic Data Centre, Cambridge,
# for testing and for the useful feedback.

"""
Step 1 - extract the fractional coordinates of the potential cages from the CSD

Input: - a GCD list (text format for list of structures used by the CCDC) obtained with e.g. ConQuest
Output: - the CIFs corresponding to the structures in the GCD list, after some preliminary checks
        - CSV files containing their fractional coordinates
        - a GCD list confirming the list of structures for which CIFs have been output

To run this script directly in the CSD Python API command shell:
- Specify the location of the GCD file and the working directory (marked with # MODIFY THIS)
- Open the CSD Python API command line
- Navigate to the desired folder (using e.g. os)
- exec(open('TDA_structural_data_preparation.py').read())

"""

## Part Ia - check the structures identified with ConQuest and obtain the clean CIFs

# NB: this requires access to the CSD Python API
# see https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/primer.html for first time use
import ccdc.molecule
from ccdc.io import EntryReader, CrystalWriter
import csv
import os

# Read list of structures identified with ConQuest
# List of structures are saved as GCD (.gcd) text files in ConQuest
path_to_GCD = 'GCD.gcd' # MODIFY THIS
csd_reader = EntryReader('CSD')
MOF_list_entries = EntryReader(path_to_GCD)

# To keep track of numbers
potential_cage_list = []
total_structures = len(MOF_list_entries)
potential_structures = 0

# This loops checks that the component that could be potentially a cage has organic parts and is not fully linear
for entry in MOF_list_entries:
    atom_list=[]
    refcode = entry.identifier
    # Check that the heaviest weight component has organic parts, and if so,
    # keep it. Otherwise, if it's not organic at all, look at the other components
    # (there should be one other in general) and check that it is organometallic.
    if 'Atom(C' in str(entry.molecule.heaviest_component.atoms):
        cage = entry.molecule.heaviest_component
    else:
        for idx, component in enumerate(list(entry.molecule.components)):
            if entry.molecule.components[idx] != entry.molecule.heaviest_component:
                if 'Atom(C' in str(entry.molecule.components[idx].atoms) and entry.molecule.components[idx].is_organometallic:
                    cage  = entry.molecule.components[idx]
    # Check that at least one atom is part of a 'cycle'
    for atom in cage.atoms:
        if atom.is_cyclic is True:
            atom_list.append(atom)
    if len(atom_list) != 0:
        potential_cage_list.append(refcode)
        potential_structures += 1
        total_structures += 1
        print("This is the", potential_structures, "th potential structure out of", total_structures)
        print("This structure is", refcode)
        entry.crystal.molecule = cage
        # Write out the corresponding cif for checking (not necessary)
        with CrystalWriter(refcode+'.cif') as cryst_writer:
            cryst_writer.write(entry.crystal)
        # Write out the final list of structures
        with open('all_potential_cages_selected.gcd', 'a+') as f:
            f.write(refcode)
            f.write("\n")
            f.close()

## Part Ib - extract the fractional coordinates from the CIFs

# Specify the directory we want the csv files to be written to
directory = 'URPATH'# MODIFY THIS

# Function to remove parenthesis in cifs
def remove_parenthesis(list_of_strings):
    '''
    Takes a list of strings that contains numbers and something like this:([0-9][0-9]).
    Remove the parentheses and what's between them
    Return a list of floats
    '''
    for i, element in enumerate(list_of_strings):
        if '(' in element:
            first = element.find('(')
            second = element.find(')')
            element = element[:-(second - first +1)]
            element = float(element)
            list_of_strings[i] = element
    return list_of_strings

# Write coordinates into CSV
for filename in os.listdir(directory):
    if filename.endswith(".cif"):
        cif_reader = EntryReader(filename)
        for cif in cif_reader:
            cif = cif_reader[0] # need to specify which datablock in the CIF we're reading
            try:
                if cif.has_3d_structure:
                    print('Extracting coordinates for ', filename)
                    atom_x_coordinate = remove_parenthesis(cif.attributes['_atom_site_fract_x'])
                    atom_y_coordinate = remove_parenthesis(cif.attributes['_atom_site_fract_y'])
                    atom_z_coordinate = remove_parenthesis(cif.attributes['_atom_site_fract_z'])
                    rows = zip(atom_x_coordinate, atom_y_coordinate, atom_z_coordinate)
                    with open(filename[:-4]+".csv", "w") as f:
                        writer = csv.writer(f)
                        for row in rows:
                            writer.writerow(row)
            except RuntimeError:
                print('Error when extracting coordinates from', filename)
