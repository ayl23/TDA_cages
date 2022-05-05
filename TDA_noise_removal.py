# Created by Aurelia Li, Adsorption and Advanced Materials Group (aam.ceb.cam.ac.uk),
# led by David Fairen-Jimenez from the Department of Chemical Engineering and
# Biotechnology, University of Cambridge.

"""
Step 3 - remove noisy structures

This script contains functions to identify noisy structures. 
Use find_similar() to identify structures similar to a given structure.
The "heatmaps" must be computed prior to using find_similar().

Input: - the results.csv file obtained previously in step 2      -  
"""
import gudhi as gd
import csv
import pandas as pd

# Process previously obtained results
results = pd.read_csv('results.csv', header = None)
results.columns=['refcode', 'persistence']
results = results.drop(results[results['persistence'] == 'xyz missing'].index)
results = results.dropna()
results = dict(zip(list(results.refcode), list(results.persistence)))

# Check that all the entries are 'complete'.
for refcode in results.keys():
    if results[refcode][-1] != ']' :
        print(refcode, 'incomplete')
        find_index = len(results[refcode])-1
        match = ', (0, (0'
        segment = results[refcode][-1]
        while match not in segment:
            find_index -= 1
            segment_new = results[refcode][find_index] + segment
            segment = segment_new
        results[refcode] = results[refcode][:find_index-1] + ')]'

# The persistence lists are read in as strings so we need to
# convert them to lists first
# But, the presence of 'inf' in the list make it hard for ast or jason to
# convert them to lists so we need to replace the places where inf appears first

for refcode in results.keys():
    try:
        string_to_convert = results[refcode]
        string_to_convert = string_to_convert.replace("(0, (0.0, inf)), ", "")
        results[refcode]=eval(string_to_convert)
    except:
        print('error with', refcode)

# look for structure wihout betti 2 numbers
def check_betti_2(results, refcode):
    '''
    This function takes a persistence diagram and checks if there are results
    corresponding to a betti number of 2.
    '''
    if results[refcode][0][0] < 2:
        return False
    else:
        return True

# Gather these structures without betti 2s into a new dict and also save the list in a csv file
no_betti_2 = {}
for refcode in results.keys():
    if check_betti_2(results, refcode) == False:
        w = csv.writer(open("no_betti_2.csv", "a"))
        w.writerow([refcode])
        no_betti_2[refcode]=results[refcode]

# We want to compare the persistence diagrams cooresponding to betti numbers 1 and 2
# We need to extract the diagrams as a list of lists type
def persistence_to_compare(results, refcode, betti):
    '''
    results - dictionary of persistence diagram previously obtained in the format of
              keys: refcodes, values: persistence diagram (list of tuples (betti, (birth, death)))
    refcode - structure we are looking at
    betti - the number we want to look at, 1 or 2
    This function returns a persistence diagram of betti number for a given structure by
    preparing the diagram in the format of a list of lists [[birth, death], [birth, death]]
    '''
    persistence = []
    i = 0
    while results[refcode][i][0] != betti and i <= len(results[refcode]):
        i += 1
    while results[refcode][i][0] == betti and i < len(results[refcode]):
        persistence.append(list(list(results[refcode][i])[1]))
        i += 1
    return persistence

# Load results obtained previsouly
structures = results.keys()

# Prepare dataframes with all the bottleneck distances
heatmap_data_1 = pd.DataFrame(index = structures, columns = structures)
heatmap_data_2 = pd.DataFrame(index = structures, columns = structures)

for refcode in structures:
    for refcode_to_compare in structures:
        if refcode not in no_betti_2.keys() and refcode_to_compare not in no_betti_2.keys():
            print(refcode, refcode_to_compare)
            print('Comparing', refcode, 'and', refcode_to_compare)
            heatmap_data_2.loc[refcode][refcode_to_compare] = gd.bottleneck_distance(persistence_to_compare(results, refcode, 2), persistence_to_compare(results, refcode_to_compare, 2))
            heatmap_data_1.loc[refcode][refcode_to_compare] = gd.bottleneck_distance(persistence_to_compare(results, refcode, 1), persistence_to_compare(results, refcode_to_compare, 1))

# Function to identify similar structures
def find_similar(heatmap, refcode, threshold):
    '''
    This function transforms the bottleneck distances stored in heatmap into
    another similarity scale:
        1 - structure is very similar (corresponds to the structure itself)
        0 - structure is completely dissimilar
    Then, the function looks at the second most similar structure apart
    from itself, and looks for the closest values to that second most similar
    structure
    Threshold corresponds to the % to which we want the structures to be similar
    to that second most similar structure.
    Returns a dict of similar structures with name: bottleneck distance
    '''
    similar = {}
    refcode_col = heatmap.loc[refcode]
    # normalising the scale of values and ranging similarity from 1 (very similar) to 0 (least similar) and take out the value corresponding to refcode
    norm_values = refcode_col.sort_values().transform(lambda x:(refcode_col.max()-x)/refcode_col.max())
    indices = norm_values.index.values.tolist()
    # after being sorted, the second element has the highest similarity value (doesn't mean necessarily they are similar in absolute terms)
    best_max = norm_values.loc[indices[1]]
    similar[indices[0]] = heatmap.loc[refcode, indices[0]]
    # take out the element corresponding to best_max
    norm_values = norm_values.iloc[1:]
    for index in indices[1:]:
        if  norm_values.loc[index] >= threshold*best_max:
            similar[index]=heatmap.loc[refcode][index]

    return similar