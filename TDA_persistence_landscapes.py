# Created by Aurelia Li, Adsorption and Advanced Materials Group (aam.ceb.cam.ac.uk),
# led by David Fairen-Jimenez from the Department of Chemical Engineering and
# Biotechnology, University of Cambridge,
# Thank you to Seth Wiggin from the Cambridge Crystallographic Data Centre, Cambridge,
# for testing and for the useful feedback.

"""
Step 2 - Compute the persistent landscapes

Input: - a GCD list e.g. the GCD list output by TDA_structural_data_preparation.py
       - the .csv files containing the structures coordinates, also output from TDA_structural_data_preparation.py
Output: - an .npy file for each structure, where the persistence data is stored
        - a results.csv file containing the persistence calculated for all the structures in a human-readable format. 
        Each row corresponds to a structure, and the persistence data is recorded in the format (Betti number, (birth, death))
        - a times.csv containing the computation time for each structure

"""
import numpy as np
import gudhi as gd
import os
import time
import gudhi.representations
import csv

# Prepare dictionaries for easier browsing of diagrams (results for diagrams,
# times for computation time and simplex_trees for simplex_trees)
results = {}
times = {}
simplex_trees = {}

# Read in the relevant list of structures
path_to_GCD = 'all_potential_cages_selected.gcd' # MODIFY THIS IF NECESSARY
gcd_list = open(path_to_GCD, 'r').read()
refcodes = gcd_list.split("\n")

for refcode in refcodes:
    if os.path.exists(refcode+".csv"):
        print('Loading', refcode)
        # Extract the corresponding coordinates from CSV
        coordinates = np.genfromtxt(refcode+".csv", delimiter=",")

        print('Now trying Rips')
        # calculate persistence (I know some structures such as ABIGEY don't work
        # hence the try) (this is actually because of memory issues)
        try:
            # Keep track of computation time
            start=time.time()
            # Perform TDA
            Rips_complex_sample = gd.RipsComplex(points = coordinates, max_edge_length=0.6)
            Rips_simplex_tree_sample = Rips_complex_sample.create_simplex_tree(max_dimension=3)
            diag_Rips = Rips_simplex_tree_sample.persistence()
            stop=time.time()
            # store persistent diagrams in results
            results[refcode]=diag_Rips
            times[refcode]=stop-start
            simplex_trees[refcode] = Rips_simplex_tree_sample
            print('Persistence for', refcode, 'calculated')
        except:
            print('Something went wrong with', refcode)
            results[refcode]='Error'
            times[refcode]='Error'
        # Calculate landscapes, L1 for Betti 1 and L2 for Betti 2
        LS = gd.representations.Landscape(resolution=100) # MODIFY RESOLUTION AS DESIRED
        print('LS calculated')
        L1 = LS.fit_transform([Rips_simplex_tree_sample.persistence_intervals_in_dimension(1)])
        print('L1 calculated')
        try:
            L2 = LS.fit_transform([Rips_simplex_tree_sample.persistence_intervals_in_dimension(2)])
            print('L2 calculated')
            L = L2
        except:
            L = np.concatenate((L1[0], L2[0]))
            print('L computed')
    else:
        results[refcode]='xyz missing'
        times[refcode]='NA'
        print('xyz file not found for', refcode)
    w = csv.writer(open("results.csv", "a", newline = ''))
    w.writerow([refcode, results[refcode]])
    t= csv.writer(open("times.csv", "a", newline = ''))
    t.writerow([refcode, times[refcode]])
    # landscape is saved as an NPY
    with open(refcode+".npy", 'wb') as f:
        np.save(f,L)
    print("Trees array saved for", refcode)
