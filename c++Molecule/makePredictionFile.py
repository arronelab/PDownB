#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 10:20:51 2022

@author: ccck67
"""
import urllib.request
import os
import re
import json
import numpy as np
import pandas as pd
import ast
from Bio import PDB, SeqIO
from Bio.SeqUtils import seq1 

def write_sh_file(project_name = input("Enter project name: "), 
                  preprocessed_data_dir = "newFitData/human_SMARCAL1",
                  fingerprint_path = "newFitData/human_SMARCAL1/fingerPrint.dat", 
                  coords_path = "newFitData/human_SMARCAL1/coordinates.dat", 
                  contacts_path = None, 
                  fixed_sections_path = None, 
                  fit_steps = input("Enter number of changes from input structure: "),
                  no_mols = input("Enter number of outputs: ")
                  ):
    """Writes the config shell script based on parameters given by the user on data upload.

    Args:
        project_name (str): [description]
        preprocessed_data_dir (str): [description]
        fingerprint_path (str): [description]
        coords_path (str, optional): [description]. Defaults to None.
        contacts_path (str, optional): [description]. Defaults to None.
        fixed_sections_path (str, optional): [description]. Defaults to None.
        fit_steps (int, optional): [description]. Defaults to 1000.

    Returns:
        output_filepath (str): path to the finished shell script
    """
    output_filepath = project_name + '_config.sh'
    output_dir = os.path.join(preprocessed_data_dir,project_name)
    with open(output_filepath, 'w+') as fout:
    	fout.write('#!/bin/bash')
    	fout.write('\nSequenceFile={}'.format(fingerprint_path))
    	if coords_path is not None:
    		fout.write('\ninitialCoordsFile={}'.format(coords_path))
    	else:
    		fout.write('\ninitialCoordsFile=none')
    	if contacts_path is not None:
    		fout.write('\npairedPredictions={}'.format(contacts_path))
    	else:
    		fout.write('\npairedPredictions=none')
    	if fixed_sections_path is not None:
    		fout.write('\nfixedsections={}'.format(fixed_sections_path))
    	else:
    		fout.write('\nfixedsections=none')
        # no options for crystal symmetry or or hydro cover in prototype
    	fout.write('\ncrystalSymmetry=none')
    	fout.write('\nwithinMonomerHydroCover=none')
    	fout.write('\nbetweenMonomerHydroCover=none')
    	fout.write('\nmaxNoFitSteps={}'.format(str(fit_steps)))
    	fout.write('\n\nmkdir {}'.format(output_dir))
    	fout.write('\n\n\nfor i in {1..'+str(no_mols)+'}')
    	fout.write('\ndo')
    	fout.write('\n    predictStructure $SequenceFile '
    	'$initialCoordsFile $pairedPredictions $fixedsections $crystalSymmetry '
    	'$withinMonomerHydroCover $betweenMonomerHydroCover '
    	'$maxNoFitSteps {}/mol$i.dat'.format(output_dir))
    	fout.write('\ndone')
    return output_filepath

if __name__ == "__main__":
    write_sh_file()
    