#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:25:08 2022

@author: ccck67
"""
import os
import pandas as pd
import numpy as np
from Bio import PDB

protein_dir = 'PDB'
root = 'PDB'

def magnitude(vector):
    return np.sqrt(vector.dot(vector))

def normalise(vector):
    return vector/magnitude(vector)

amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", 
               "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", 
               "TYR", "VAL"]

def strip_water(chain):
    return [resi for resi in chain if resi.get_resname() in amino_acids]

### chunker splits a list in subsections of consecutive elements
### e.g. chunker([1,2,3,4,5], 3) returns [[1,2,3], [2,3,4], [3,4,5]]

def chunker(list_in, chunk_size):
    output = []
    pos = 0
    for item in list_in:
        chunk = list_in[pos:pos+chunk_size]
        if len(chunk) == chunk_size:
            output.append(tuple(list_in[pos:pos+chunk_size]))
            pos+=1
        else:
            pass
    return output

### We define the curvature as the scalar product of consecutive unit
### tangent vectors T_i . T_{i+1}

def get_curvature(chunk):
    Ti = normalise(chunk[1]['CA'].get_coord() - chunk[0]['CA'].get_coord())
    Tip1 = normalise(chunk[2]['CA'].get_coord() - chunk[1]['CA'].get_coord())
    curvature = np.dot(Ti,Tip1)
    return curvature

### We define the torsion as the sclar product of the plane normals given by
### the cross product of pairs of consecutive unit tangent vectors.
### i.e (T_i x T_{i+1}) . (T_{i+1} x T_{i+2})

def get_torsion(chunk):
    Ti = normalise(chunk[1]['CA'].get_coord() - chunk[0]['CA'].get_coord())
    Tip1 = normalise(chunk[2]['CA'].get_coord() - chunk[1]['CA'].get_coord())
    Tip2 = normalise(chunk[3]['CA'].get_coord() - chunk[2]['CA'].get_coord())
    Ni = np.cross(Ti,Tip1)
    Nip1 = np.cross(Tip1,Tip2)
    torsion = np.dot(Ni,Nip1)
    return torsion

### The following will run through the PDB files, and compute the curvature 
### and torsion along the length of each protein chain.

def get_ct_distribution(): 
    output_dicts_list = []
    total_chains = 0
    for fin in os.listdir(protein_dir):
        try:
            name = fin.split('//')[-1][:-4]
            parser = PDB.PDBParser(QUIET=True)  ### Warnings are silenced - careful
            structure = parser.get_structure(name, os.path.join(root,fin))
            for model in structure:
                for chain in model:
                    total_chains += 1
                    chain_name = name + '_chain_' + chain.get_id()
                    chain = strip_water(chain)
                    chunked_chain = chunker(chain,4)
                    chunk_no = 0
                    for subsection in chunked_chain:
                        try:
                            to_add = {
                                'Chain' : chain_name,
                                'Subsection' : chunk_no,
                                'Curvature' : get_curvature(subsection),
                                'Torsion' : get_torsion(subsection)
                                }
                            output_dicts_list.append(to_add)
                            print(to_add)
                            chunk_no+=1
                        except:
                            chunk_no+=1
                            pass
        except:
            pass
        output = pd.DataFrame(output_dicts_list)
        output.to_csv('curvature_torsion_analysis.csv', index=False)
        return output.head()
    
print(get_ct_distribution())
                    
                    
            






