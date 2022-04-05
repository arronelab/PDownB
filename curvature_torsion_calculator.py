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

def get_curvature(chunk):
    # get midpoints #
    mp1 = (chunk[1]['CA'].get_coord() + chunk[0]['CA'].get_coord())/2
    mp2 = (chunk[2]['CA'].get_coord() + chunk[1]['CA'].get_coord())/2
    mp3 = (chunk[3]['CA'].get_coord() + chunk[2]['CA'].get_coord())/2
    # get sin(theta) from ||u x v || = ||u||||v||sin(theta) #
    u = mp1 - mp2
    v = mp2 - mp3
    cross = np.cross(u,v)
    magcross = magnitude(cross)
    magu = magnitude(u)
    magv = magnitude(v)
    sintheta = magcross/(magu*magv)
    curvature = (2*np.absolute(sintheta))/magu
    return curvature

def get_torsion(chunk):
    u = chunk[1]['CA'].get_coord() - chunk[0]['CA'].get_coord()
    v = chunk[2]['CA'].get_coord() - chunk[1]['CA'].get_coord()
    w = chunk[3]['CA'].get_coord() - chunk[2]['CA'].get_coord()
    n1 = normalise(np.cross(u,v))
    n2 = normalise(np.cross(v,w))
    theta = np.arccos(np.dot(n1,n2))
    if np.dot(np.cross(n1,n2),v) < 0:
        theta*=-1
    seglength = (magnitude(u)+magnitude(v)+magnitude(w))/3
    torsion = (2/length)*np.sin(theta/2)))
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
                    chunked_chain = chunker(chain,3)
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
                            chunk_no+=1
                        except:
                            chunk_no+=1
                            pass
        except:
            pass
        output = pd.DataFrame(output_dicts_list)
        output.to_csv('curvature_torsion_analysis.csv', index=False)
        return output.head()
                    
                    
            






