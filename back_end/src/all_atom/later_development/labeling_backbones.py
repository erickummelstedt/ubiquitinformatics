#from openbabel import pybel
import os
import pandas as pd
### this should all go before the bonding and atom changes...
## as bonding and atom changes delete atoms..

from scipy.spatial.transform import Rotation as R
import numpy as np
import pandas as pd
import math
import json
import copy
import sys

home_dir = os.path.expanduser('~')
local_path = home_dir + '/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/'

sys.path.insert(0, local_path + 'aa_functions_folder')
sys.path.insert(0, local_path + '.data')

def create_mol2_from_smiles(smiles_string):
    mol = pybel.readstring("smi", smiles_string)
    mol.addh()
    mol.make3D()
    #os.remove(smiles_string + ".txt")
    mol.write("mol2", ".data/" + smiles_string + ".mol2", overwrite=True) 

## add notes and change the column management so that the type of each column matches the data file... 
def pull_out_mol(df):
    molIndex = df[df['#']=="@<TRIPOS>MOLECULE"].index[0]
    atomIndex = df[df['#']=="@<TRIPOS>ATOM"].index[0]
    if df['#'].isin(["@<TRIPOS>UNITY_ATOM_ATTR"]).any():
        unityAtomIndex = df[df['#']=="@<TRIPOS>UNITY_ATOM_ATTR"].index[0]
    bondIndex = df[df['#']=="@<TRIPOS>BOND"].index[0]
    if df['#'].isin(["@<TRIPOS>SUBSTRUCTURE"]).any():
        subsIndex = df[df['#']=="@<TRIPOS>SUBSTRUCTURE"].index[0]
    else: 
        subsIndex = df.tail(1).index[0]+1
    ## input is df... output is a file with all the atoms
    molTRIPODdf = df.loc[molIndex+1:atomIndex-1]
    molTRIPODdf = molTRIPODdf.reset_index()
    if 'index' in molTRIPODdf.columns:
        molTRIPODdf = molTRIPODdf.drop('index',axis=1)
    molTRIPODdf
    return molTRIPODdf

def pull_out_atoms(df):
    molIndex = df[df['#']=="@<TRIPOS>MOLECULE"].index[0]
    atomIndex = df[df['#']=="@<TRIPOS>ATOM"].index[0]
    if df['#'].isin(["@<TRIPOS>UNITY_ATOM_ATTR"]).any():
        unityAtomIndex = df[df['#']=="@<TRIPOS>UNITY_ATOM_ATTR"].index[0]
    bondIndex = df[df['#']=="@<TRIPOS>BOND"].index[0]
    if df['#'].isin(["@<TRIPOS>SUBSTRUCTURE"]).any():
        subsIndex = df[df['#']=="@<TRIPOS>SUBSTRUCTURE"].index[0]
    else: 
        subsIndex = df.tail(1).index[0]+1
    ## input is df... output is a file with all the atoms
    
    ## in case there is unity atom atrr
    if df['#'].isin(["@<TRIPOS>UNITY_ATOM_ATTR"]).any():
        atomsTRIPODdf = df.loc[atomIndex+1:unityAtomIndex-1]
    else:
        atomsTRIPODdf = df.loc[atomIndex+1:bondIndex-1]
    for x in range(100):
        atomsTRIPODdf.loc[:,'#'] = atomsTRIPODdf.loc[:,'#'].astype(str).str.replace('  ',' ')
    atomsTRIPODdf.loc[:,'#'] = atomsTRIPODdf.loc[:,'#'].str.strip(' ')
    atomsTRIPODdf = atomsTRIPODdf['#'].str.split(' ', expand=True)
    ## check the len of the columsn in the mol2 file... sometimes open babel doesnt add it
    if len(atomsTRIPODdf.columns) == 9:
        atomsTRIPODdf[9]= 'None'
    atomColumnNames = 'atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge', 'status_bit'
    atomsTRIPODdf.columns = atomColumnNames
    atomsTRIPODdf = atomsTRIPODdf.reset_index()
    if 'index' in atomsTRIPODdf.columns:
        atomsTRIPODdf = atomsTRIPODdf.drop('index',axis=1)
    return atomsTRIPODdf

def pull_out_bonds(df):
    molIndex = df[df['#']=="@<TRIPOS>MOLECULE"].index[0]
    atomIndex = df[df['#']=="@<TRIPOS>ATOM"].index[0]
    if df['#'].isin(["@<TRIPOS>UNITY_ATOM_ATTR"]).any():
        unityAtomIndex = df[df['#']=="@<TRIPOS>UNITY_ATOM_ATTR"].index[0]
    bondIndex = df[df['#']=="@<TRIPOS>BOND"].index[0]
    if df['#'].isin(["@<TRIPOS>SUBSTRUCTURE"]).any():
        subsIndex = df[df['#']=="@<TRIPOS>SUBSTRUCTURE"].index[0]
    else: 
        subsIndex = df.tail(1).index[0]+1
    ## input is df... output is a file with all the bonds
    bondsTRIPODdf = df.loc[bondIndex+1:subsIndex-1]
    for x in range(100):
        bondsTRIPODdf.loc[:,'#'] = bondsTRIPODdf.loc[:,'#'].astype(str).str.replace('  ',' ')
    #TRIPODbondsdf = TRIPODbondsdf['#'].str.split(' ', expand=True)
    bondsTRIPODdf.loc[:,'#'] = bondsTRIPODdf.loc[:,'#'].str.strip(' ')
    bondsTRIPODdf = bondsTRIPODdf['#'].str.split(' ', expand=True)
    if len(bondsTRIPODdf.columns) == 4:
        bondsTRIPODdf[5] =''
    bondColumnNames = 'bond_id', 'origin_atom_id', 'target_atom_id', 'bond_type', 'status_bits'
    bondsTRIPODdf.columns = bondColumnNames
    bondsTRIPODdf = bondsTRIPODdf.reset_index()
    if 'index' in bondsTRIPODdf.columns:
        bondsTRIPODdf = bondsTRIPODdf.drop('index',axis=1)
    bondsTRIPODdf
    return bondsTRIPODdf

## df = download from mol2 file..
def get_BONDS_ATOMS(df):
    ### track ubi number and conjugation side... 
    # initial things
    allAtoms = pull_out_atoms(df)
    allBonds = pull_out_bonds(df)

    allAtoms.loc[:,'atom_id'] = allAtoms.loc[:,'atom_id'].astype(int)
    allAtoms.loc[:,'subst_id'] = allAtoms.loc[:,'subst_id'].astype(int)
    allAtoms.loc[:,'x']= allAtoms.loc[:,'x'].astype(float)
    allAtoms.loc[:,'y']= allAtoms.loc[:,'y'].astype(float)
    allAtoms.loc[:,'z']= allAtoms.loc[:,'z'].astype(float)

    allBonds.loc[:,'origin_atom_id']= allBonds.loc[:,'origin_atom_id'].astype(int)
    allBonds.loc[:,'target_atom_id']= allBonds.loc[:,'target_atom_id'].astype(int)
    allBonds.loc[:,'bond_id']= allBonds.loc[:,'bond_id'].astype(int)
    return allAtoms, allBonds



## https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space

### this is just a rotation matrix... not a scaling one...
### def give both vectors the same starting point...
def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v): #if not all zeros then 
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    else:
        return np.eye(3) #cross of all zeros only occurs on identical directions
    
## function that finds atom_types and atom_ids of the bonded atoms..
## input; atom_id, atoms and bonds
## output; atom_id's of bonded atoms, atom_types of bonded atoms
def find_bonded_atoms(atom_id, atoms, bonds): 
    ## create empty list of atom_ids
    neighbouring_atom_ids = []
    bonds_of_atomid = bonds[(bonds['origin_atom_id']== atom_id)|(bonds['target_atom_id']==atom_id)]
    neighbouring_atom_ids = list(bonds_of_atomid.target_atom_id)+list(bonds_of_atomid.origin_atom_id)
    neighbouring_atom_ids = [ele for ele in neighbouring_atom_ids if ele != atom_id]
    ## if there is not 2H and one C break...
    neighbouring_atom_types = []
    for id in neighbouring_atom_ids:
        atom_type_index = atoms[atoms['atom_id']==id].atom_type.index[0]
        infunch_atom_type = atoms.at[atom_type_index, 'atom_type']
        ## if it has 
        neighbouring_atom_types.append(infunch_atom_type)
    return neighbouring_atom_ids, neighbouring_atom_types 

## insert concat
def reformat_atom_bond_ids(atoms, bonds): 
    atoms['new_atom_id'] = range(1,len(atoms)+1)
    atomID_Dict = pd.Series(atoms.new_atom_id.values,index=atoms.atom_id).to_dict()

    #Using exhaustive map() dict
    atoms["atom_id"]=atoms["atom_id"].map(atomID_Dict)
    bonds['bond_id'] = range(1,len(bonds)+1)

    bonds["origin_atom_id"]=bonds["origin_atom_id"].map(atomID_Dict)
    bonds["target_atom_id"]=bonds["target_atom_id"].map(atomID_Dict)
    atoms=atoms.reset_index()
    bonds=bonds.reset_index()
    atoms=atoms.drop('index', axis=1)
    atoms=atoms.drop('new_atom_id', axis=1)
    bonds=bonds.drop('index', axis=1)

    atoms.loc[:,'x'] = atoms.loc[:,'x'].astype(float)
    atoms.loc[:,'y'] = atoms.loc[:,'y'].astype(float)
    atoms.loc[:,'z'] = atoms.loc[:,'z'].astype(float)

    atoms['x'] = atoms['x'].apply(lambda x: format(x, '.4f'))
    atoms['y'] = atoms['y'].apply(lambda x: format(x, '.4f'))
    atoms['z'] = atoms['z'].apply(lambda x: format(x, '.4f'))
    return atoms, bonds


## only call this function if the N-terminus has 2 H's
def mol2_Nterminus_H3addition(N_atom_id, atoms, bonds):
    ## copy base_atoms & base_bonds
    base_atoms = atoms.copy()
    base_bonds = bonds.copy()

    # ALERT: this will be input in the def variables
    N_atoms_index = base_atoms[base_atoms['atom_id'] == N_atom_id].index[0]

    ## FINDING COORDINATES 
    ## find position of two H's
    ## take the midpoint between them from H3 which we will call HE
    nitro_atom_ids, nitro_atom_types = find_bonded_atoms(N_atom_id, atoms, bonds)
    
    ## find c3 index IN LIST and use it to get atom_id for C.3
    c3_list_index = nitro_atom_types.index('C.3')
    c3_atom_id = nitro_atom_ids[c3_list_index]
    c3_atoms_index = base_atoms[base_atoms['atom_id'] == c3_atom_id].index[0]
    ## must delete from both lists
    del nitro_atom_types[c3_list_index]
    del nitro_atom_ids[c3_list_index]

    ## find H1 index IN LIST and use it to get atom_id for H
    H1_list_index = nitro_atom_types.index('H')
    H1_atom_id = nitro_atom_ids[H1_list_index]
    H1_atoms_index = base_atoms[base_atoms['atom_id'] == H1_atom_id].index[0]
    ## must delete from both lists
    del nitro_atom_types[H1_list_index]
    del nitro_atom_ids[H1_list_index]

    ## find H1 index IN LIST  and use it to get atom_id for H
    H2_list_index = nitro_atom_types.index('H')
    H2_atom_id = nitro_atom_ids[H2_list_index]
    H2_atoms_index = base_atoms[base_atoms['atom_id'] == H2_atom_id].index[0]
    ## must delete from both lists
    del nitro_atom_types[H2_list_index]
    del nitro_atom_ids[H2_list_index]

    ## pull bonding information for N-H bonds
    N_H1_index = base_bonds[((base_bonds['origin_atom_id'] == N_atom_id)&(base_bonds['target_atom_id'] == H1_atom_id))|((base_bonds['origin_atom_id'] == H1_atom_id)&(base_bonds['target_atom_id'] == N_atom_id))].index[0]
    N_H2_index = base_bonds[((base_bonds['origin_atom_id'] == N_atom_id)&(base_bonds['target_atom_id'] == H2_atom_id))|((base_bonds['origin_atom_id'] == H2_atom_id)&(base_bonds['target_atom_id'] == N_atom_id))].index[0]

    # find the coordinates of the atoms of interest
    x_N ,y_N, z_N  = base_atoms.at[N_atoms_index, 'x'],base_atoms.at[N_atoms_index, 'y'],base_atoms.at[N_atoms_index, 'z']
    x_C3 ,y_C3, z_C3  = base_atoms.at[c3_atoms_index, 'x'],base_atoms.at[c3_atoms_index, 'y'],base_atoms.at[c3_atoms_index, 'z']
    x_H1 ,y_H1, z_H1  = base_atoms.at[H1_atoms_index, 'x'],base_atoms.at[H1_atoms_index, 'y'],base_atoms.at[H1_atoms_index, 'z']
    x_H2 ,y_H2, z_H2  = base_atoms.at[H2_atoms_index, 'x'],base_atoms.at[H2_atoms_index, 'y'],base_atoms.at[H2_atoms_index, 'z']
    

    # build coordinates for H3 which you can just use the minus of the the N-C3 vector 
    # find length of bond... 
    x_mid_H1_H2 = (x_H1+x_H2)/2
    y_mid_H1_H2 = (y_H1+y_H2)/2
    z_mid_H1_H2 = (z_H1+z_H2)/2

    ## take the vector to minus from the opposite from 
    x_adjust_vec = x_N - x_mid_H1_H2
    y_adjust_vec = y_N - y_mid_H1_H2
    z_adjust_vec = z_N - z_mid_H1_H2

    # find the location of the coordinates opposite to the c3 divide this vector by two...
    x_C3_180 = x_N + (x_N-x_C3)
    y_C3_180 = y_N + (y_N-y_C3)
    z_C3_180 = z_N + (z_N-z_C3)

    ## addd the adjustment vector to the coordinates opposite to the c3
    x_H3 = x_C3_180 + x_adjust_vec
    y_H3 = y_C3_180 + y_adjust_vec
    z_H3 = z_C3_180 + z_adjust_vec

    # build new entries in bonds and atoms...    
    ## make the new atom id 0...
    # change subst_id - '1', 
    # subst_name - '<1>', 
    # find all the info and add it in
    # change  index = 17, add in after second H 
    HE_atom_type = 'HE'
    HE_subst_id = base_atoms.at[H2_atoms_index, 'subst_id']
    HE_subst_name = base_atoms.at[H2_atoms_index, 'subst_name']
    HE_charge = 0
    HE_status_bit = 'None'
    HE_atom_id = 0
    base_atoms = pd.DataFrame(np.insert(base_atoms.values, (H2_atoms_index+1), values=[HE_atom_id, HE_atom_type, x_H3, y_H3, z_H3, 'H', HE_subst_id, HE_subst_name, HE_charge, HE_status_bit], axis=0), columns=base_atoms.columns)
    # add in after second H..
    N_HE_bond_id = 0
    N_HE_bond_type = 1
    N_HE_status_bits = ''
    base_bonds = pd.DataFrame(np.insert(base_bonds.values, (N_H2_index+1), values=[N_HE_bond_id, N_atom_id, HE_atom_id, N_HE_bond_type, N_HE_status_bits], axis=0), columns=base_bonds.columns)
    
    ## reformat atom_ids and bond_ids 
    base_atoms, base_bonds = reformat_atom_bond_ids(base_atoms, base_bonds)

    return base_atoms, base_bonds

#probably good to round certrain columns to particular decimal places
## 4 is all that is needed 

### build mol2 file and check it out in pymol...
### reverse engineer the mol2 file... 
## check it out in pymol 
def push_togther_atoms(df):
    ## input is df... output is a file with all the atoms
    df1 = df.copy()
    df2 = pd.DataFrame()
    df2['atom_id'] = df1['atom_id']
    df2["C"] = ""
    for x in list(df1.columns):
        df1[x] = df1[x].astype(str)
        if x == 'atom_type' or x == 'atom_name' or x == 'subst_name':
            df1[x] = df1[x].apply(lambda x: "{:<9}".format(x[:9]))
        if x == 'atom_id':
            df1[x] = df1[x].apply(lambda x: "{:>7}".format(x[:7]))
        else:    
            df1[x] = df1[x].apply(lambda x: "{:>9}".format(x[:9]))
        df2["C"] = df2["C"] + ' ' + df1[x]
    
    df2= df2.drop('atom_id', axis = 1)
    return df2

def push_togther_bonds(df):
    ## input is df... output is a file with all the atoms
    df1 = df.copy()
    df2 = pd.DataFrame()
    df2['bond_id'] = df1['bond_id']
    df2["C"] = ""
    for x in list(df1.columns):
        df1[x] = df1[x].astype(str)
        if x == 'bond_type':
            df1[x] = df1[x].apply(lambda x: "{:<9}".format(x[:9]))
        else:    
            df1[x] = df1[x].apply(lambda x: "{:>9}".format(x[:9]))
        df2["C"] = df2["C"] + ' ' + df1[x]
    
    df2= df2.drop('bond_id', axis = 1)
    return df2

def write_mol2_local(final_ubiAllAtoms, final_ubiAllBonds, PDB_ID):
    maxAtomID = max(final_ubiAllAtoms['atom_id'])
    maxBondID = max(final_ubiAllBonds['bond_id'])
    atoms_combined = push_togther_atoms(final_ubiAllAtoms)
    bonds_combined = push_togther_bonds(final_ubiAllBonds)

    mol2Header = pd.DataFrame({'C': ['#',
                                    '# Created by the Jeffrey Bode Group of ETH Zurich',
                                    '# Software Developed by ERIC KUMMELSTEDT', 
                                    '#', 
                                    '#', 
                                    '@<TRIPOS>MOLECULE', 
                                    str(PDB_ID), 
                                    str(maxAtomID) + ' ' + str(maxBondID),
                                    'SMALL',
                                    'NO_CHARGES',
                                    '@<TRIPOS>ATOM'
                                    ]})

    pandas_for_mol = pd.concat([mol2Header, atoms_combined, pd.DataFrame({'C': ['@<TRIPOS>BOND']}), bonds_combined, pd.DataFrame({'C': ['@<TRIPOS>SUBSTRUCTURE', '0 pdb 0 ****']}                                                                                                                                                                                                                                                                                             )],ignore_index=True)
    np.savetxt(local_path + ".data/mol2_files/amino_acids/" + PDB_ID + ".txt", pandas_for_mol, delimiter="", fmt="%s") 
    os.rename(local_path + ".data/mol2_files/amino_acids/" + PDB_ID + ".txt",local_path + ".data/mol2_files/amino_acids/" + PDB_ID + ".mol2")
    ## add the number of bonds and number of atoms... 3rd line of the mol2 file...
    ### max of each value....  

## function that finds atom_types and atom_ids of the bonded atoms..
## input; atom_id, atoms and bonds
## output; atom_id's of bonded atoms, atom_types of bonded atoms
def find_bonded_atoms(atom_id, atoms, bonds): 
    ## create empty list of atom_ids
    neighbouring_atom_ids = []
    bonds_of_atomid = bonds[(bonds['origin_atom_id']== atom_id)|(bonds['target_atom_id']==atom_id)]
    neighbouring_atom_ids = list(bonds_of_atomid.target_atom_id)+list(bonds_of_atomid.origin_atom_id)
    neighbouring_atom_ids = [ele for ele in neighbouring_atom_ids if ele != atom_id]
    ## if there is not 2H and one C break...
    neighbouring_atom_types = []
    for id in neighbouring_atom_ids:
        atom_type_index = atoms[atoms['atom_id']==id].atom_type.index[0]
        infunch_atom_type = atoms.at[atom_type_index, 'atom_type']
        ## if it has 
        neighbouring_atom_types.append(infunch_atom_type)
    return neighbouring_atom_ids, neighbouring_atom_types 

   
## just take the first H in bonding
## just take the first H in bonding
def check_atoms_on_nitro(input_atom_types, aa_position = 'individual'):

    ## should probably test this lol
    ## second half of all of this is proline....things with two c.3
    if aa_position == 'individual':
        desired_atom_types_list= [['C.3', 'H', 'H'], ['C.3', 'H', 'H', 'H'], ['C.3', 'C.3', 'H'], ['C.3', 'C.3', 'H', 'H']] ## create for loop
    elif aa_position  == 'n_terminus':
        desired_atom_types_list= [['C.3', 'H', 'H'], ['C.3', 'H', 'H', 'H'], ['C.3', 'C.3', 'H'], ['C.3', 'C.3', 'H', 'H']] ## create for loop
    ## when the nitrogen is bonded it doesn't have the extra hydrogen
    elif aa_position  == 'c_terminus':
        desired_atom_types_list= [['C.3', 'C.2', 'H'], ['C.3', 'C.3', 'C.2']]  ## create for loop
    elif aa_position  == 'double_bonded_peptide':
        desired_atom_types_list= [['C.3', 'C.2', 'H'], ['C.3', 'C.3', 'C.2']]  ## create for loop

    for desired_atom_types in desired_atom_types_list: 
        input_copy = input_atom_types.copy()
        desired_copy = desired_atom_types.copy()
        ## check len of both lists is the same
        if len(desired_atom_types) == len(input_atom_types):
            for i in desired_atom_types:
                removable_atom = desired_copy.pop()
                if removable_atom in input_copy:
                    input_copy.remove(removable_atom)
        else:
            outcome = False

        if input_copy == [] and input_atom_types != []:
            outcome = True
            break
        else:
            outcome = False
    return outcome == True

def check_atoms_on_c3(input_atom_types, aa_position = 'individual'):
    ## should probably test this lol
    if aa_position == 'individual':
        desired_atom_types_list= [['C.2', 'H', 'N.3'], ['C.2', 'H', 'N.4']] ## create for loop
    elif aa_position  == 'n_terminus':
        desired_atom_types_list= [['C.2', 'H', 'N.3'], ['C.2', 'H', 'N.4']] ## create for loop
    elif aa_position  == 'c_terminus':
        desired_atom_types_list= [['C.2', 'H', 'N.am']]  ## create for loop
    elif aa_position  == 'double_bonded_peptide':
        desired_atom_types_list= [['C.2', 'H', 'N.am']]  ## create for loop

    for desired_atom_types in desired_atom_types_list:     
        input_copy = input_atom_types.copy()
        desired_copy = desired_atom_types.copy()
        ## check len of both lists is the same
        if (len(desired_atom_types)+1) == len(input_atom_types):
            for i in desired_atom_types:
                removable_atom = desired_copy.pop()
                if removable_atom in input_copy:
                    input_copy.remove(removable_atom)
        else:
            outcome = False

        if len(input_copy) == 1 and input_atom_types != []:
            outcome = True
            break
        else:
            outcome = False
    return outcome == True

def check_atoms_on_c2(input_atom_types, aa_position = 'individual'):
    ## should probably test this lol
    if aa_position == 'individual':
        desired_atom_types_list= [['O.2', 'O.3', 'C.3'], ['O.co2', 'O.co2', 'O.co2', 'C.3']] ## create for loop
    elif aa_position  == 'n_terminus':
        desired_atom_types_list= [['O.2', 'N.am', 'C.3'], ['O.2', 'O.3', 'C.3']] ## the second is complicated...
    elif aa_position  == 'c_terminus':
        desired_atom_types_list= [['O.2', 'O.3', 'C.3'], ['O.co2', 'O.co2', 'O.co2', 'C.3']]  ## create for loop
    elif aa_position  == 'double_bonded_peptide':
        desired_atom_types_list= [['O.2', 'N.am', 'C.3']]  ## create for loop

    for desired_atom_types in desired_atom_types_list:     
        input_copy = input_atom_types.copy()
        desired_copy = desired_atom_types.copy()
        ## check len of both lists is the same
        if len(desired_atom_types) == len(input_atom_types):
            for i in desired_atom_types:
                removable_atom = desired_copy.pop()
                if removable_atom in input_copy:
                    input_copy.remove(removable_atom)
        else:
            outcome = False

        if input_copy == [] and input_atom_types != []:
            outcome = True
            break
        else:
            outcome = False
    return outcome == True

## nitro
def fix_bonding_nitrogen(nitrogen, atoms, bonds):
    nitro_atom_ids, nitro_atom_types = find_bonded_atoms(nitrogen, atoms, bonds)
    nitrogen_H_atom_id_list = []
    ## fix this....
    while 'H' in nitro_atom_types:
        H_index = nitro_atom_types.index('H')
        H_atom_id = nitro_atom_ids[H_index]
        nitrogen_H_atom_id_list.append(H_atom_id)
        del nitro_atom_types[H_index]
        del nitro_atom_ids[H_index]

    if len(nitrogen_H_atom_id_list) == 3: 
        ## while H1 is 
        Ha_removed_atom_id = nitrogen_H_atom_id_list[-3]
        ## H2 is just removed 
        Hb_removed_atom_id = nitrogen_H_atom_id_list[-2]
        ## H3 turns in HZ1
        Hc_changed_atom_id = nitrogen_H_atom_id_list[-1]
        ## remove Hb is there are 3 
        atoms, bonds = remove_atom_from_bonds_atoms(Hb_removed_atom_id, atoms, bonds)

    elif len(nitrogen_H_atom_id_list) == 2: 
        Ha_removed_atom_id = nitrogen_H_atom_id_list[-2]
        ## H3 turns in HZ1 
        Hc_changed_atom_id = nitrogen_H_atom_id_list[-1]
    
    Hc_changed_index = atoms[atoms['atom_id']==Hc_changed_atom_id].index[0]
    atoms.loc[Hc_changed_index,'atom_name'] = 'HZ'

    return Ha_removed_atom_id, atoms, bonds


## carry out general changes here that will always bonds that are broken bonds that are formed... 
## then create any extra changes you might need to create...
def remove_atom_from_bonds_atoms(atom_id, atoms, bonds): 
    ## bond removal
    removalable_bonds_index_df = bonds[(bonds['origin_atom_id'] == atom_id)|(bonds['target_atom_id'] == atom_id)]
    removalable_bonds_index_list = list(removalable_bonds_index_df.index)
    for idx in removalable_bonds_index_list: 
        bonds = bonds.drop(index=idx)
    
    ## atom removal 
    atom_index = atoms[atoms['atom_id'] == atom_id].index[0]
    atoms = atoms.drop(index=atom_index)

    return atoms, bonds


## if the amino acid is bonded then the atom types will be different.. for N and C.2 
## N will have one H and one C.2 and one C.3
## can also be ['C.3', 'H', 'C.2']
## can also be ['C.3', 'H', 'H', 'C.2']
## proline-esque amino acid..
## c-terminus amino acid
## to check if it is the N or C terminus check the the subst_id number....
## input should be single 
## find/fix peptide backbone...
## the C terminus and N terminus should always look the same
## N terminus has three H and one N 
## DONE C terminus has an OXT/O.2 and O/O.3
## peptide has N amide bond and O/O.2 which means OXT should change to O

## finding max/min subst_id happens outside find_peptide_backbone
## adapt for proline... 
 
def find_peptide_backbone(atoms, bonds, subst_id = 1, aa_position = 'individual'):
    nitrogens_df = atoms[((atoms['atom_type']=='N.3')|(atoms['atom_type']=='N.4')|(atoms['atom_type']=='N.am'))&(atoms['subst_id']==subst_id)]
    nitrogens_atom_id_list = list(nitrogens_df.atom_id)
    backbone_atom_id_list = []
    nitrogen_H_atom_id_list = []
    o3_atom_id_list = []
    co2_atom_id_list = []
    peptide_backboneQ = False
    stop = False
    print(nitrogens_atom_id_list)
    for nitro_atom_id in nitrogens_atom_id_list: 
        co2_atom_id_list = []
        backbone_atom_id_list = []
        nitrogen_H_atom_id_list = []
        o3_atom_id_list = []
        nitro_atom_ids, nitro_atom_types = find_bonded_atoms(nitro_atom_id, atoms, bonds)
        ## check_list_of_atoms gives a boolean
        ## if there is ['C.3', 'H', 'H']
        ## can also be ['C.3', 'H', 'H', 'H']
        print('nitro_atom_types' + str(nitro_atom_types))
        c3_atom_id_list = []
        if check_atoms_on_nitro(nitro_atom_types, aa_position):
            ## figure out backbone_atom_id_list            
            # pull the c3 indexes and 
            # work yourself through the list of C.3 
            backbone_atom_id_list.append(nitro_atom_id)
            while ('C.3' in nitro_atom_types):
                c3_index = nitro_atom_types.index('C.3')
                c3_atom_id_list.append(nitro_atom_ids[c3_index])
                del nitro_atom_types[c3_index]
                del nitro_atom_ids[c3_index]
                print('c3_atom_id_list' + str(c3_atom_id_list))

            #counter = 0
            #while 'H' in nitro_atom_types:
            #    H_index = nitro_atom_types.index('H') + counter
            #    H_atom_id = nitro_atom_ids[H_index]
            #    nitrogen_H_atom_id_list.append(H_atom_id)
            #    nitro_atom_types.remove('H')
            #    counter = counter + 1 
            
            while 'H' in nitro_atom_types:
                H_index = nitro_atom_types.index('H')
                H_atom_id = nitro_atom_ids[H_index]
                nitrogen_H_atom_id_list.append(H_atom_id)
                del nitro_atom_types[H_index]
                del nitro_atom_ids[H_index]
                
            ## copy the backbone so it can be used each time in the for loop
            backbone_atom_id_list_copy = backbone_atom_id_list.copy()

            ## for loop of C3
            for c3_atom_id in c3_atom_id_list:
                backbone_atom_id_list = backbone_atom_id_list_copy.copy()
                c3_atom_ids, c3_atom_types = find_bonded_atoms(c3_atom_id, atoms, bonds)
                print('c3_atom_types' + str(c3_atom_types))
                if check_atoms_on_c3(c3_atom_types, aa_position):
                    backbone_atom_id_list.append(c3_atom_id)
                    ## the c2 index becomes complicated its a list because there can be two c.2's
                    c3_atom_types_copy = c3_atom_types.copy()
                    
                    ## can turn this into while loop
                    print('c3_atom_types_copy' + str(c3_atom_types_copy))
                    counter = 0
                    while 'C.2' in c3_atom_types_copy:
                        ## why + counter
                        c2_index = c3_atom_types_copy.index('C.2') + counter
                        c2_atom_id = c3_atom_ids[c2_index]
                        c2_atom_ids, c2_atom_types = find_bonded_atoms(c2_atom_id, atoms, bonds)
                        print('c2_atom_types before change: ' + str(c2_atom_types))
                        if check_atoms_on_c2(c2_atom_types, aa_position):
                            ## THIS PART OF THE FUNCTION CHANGES THE CTERMINUS TO THE DESIRED CTERMINUS...
                            ## REMOVING O.co2 and putting in its place the O.2 & O.3 
                            ## in O.co2 change it to O.2 -- OXT then in the reaction you lose the O.2
                            ## change it to O.3 - 0 
                            ## repeat c2_atom_ids, c2_atom_types = find_bonded_atoms(c2_atom_id, atoms, bonds)
                            ## to see how many c.O2 remain...
                            ## the C terminus and N terminus should always look the same
                            ## N terminus has three H and one N 
                            ## C terminus has an OXT/O.2 and O/O.3
                            ## peptide has N amide bond and O/O.2 which means OXT should change to O
                            counter1 = 0
                            while 'O.co2' in c2_atom_types:
                                c2_atom_ids, c2_atom_types = find_bonded_atoms(c2_atom_id, atoms, bonds)
                                if counter1 == 0: 
                                    ## here we pull the index from the list not the atoms dataframe
                                    Oco2_index = c2_atom_types.index('O.co2')
                                    Oco2_atom_id = c2_atom_ids[Oco2_index]
                                    Oco2_atom_ids, Oco2_atom_types = find_bonded_atoms(Oco2_atom_id, atoms, bonds)
                                    ## if len = 2 and count() = 2 delete second bonded atom...
                                    ## bonding changes fixed...
                                    ## if there is a double bond to the same atom then drop the double bond..
                                    if (len(Oco2_atom_ids) == 2) & (len(set(Oco2_atom_ids)) ==1):
                                        ar_data_row = bonds[(bonds['origin_atom_id'] == c2_atom_id) & (bonds['target_atom_id']== Oco2_atom_id) & (bonds['bond_type']=='ar')]
                                        ar_index_to_drop = ar_data_row.index[0]
                                        bonds = bonds.drop([ar_index_to_drop])
                                    ## change the atom_id value...
                                    ## to change the atom_id fiund the index of the Oco2_atom_id
                                    Oco2_atom_index = atoms[atoms.loc[:,'atom_id'] == Oco2_atom_id].index[0]
                                    atoms.at[Oco2_atom_index, 'atom_type'] = 'O.2'
                                    atoms.at[Oco2_atom_index, 'atom_name'] = 'OXT'
                                    ## check the bonded atoms... for H
                                    ## first find the doubly bonded O.co2 and delete
                                if counter1 == 1: 
                                    ## here we pull the index from the list not the atoms dataframe
                                    Oco2_index = c2_atom_types.index('O.co2')
                                    Oco2_atom_id = c2_atom_ids[Oco2_index]
                                    Oco2_atom_ids, Oco2_atom_types = find_bonded_atoms(Oco2_atom_id, atoms, bonds)
                                    ## if len = 2 and count() = 2 delete second bonded atom...
                                    ## bonding changes fixed...
                                    ## if there is a double bond to the same atom then drop the double bond..
                                    if (len(Oco2_atom_ids) == 2) & (len(set(Oco2_atom_ids)) ==1):
                                        ar_data_row = bonds[(bonds['origin_atom_id'] == c2_atom_id) & (bonds['target_atom_id']== Oco2_atom_id) & (bonds['bond_type']=='ar')]
                                        ar_index_to_drop = ar_data_row.index[0]
                                        bonds = bonds.drop([ar_index_to_drop])
                                    ## change the atom_id value...
                                    ## to change the atom_id fiund the index of the Oco2_atom_id
                                    Oco2_atom_index = atoms[atoms.loc[:,'atom_id'] == Oco2_atom_id].index[0]
                                    atoms.at[Oco2_atom_index, 'atom_type'] = 'O.3'
                                    atoms.at[Oco2_atom_index, 'atom_name'] = 'O'
                                    ## check the bonded atoms... for H
                                    ## first find the doubly bonded O.co2 and delete
                                if ('O.co2' in c2_atom_types) == False: 
                                    break
                                counter1 = counter1 + 1 


                        c2_atom_ids, c2_atom_types = find_bonded_atoms(c2_atom_id, atoms, bonds)
                        print('c2_atom_types after change: ' + str(c2_atom_types))
                        if check_atoms_on_c2(c2_atom_types, aa_position):
                            backbone_atom_id_list.append(c2_atom_id)
                            if 'O.2' in c2_atom_types:
                                o2_index = c2_atom_types.index('O.2')
                                o2_atom_id = c2_atom_ids[o2_index]
                                backbone_atom_id_list.append(o2_atom_id)
                            if 'O.3' in c2_atom_types:
                                o3_index = c2_atom_types.index('O.3')
                                o3_atom_id = c2_atom_ids[o3_index]
                                o3_atom_id_list.append(o3_atom_id)
                            stop = True
                            break
                        c3_atom_types_copy.remove('C.2')
                        ## test that the break works
                
                    if stop == True:
                        break
                    print('break_test1')
                
            if stop == True:
                break
            print('break_test')

    if len(backbone_atom_id_list) == 4:
        ## pull the nitrogen atom and make sure it has 3 H if not add the 3 N
        ## if statement to check it is n_terminus
        for i in backbone_atom_id_list:
            new_index = atoms[atoms['atom_id'] == i].index[0]
            atoms.loc[new_index,'status_bit'] = 'BACKBONE'
            peptide_backboneQ = True
        
        ## this function goes last as the atom_ids will change 
        if (aa_position == 'individual')|(aa_position == 'n_terminus'): 
            nitro_atom_id = backbone_atom_id_list[0]
            nitro_atom_ids, nitro_atom_types = find_bonded_atoms(nitro_atom_id, atoms, bonds)
            ## throw error if nitro atoms does not equal 2 or 3...
            if nitro_atom_types.count('H') == 2: 
                print('protonating N-terminus')
                atoms, bonds = mol2_Nterminus_H3addition(nitro_atom_id, atoms, bonds)
                ## calls itself but this time it should have the extra H so the if statement wil be skipped
                atoms, bonds, result_dictionary = find_peptide_backbone(atoms, bonds, subst_id, aa_position)

                peptide_backboneQ = result_dictionary['peptide_backboneQ']
                backbone_atom_id_list = result_dictionary['backbone_atom_id_list']
                nitrogen_H_atom_id_list = result_dictionary['nitrogen_H_atom_id_list']
                o3_atom_id_list = result_dictionary['o3_atom_id_list']
                co2_atom_id_list = result_dictionary['co2_atom_id_list']
    
    print('final backbone_atom_id_list: ' + str(backbone_atom_id_list))
    result_dictionary = {'peptide_backboneQ': peptide_backboneQ == True, 
                         'backbone_atom_id_list': backbone_atom_id_list, 
                         'nitrogen_H_atom_id_list': nitrogen_H_atom_id_list, 
                         'o3_atom_id_list': o3_atom_id_list, 
                         'co2_atom_id_list': co2_atom_id_list}
        
    ## other than atoms and bonds change the output to a dictionary
    return atoms, bonds, result_dictionary