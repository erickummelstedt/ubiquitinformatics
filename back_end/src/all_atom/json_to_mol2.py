import os
import sys
import math
import json
import copy
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.spatial.transform import Rotation as R

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))

from src.main import \
    validate_protein_keys, \
    validate_branching_sites, \
    find_branching_site, \
    process_current_protein, \
    add_max_chain_number

from src.utils.logging_utils import \
    log_protein_details, \
    log_branching_details, \
    log_end_of_branching, \
    log_end_of_protein   

from src.utils.utils import \
    convert_json_to_dict

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

def rotation_transformation_of_new_PG(ubiquitinAllAtoms, ubiquitinAllBonds, single_AllAtoms, single_AllBonds, ubiquitinConjugationSite, conjugatingAtomsStart, conjugatingAtomsEnd, protecting_group='SMAC'):
    # define the new ubiquitin... and find the lysine of choice... 
    ## copy so that the value isnt assigned...
    ## this is be a starting place.. 

    ## adaptable atoms is a copy of the ubiquitin of choice...
    ## this will change to be a subset of the original ubiquitin...
    ## adaptable ubi is the ubi unit being added...
    ## this will not be the case...
    full_conjugatingAtoms = ubiquitinAllAtoms.copy()
    full_conjugatingUbiBonds = ubiquitinAllBonds.copy()

    ## find the lysines of interest     
    conjugatingAtoms = full_conjugatingAtoms[(conjugatingAtomsStart <= ubiquitinAllAtoms['subst_id']) &  (conjugatingAtomsEnd  >= ubiquitinAllAtoms['subst_id'])]
    conjugatingUbiBonds = full_conjugatingUbiBonds

    ## adaptable atoms is a copy of the ubiquitin of choice...
    ## this will change to be a subset of the original ubiquitin...
    ## adaptable ubi is the ubi unit being added...
    newAtoms = single_AllAtoms.copy()
    newBonds = single_AllBonds.copy()

    ## FIRST     origin the centre of mass...
    ## a transform matrix that aligns two vectors... does it do scaling and transformation..??
    # minus the centre of mass from all the positions, 
    # this zeros the protein on the centre of mass and the coordinates become the origin... and the rotation occurs around the origin..
    ## TRANSFORM ALL THE COORDINATES 
    ## end of day... lets get the numbers reordered. done....
    ## start looking at the transformation...  

    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...

    ## ALERTTTT: the centre of mass here should be of the ubiquitin that is going to be added 
    ## when you transfer coordiantes there may be misalignment.. be careful here..
    ## pull rotational and transformational matrices from the coordaintes of the ubi of interest 
    # and apply to the new ubi
    ## you can do this from the origin..
    ## makes things easier... just find the coordinates of interest...
    ## try this with the trimer..


    # ubi of interest... find the aa numbers that matter...
    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...
    conjugatingAtoms['c'] = conjugatingAtoms['atom_name'].str[:1]
    atomID_mass_Dict = {'N':'14', 'C': '12', 'O':'16', 'H':'1', 'S': '32', 'L': '17'}
    conjugatingAtoms['c'] = conjugatingAtoms['c'].map(atomID_mass_Dict)
    conjugatingAtoms.loc[:,'c'] = conjugatingAtoms.loc[:,'c'].astype(int)
    conjugatingAtoms.loc[:,'x'] = conjugatingAtoms.loc[:,'x'].astype(float)
    conjugatingAtoms.loc[:,'y'] = conjugatingAtoms.loc[:,'y'].astype(float)
    conjugatingAtoms.loc[:,'z'] = conjugatingAtoms.loc[:,'z'].astype(float)
    ## finding centre of mass
    conjugatingAtoms['cx'] = conjugatingAtoms['c']*conjugatingAtoms['x']
    conjugatingAtoms['cy'] = conjugatingAtoms['c']*conjugatingAtoms['y']
    conjugatingAtoms['cz'] = conjugatingAtoms['c']*conjugatingAtoms['z']
    x_COM=conjugatingAtoms['cx'].sum()/conjugatingAtoms['c'].sum()
    y_COM=conjugatingAtoms['cy'].sum()/conjugatingAtoms['c'].sum()
    z_COM=conjugatingAtoms['cz'].sum()/conjugatingAtoms['c'].sum()
    ## drop the new columns
    conjugatingAtoms=conjugatingAtoms.drop('c', axis=1)
    conjugatingAtoms=conjugatingAtoms.drop('cx', axis=1)
    conjugatingAtoms=conjugatingAtoms.drop('cy', axis=1)
    conjugatingAtoms=conjugatingAtoms.drop('cz', axis=1)

    # new ubiquitin... 
    ## this could also be calculated from the conjugating ubiquitin
    ## there might be a slight difference in centre of mass... as the OXT is missing...
    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...
    newAtoms['c'] = newAtoms['atom_name'].str[:1]
    atomID_mass_Dict = {'N':'14', 'C': '12', 'O':'16', 'H':'1', 'S': '32', 'L': '17'}
    newAtoms.loc[:,'c'] = newAtoms.loc[:,'c'].map(atomID_mass_Dict)
    newAtoms.loc[:,'c'] = newAtoms.loc[:,'c'].astype(int)
    newAtoms.loc[:,'x'] = newAtoms.loc[:,'x'].astype(float)
    newAtoms.loc[:,'y'] = newAtoms.loc[:,'y'].astype(float)
    newAtoms.loc[:,'z'] = newAtoms.loc[:,'z'].astype(float)
    
    newAtoms['cx'] = newAtoms['c']*newAtoms['x']
    newAtoms['cy'] = newAtoms['c']*newAtoms['y']
    newAtoms['cz'] = newAtoms['c']*newAtoms['z']
    x_NEW_COM=newAtoms['cx'].sum()/newAtoms['c'].sum()
    y_NEW_COM=newAtoms['cy'].sum()/newAtoms['c'].sum()
    z_NEW_COM=newAtoms['cz'].sum()/newAtoms['c'].sum()
    ## drop the new columns
    newAtoms=newAtoms.drop('c', axis=1)
    newAtoms=newAtoms.drop('cx', axis=1)
    newAtoms=newAtoms.drop('cy', axis=1)
    newAtoms=newAtoms.drop('cz', axis=1)

    ## working glycine belongs to the new ubiquitin... never changes...
    ## need to find centre of mass for this molecule...
    ### find the sc terminal glycine...
    index_PG_Cl= newAtoms[newAtoms['atom_type'] == 'LG'].index[0]
    atomID_PG_Cl = newAtoms.at[index_PG_Cl, 'atom_id']
    # check the bonding... 
    Cl_C_index = newBonds[(newBonds['origin_atom_id'] == atomID_PG_Cl)|(newBonds['target_atom_id']== atomID_PG_Cl)].index[0]
    # atomID_PG_C equals the atom id from the value that is not Cl
    if newBonds.at[Cl_C_index, 'origin_atom_id'] == atomID_PG_Cl: 
        atomID_PG_C = newBonds.at[Cl_C_index, 'target_atom_id']
    else: 
        atomID_PG_C = newBonds.at[Cl_C_index, 'origin_atom_id']
    index_PG_C= newAtoms[newAtoms['atom_id'] == atomID_PG_C].index[0]

    ## may need to look at the vector of each bond... so the C_OXT 
    ## and the NZ3 -- HZ3
    ### working lysine belongs to the posiition of conjugation...
    ## this can be on one of the ubiquitins of the multimer
    workingLysine = conjugatingAtoms[conjugatingAtoms['subst_id'] == ubiquitinConjugationSite]
    index_LYS_NZ = workingLysine[workingLysine['atom_name'] == 'NZ'].index[0]
    index_LYS_HZ3 = workingLysine[workingLysine['atom_name'] == 'HZ3'].index[0]
    index_LYS_C = workingLysine[workingLysine['atom_name'] == 'C'].index[0]

    ## Centre the Coordinates on the centre of mass (COM)...
    ## probably dont need this...
    conjugatingAtoms['x1'] = conjugatingAtoms['x']-x_COM
    conjugatingAtoms['y1'] = conjugatingAtoms['y']-y_COM
    conjugatingAtoms['z1'] = conjugatingAtoms['z']-z_COM
    conjugatingAtoms['nparray'] = conjugatingAtoms[['x1', 'y1', 'z1']].values.tolist()
    
    ## Centre the Coordinates on the centre of mass (COM)...
    newAtoms['x1'] = newAtoms['x']-x_NEW_COM
    newAtoms['y1'] = newAtoms['y']-y_NEW_COM
    newAtoms['z1'] = newAtoms['z']-z_NEW_COM
    newAtoms['nparray'] = newAtoms[['x1', 'y1', 'z1']].values.tolist()

    ## find the coordinates of glysine carbon and lysine carbon
    ## need to decide whether this is done by bond vector or done by the centre of mass to central atom...
    x_PG_Cl,y_PG_Cl,z_PG_Cl  = newAtoms.at[index_PG_Cl, 'x'],newAtoms.at[index_PG_Cl, 'y'],newAtoms.at[index_PG_Cl, 'z']
    x_PG_C,y_PG_C,z_PG_C  = newAtoms.at[index_PG_C, 'x'],newAtoms.at[index_PG_C, 'y'],newAtoms.at[index_PG_C, 'z']
    
    x_LYS_HZ3,y_LYS_HZ3,z_LYS_HZ3  = workingLysine.at[index_LYS_HZ3, 'x'],workingLysine.at[index_LYS_HZ3, 'y'],workingLysine.at[index_LYS_HZ3, 'z']
    x_LYS_NZ,y_LYS_NZ,z_LYS_NZ  = workingLysine.at[index_LYS_NZ, 'x'],workingLysine.at[index_LYS_NZ, 'y'],workingLysine.at[index_LYS_NZ, 'z']
    x_LYS_C,y_LYS_C,z_LYS_C  = workingLysine.at[index_LYS_C, 'x'],workingLysine.at[index_LYS_C, 'y'],workingLysine.at[index_LYS_C, 'z']

    ## SECOND    find your rotation matrix.. 
    ## maybe need to minus the points here... flip the coordinates...
    ### here flip the glysine coordinates as we are looking for a smaller rotation... 
    ## we are rotating around the centree of mass
    ## because we are going to the one 180 position of the final position...
    
    ## two different centre of mass... glycine moves v little... 
    ## these vectors are important!!!!
    # bond replacement
    vec1 = np.array([x_PG_C-x_PG_Cl, y_PG_C-y_PG_Cl, z_PG_C-z_PG_Cl])
    # carbon to centre of mass replacement
    # vec1 = np.array([x_NEW_COM-x_GLY_C, y_NEW_COM-y_GLY_C, z_NEW_COM-z_GLY_C])

    ## lysine moves much more...
    # bond replacement
    vec2 = np.array([x_LYS_HZ3-x_LYS_NZ, y_LYS_HZ3-y_LYS_NZ, z_LYS_HZ3-z_LYS_NZ])
    # carbon to centre of mass replacement
    # vec2 = np.array([x_LYS_C-x_COM, y_LYS_C-y_COM, z_LYS_C-z_COM])


    ## find the rotation matrix...
    rotationVector = rotation_matrix_from_vectors(vec1, vec2)
    r = R.from_matrix(rotationVector)
    r.as_matrix().shape
    recreatedVector = r.apply(vec1)

    ## check that the rotation occured around the centre of mass... 
    ## a scaling matrix should equalise the two vectors.. vec2 and vec4 here 
    ## so this works beautifully....
    vec = recreatedVector
    testvec1 = vec[0]**2 + vec[1]**2 + vec[2]**2
    testvec2 = vec2[0]**2 + vec2[1]**2 + vec2[2]**2
    ratio = math.sqrt(testvec1)/math.sqrt(testvec2)
    recreatedVector.round(decimals=2) == (vec2*ratio).round(decimals=2)

    ## THIRD    ROTATIONAL TRANSFORMATION OF ALL POINTS on the new ubiquitin...
    newAtoms['np_rotated'] = newAtoms['nparray'].map(lambda x: r.apply(np.array(x)))
    newAtoms['x2'] = newAtoms['np_rotated'].map(lambda x: x[0])
    newAtoms['y2'] = newAtoms['np_rotated'].map(lambda x: x[1])
    newAtoms['z2'] = newAtoms['np_rotated'].map(lambda x: x[2])

    ## FOURTH    return the points to the centre of mass of the ubiquitin of interest...
    ## apply centre of mass transformation from ubi of interest to the new ubi
    newAtoms['x3'] = newAtoms['x2']+x_COM
    newAtoms['y3'] = newAtoms['y2']+y_COM
    newAtoms['z3'] = newAtoms['z2']+z_COM

    ### FIFTH    move the glysine on C of new ubiquitin to the HZ3 position...

    ## and then push all the points to where the HZ3 was... from the new coordinate...
    ## FIFTH    TRANSFORM ALL THE COORDINATES to the HZ3 from the C OF THE ROTATED GLYCINE...
    ### get the coordinates of HZ3 on the lysine of interest
    # use adapted as the HZ3 is still there... adapted HZ3 to startubi C...

    # lysine .. HZ3 === this disappears
    ## find the index of the HZ3 on the lysine of the starting ubiquitin... index to be deleted
    ## can still pul HZ3 coordiantes from adaptableUbiAtoms as the x, y, z coordinates havent been changed yet
    x_HZ3, y_HZ3, z_HZ3 = conjugatingAtoms.at[index_LYS_HZ3, 'x'], conjugatingAtoms.at[index_LYS_HZ3, 'y'], conjugatingAtoms.at[index_LYS_HZ3, 'z']

    ### get the coordinates of C on glycine on the rotated carbon... working ubi...
    x_C, y_C, z_C = newAtoms.at[index_PG_C, 'x3'], newAtoms.at[index_PG_C, 'y3'], newAtoms.at[index_PG_C, 'z3']

    ### get the transformation matrix.. 
    x_transform, y_transform, z_transform = x_HZ3 - x_C, y_HZ3 - y_C, z_HZ3 - z_C

    ### apply the transformation matrix... 
    newAtoms['x4'] = newAtoms['x3']+x_transform
    newAtoms['y4'] = newAtoms['y3']+y_transform
    newAtoms['z4'] = newAtoms['z3']+z_transform

    ### make x4, y4, z4 the coordinates of the adaptable protein...
    newAtoms['x'] = newAtoms['x4']
    newAtoms['y'] = newAtoms['y4']
    newAtoms['z'] = newAtoms['z4']

    ## drop everything...
    newAtoms = newAtoms.drop(['x1', 'x2', 'x3', 'x4', 'y1', 'y2', 'y3', 'y4', 'z1', 'z2', 'z3', 'z4', 'np_rotated', 'nparray'], axis=1)
    conjugatingAtoms = conjugatingAtoms.drop(['x1', 'y1', 'z1', 'nparray'], axis=1)
    return newAtoms, full_conjugatingAtoms, newBonds, full_conjugatingUbiBonds

## conjugated lysine number... using K63 for now... 
## ATOMIC CHANGES...
## NZ 
# 1006 - NZ === this will turn into a N.am
# 1014 - HZ3 === this disappears

def numbering_for_bonds_and_atoms_PG(newAtoms, conjugatingAtoms, newBonds, conjugatingUbiBonds, ubiquitinConjugationSite, protecting_group):
    
    # lysine .. NZ === atom_type will turn into a N.am
    ## this can be any of the lysines... 48 or 63 
    workingLysine = conjugatingAtoms[conjugatingAtoms['subst_id']==ubiquitinConjugationSite]
    # lysine .. NZ === atom_type will turn into a N.am
    ## find the index of the HZ3... index to be deleted
    indexNZ = workingLysine[workingLysine['atom_name'] == 'NZ'].index[0]
    atomID_NZ = conjugatingAtoms.at[indexNZ, 'atom_id']
    # lysine .. HZ3 === this disappears
    ## find the index of the HZ3... index to be deleted
    indexHZ3 = workingLysine[workingLysine['atom_name'] == 'HZ3'].index[0]
    atomID_HZ3 = conjugatingAtoms.at[indexHZ3, 'atom_id']

    # 1006 - NZ === this will turn into a N.am
    # 1014 - HZ3 === this disappears
    conjugatingAtoms.at[indexNZ, 'atom_type'] = 'N.am'
    conjugatingAtoms = conjugatingAtoms.drop(index=indexHZ3)


    ## adaptableUbiAtoms is here the new ubi...
    ## bonding changes.. 
    ## change..adaptable bonds... & atoms... add 2000.. this could be anything... 
    ## atom id changes...
    newAtoms['atom_id'] = newAtoms['atom_id'] + 1000000
    ## origin_atom_id and target_atom_id change
    newBonds['origin_atom_id'] = newBonds['origin_atom_id'] + 1000000
    newBonds['target_atom_id'] = newBonds['target_atom_id'] + 1000000
    newBonds['bond_id'] = newBonds['bond_id'] + 1000000


    ## conjugated glycine number... always using G76 for now... 
    ## ATOMIC CHANGES...
    # C
    # 1226 - C === nothing changes
    # 1228 - OXT === 

    # glycine .. OXT === this disappears
    ## this can be any of the lysines...
    ## find the index of the OXT... index to be deleted

    # indexC = newAtoms[newAtoms['atom_name'] == 'C'].index[0]
    # atomID_C = newUbiAtoms.at[indexC, 'atom_id']
    # indexOXT = workingCterminusGlycine[workingCterminusGlycine['atom_name'] == 'OXT'].index[0]
    # atomID_OXT = newUbiAtoms.at[indexOXT, 'atom_id']

    index_PG_Cl= newAtoms[newAtoms['atom_type'] == 'LG'].index[0]
    atomID_PG_Cl = newAtoms.at[index_PG_Cl, 'atom_id']
    # check the bonding... 
    Cl_C_index = newBonds[(newBonds['origin_atom_id'] == atomID_PG_Cl)|(newBonds['target_atom_id']== atomID_PG_Cl)].index[0]
    # atomID_PG_C equals the atom id from the value that is not Cl
    if newBonds.at[Cl_C_index, 'origin_atom_id'] == atomID_PG_Cl: 
        atomID_PG_C = newBonds.at[Cl_C_index, 'target_atom_id']
    else: 
        atomID_PG_C = newBonds.at[Cl_C_index, 'origin_atom_id']
    index_PG_C= newAtoms[newAtoms['atom_id'] == atomID_PG_C].index[0]

    # 1226 - C === nothing changes
    # 1228 - OXT === 
    ## drop ATOMS
    newAtoms= newAtoms.drop(index=index_PG_Cl)

    ## BONDING CHANGES...
    # on the lysine loss of the NZ - HZ3 bond ==== startUbiBonds
    NZ_HZ3_index = conjugatingUbiBonds[(conjugatingUbiBonds['origin_atom_id'] == atomID_NZ) & (conjugatingUbiBonds['target_atom_id']== atomID_HZ3)].index[0]
    # on the glycine loss of the C - OXT bond ==== adaptableUbiBonds
    #C_OXT_index = newUbiBonds[(newUbiBonds['origin_atom_id'] == atomID_C) & (newUbiBonds['target_atom_id']== atomID_OXT)].index[0]

    ## drop the BONDS only the C_OXT BOND AS THE NZ_HZ3 get converted below
    newBonds= newBonds.drop(index=Cl_C_index) 

    ## addition of all the adaptable bonds... see what another peptide bond does...
    ## you can actually just replace the C_OXT bond... replace the index...
    ## its complicated because you need to think of this backwards... so you need to replace the NZ_HZ3_index
    # converting the NZ_HZ3 into a NZ_C bond...
    conjugatingUbiBonds.loc[NZ_HZ3_index,'target_atom_id'] = atomID_PG_C
    conjugatingUbiBonds.loc[NZ_HZ3_index,'bond_type'] = 'am'

    conjugatingUbiBonds.loc[NZ_HZ3_index]
    ## add the C NZ bond...
    ## last bit will be the library of number changes as you add a ubiquitin... 

    ## reorder numbers in subst_id & subst_name
    # adaptableAA = newAtoms['subst_name'].str[:3]
    ## this number is MESSY figure it out
    if protecting_group == 'SMAC':
        adaptableAA = "SMA"
    elif protecting_group == 'ABOC':
        adaptableAA = "ABO"

    ## ALERTTTTT here the value will need to be consdiered this is tough...
    ## add list that knows the number of amino acids + side chains 

    ### CAREFUL THIS IS DIFFERENT FOR PROTEIN VS PROTECTING GROUP...
    ### protecting group you are numbering based on the amino acid position... 
    ### protein you are adding to the number of amino acids..
    newAtoms['subst_id'] = ubiquitinConjugationSite +1000000
    newAtoms.loc[:,'subst_id']= newAtoms.loc[:,'subst_id'].astype(str)
    ## protecting groups have their unique identifier + the amino acid number they are attached to
    #newAtoms['subst_name'] = adaptableAA + newAtoms['subst_id']
    newAtoms['subst_name'] = adaptableAA + '_' + str(ubiquitinConjugationSite)
    newAtoms.loc[:,'subst_id']= newAtoms.loc[:,'subst_id'].astype(int)

    ## reorder numbers.. and  create a new dictionary
    ## add a million to bond id... 
    ## reorder bond_id and atom_id and change... atom_id, origin_atom_id, and taarget_atom_id
    ## then map the values.. for bonds and atoms... 
    # https://sparkbyexamples.com/pandas/pandas-remap-values-in-column-with-a-dictionary-dict/

    ## atoms && bonds
    concatUbiAtoms = pd.concat([conjugatingAtoms,newAtoms])
    concatUbiAtoms['new_atom_id'] = range(1,len(concatUbiAtoms)+1)
    atomID_Dict = pd.Series(concatUbiAtoms.new_atom_id.values,index=concatUbiAtoms.atom_id).to_dict()

    #Using exhaustive map() dict
    concatUbiAtoms["atom_id"]=concatUbiAtoms['atom_id'].map(atomID_Dict)
    concatUbiBonds = pd.concat([conjugatingUbiBonds,newBonds])
    concatUbiBonds['bond_id'] = range(1,len(concatUbiBonds)+1)
    concatUbiBonds["origin_atom_id"]=concatUbiBonds['origin_atom_id'].map(atomID_Dict)
    concatUbiBonds["target_atom_id"]=concatUbiBonds['target_atom_id'].map(atomID_Dict)
    concatUbiAtoms=concatUbiAtoms.reset_index()
    concatUbiBonds=concatUbiBonds.reset_index()
    concatUbiAtoms=concatUbiAtoms.drop('index', axis=1)
    concatUbiAtoms=concatUbiAtoms.drop('new_atom_id', axis=1)
    concatUbiBonds=concatUbiBonds.drop('index', axis=1)

    concatUbiAtoms.loc[:,'x'] = concatUbiAtoms.loc[:,'x'].astype(float)
    concatUbiAtoms.loc[:,'y'] = concatUbiAtoms.loc[:,'y'].astype(float)
    concatUbiAtoms.loc[:,'z'] = concatUbiAtoms.loc[:,'z'].astype(float)
    
    concatUbiAtoms['x'] = concatUbiAtoms['x'].apply(lambda x: format(x, '.4f'))
    concatUbiAtoms['y'] = concatUbiAtoms['y'].apply(lambda x: format(x, '.4f'))
    concatUbiAtoms['z'] = concatUbiAtoms['z'].apply(lambda x: format(x, '.4f'))
    return concatUbiAtoms, concatUbiBonds
    
### CHANGES NEED TO BE MADE... the the rotation should be taken from the start ubiquitin
## as in the one which will be added tooo...

### this should all go before the bonding and atom changes...
## as bonding and atom changes delete atoms..


def rotation_transformation_of_new_ubiquitin(ubiquitinAllAtoms, ubiquitinAllBonds, single_ubiAllAtoms, single_ubiAllBonds, ubiquitinConjugationSite, conjugatingAtomsStart, conjugatingAtomsEnd):
    # define the new ubiquitin... and find the lysine of choice... 
    ## copy so that the value isnt assigned...
    ## this is be a starting place.. 

    ## adaptable atoms is a copy of the ubiquitin of choice...
    ## this will change to be a subset of the original ubiquitin...
    ## adaptable ubi is the ubi unit being added...
    ## this will not be the case...
    full_conjugatingAtoms = ubiquitinAllAtoms.copy()
    full_conjugatingUbiBonds = ubiquitinAllBonds.copy()

    ## find the lysines of interest     
    conjugatingAtoms = full_conjugatingAtoms[(conjugatingAtomsStart <= ubiquitinAllAtoms['subst_id']) &  (conjugatingAtomsEnd  >= ubiquitinAllAtoms['subst_id'])]
    conjugatingUbiBonds = full_conjugatingUbiBonds

    ## adaptable atoms is a copy of the ubiquitin of choice...
    ## this will change to be a subset of the original ubiquitin...
    ## adaptable ubi is the ubi unit being added...
    newUbiAtoms = single_ubiAllAtoms.copy()
    newUbiBonds = single_ubiAllBonds.copy()

    ## FIRST     origin the centre of mass...
    ## a transform matrix that aligns two vectors... does it do scaling and transformation..??
    # minus the centre of mass from all the positions, 
    # this zeros the protein on the centre of mass and the coordinates become the origin... and the rotation occurs around the origin..
    ## TRANSFORM ALL THE COORDINATES 
    ## end of day... lets get the numbers reordered. done....
    ## start looking at the transformation...  

    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...

    ## ALERTTTT: the centre of mass here should be of the ubiquitin that is going to be added 
    ## when you transfer coordiantes there may be misalignment.. be careful here..
    ## pull rotational and transformational matrices from the coordaintes of the ubi of interest 
    # and apply to the new ubi
    ## you can do this from the origin..
    ## makes things easier... just find the coordinates of interest...
    ## try this with the trimer..


    # ubi of interest... find the aa numbers that matter...
    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...
    conjugatingAtoms['c'] = conjugatingAtoms['atom_name'].str[:1]
    atomID_mass_Dict = {'N':'14', 'C': '12', 'O':'16', 'H':'1', 'S': '32'}
    conjugatingAtoms['c'] = conjugatingAtoms['c'].map(atomID_mass_Dict)
    conjugatingAtoms.loc[:,'c'] = conjugatingAtoms.loc[:,'c'].astype(int)
    conjugatingAtoms.loc[:,'x'] = conjugatingAtoms.loc[:,'x'].astype(float)
    conjugatingAtoms.loc[:,'y'] = conjugatingAtoms.loc[:,'y'].astype(float)
    conjugatingAtoms.loc[:,'z'] = conjugatingAtoms.loc[:,'z'].astype(float)
    ## finding centre of mass
    conjugatingAtoms['cx'] = conjugatingAtoms['c']*conjugatingAtoms['x']
    conjugatingAtoms['cy'] = conjugatingAtoms['c']*conjugatingAtoms['y']
    conjugatingAtoms['cz'] = conjugatingAtoms['c']*conjugatingAtoms['z']
    x_COM=conjugatingAtoms['cx'].sum()/conjugatingAtoms['c'].sum()
    y_COM=conjugatingAtoms['cy'].sum()/conjugatingAtoms['c'].sum()
    z_COM=conjugatingAtoms['cz'].sum()/conjugatingAtoms['c'].sum()
    ## drop the new columns
    conjugatingAtoms=conjugatingAtoms.drop('c', axis=1)
    conjugatingAtoms=conjugatingAtoms.drop('cx', axis=1)
    conjugatingAtoms=conjugatingAtoms.drop('cy', axis=1)
    conjugatingAtoms=conjugatingAtoms.drop('cz', axis=1)

    # new ubiquitin... 
    ## this could also be calculated from the conjugating ubiquitin
    ## there might be a slight difference in centre of mass... as the OXT is missing...
    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...
    newUbiAtoms['c'] = newUbiAtoms['atom_name'].str[:1]
    atomID_mass_Dict = {'N':'14', 'C': '12', 'O':'16', 'H':'1', 'S': '32'}
    newUbiAtoms['c'] = newUbiAtoms['c'].map(atomID_mass_Dict)
    newUbiAtoms.loc[:,'c'] = newUbiAtoms.loc[:,'c'].astype(int)
    newUbiAtoms.loc[:,'x'] = newUbiAtoms.loc[:,'x'].astype(float)
    newUbiAtoms.loc[:,'y'] = newUbiAtoms.loc[:,'y'].astype(float)
    newUbiAtoms.loc[:,'z'] = newUbiAtoms.loc[:,'z'].astype(float)
    
    newUbiAtoms['cx'] = newUbiAtoms['c']*newUbiAtoms['x']
    newUbiAtoms['cy'] = newUbiAtoms['c']*newUbiAtoms['y']
    newUbiAtoms['cz'] = newUbiAtoms['c']*newUbiAtoms['z']
    x_NEW_COM=newUbiAtoms['cx'].sum()/newUbiAtoms['c'].sum()
    y_NEW_COM=newUbiAtoms['cy'].sum()/newUbiAtoms['c'].sum()
    z_NEW_COM=newUbiAtoms['cz'].sum()/newUbiAtoms['c'].sum()
    ## drop the new columns
    newUbiAtoms=newUbiAtoms.drop('c', axis=1)
    newUbiAtoms=newUbiAtoms.drop('cx', axis=1)
    newUbiAtoms=newUbiAtoms.drop('cy', axis=1)
    newUbiAtoms=newUbiAtoms.drop('cz', axis=1)

    ## working glycine belongs to the new ubiquitin... never changes...
    ## need to find centre of mass for this molecule...
    ### find the sc terminal glycine...
    workingCterminus_max = newUbiAtoms[newUbiAtoms['status_bit']=='BACKBONE']['subst_id'].max()
    workingCterminusGlycine = newUbiAtoms[newUbiAtoms['subst_id'] == workingCterminus_max]    ## find the index of the OXT... index to be deleted
    index_GLY_OXT= workingCterminusGlycine[workingCterminusGlycine['atom_name'] == 'OXT'].index[0]
    index_GLY_C= workingCterminusGlycine[workingCterminusGlycine['atom_name'] == 'C'].index[0]

    ## may need to look at the vector of each bond... so the C_OXT 
    ## and the NZ3 -- HZ3

    ### working lysine belongs to the posiition of conjugation...
    ## this can be on one of the ubiquitins of the multimer
    workingLysine = conjugatingAtoms[conjugatingAtoms['subst_id'] == ubiquitinConjugationSite]

    print(conjugatingAtoms)
    print(workingLysine)
    print(ubiquitinConjugationSite)
    index_LYS_NZ = workingLysine[workingLysine['atom_name'] == 'NZ'].index[0]
    index_LYS_HZ3 = workingLysine[workingLysine['atom_name'] == 'HZ3'].index[0]
    index_LYS_C = workingLysine[workingLysine['atom_name'] == 'C'].index[0]

    ## Centre the Coordinates on the centre of mass (COM)...
    ## probably dont need this...
    conjugatingAtoms['x1'] = conjugatingAtoms['x']-x_COM
    conjugatingAtoms['y1'] = conjugatingAtoms['y']-y_COM
    conjugatingAtoms['z1'] = conjugatingAtoms['z']-z_COM
    conjugatingAtoms['nparray'] = conjugatingAtoms[['x1', 'y1', 'z1']].values.tolist()
    
    ## Centre the Coordinates on the centre of mass (COM)...
    newUbiAtoms['x1'] = newUbiAtoms['x']-x_NEW_COM
    newUbiAtoms['y1'] = newUbiAtoms['y']-y_NEW_COM
    newUbiAtoms['z1'] = newUbiAtoms['z']-z_NEW_COM
    newUbiAtoms['nparray'] = newUbiAtoms[['x1', 'y1', 'z1']].values.tolist()

    ## find the coordinates of glysine carbon and lysine carbon
    ## need to decide whether this is done by bond vector or done by the centre of mass to central atom...
    x_GLY_OXT,y_GLY_OXT,z_GLY_OXT  = workingCterminusGlycine.at[index_GLY_OXT, 'x'],workingCterminusGlycine.at[index_GLY_OXT, 'y'],workingCterminusGlycine.at[index_GLY_OXT, 'z']
    x_GLY_C,y_GLY_C,z_GLY_C  = workingCterminusGlycine.at[index_GLY_C, 'x'],workingCterminusGlycine.at[index_GLY_C, 'y'],workingCterminusGlycine.at[index_GLY_C, 'z']
    x_LYS_HZ3,y_LYS_HZ3,z_LYS_HZ3  = workingLysine.at[index_LYS_HZ3, 'x'],workingLysine.at[index_LYS_HZ3, 'y'],workingLysine.at[index_LYS_HZ3, 'z']
    x_LYS_NZ,y_LYS_NZ,z_LYS_NZ  = workingLysine.at[index_LYS_NZ, 'x'],workingLysine.at[index_LYS_NZ, 'y'],workingLysine.at[index_LYS_NZ, 'z']
    x_LYS_C,y_LYS_C,z_LYS_C  = workingLysine.at[index_LYS_C, 'x'],workingLysine.at[index_LYS_C, 'y'],workingLysine.at[index_LYS_C, 'z']

    ## SECOND    find your rotation matrix.. 
    ## maybe need to minus the points here... flip the coordinates...
    ### here flip the glysine coordinates as we are looking for a smaller rotation... 
    ## we are rotating around the centree of mass
    ## because we are going to the one 180 position of the final position...
    
    ## two different centre of mass... glycine moves v little... 
    ## these vectors are important!!!! 
    # bond replacement
    #vec1 = np.array([x_GLY_C-x_GLY_OXT, y_GLY_C-y_GLY_OXT, z_GLY_C-z_GLY_OXT])
    # carbon to centre of mass replacement
    vec1 = np.array([x_NEW_COM-x_GLY_C, y_NEW_COM-y_GLY_C, z_NEW_COM-z_GLY_C])

    ## lysine moves much more...
    # bond replacement
    #vec2 = np.array([x_LYS_HZ3-x_LYS_NZ, y_LYS_HZ3-y_LYS_NZ, z_LYS_HZ3-z_LYS_NZ])
    # carbon to centre of mass replacement
    vec2 = np.array([x_LYS_C-x_COM, y_LYS_C-y_COM, z_LYS_C-z_COM])


    ## find the rotation matrix...
    rotationVector = rotation_matrix_from_vectors(vec1, vec2)
    r = R.from_matrix(rotationVector)
    r.as_matrix().shape
    recreatedVector = r.apply(vec1)

    ## check that the rotation occured around the centre of mass... 
    ## a scaling matrix should equalise the two vectors.. vec2 and vec4 here 
    ## so this works beautifully....
    vec = recreatedVector
    testvec1 = vec[0]**2 + vec[1]**2 + vec[2]**2
    testvec2 = vec2[0]**2 + vec2[1]**2 + vec2[2]**2
    ratio = math.sqrt(testvec1)/math.sqrt(testvec2)
    recreatedVector.round(decimals=2) == (vec2*ratio).round(decimals=2)



    ## THIRD    ROTATIONAL TRANSFORMATION OF ALL POINTS on the new ubiquitin...
    newUbiAtoms['np_rotated'] = newUbiAtoms['nparray'].map(lambda x: r.apply(np.array(x)))
    newUbiAtoms['x2'] = newUbiAtoms['np_rotated'].map(lambda x: x[0])
    newUbiAtoms['y2'] = newUbiAtoms['np_rotated'].map(lambda x: x[1])
    newUbiAtoms['z2'] = newUbiAtoms['np_rotated'].map(lambda x: x[2])


    ## FOURTH    return the points to the centre of mass of the ubiquitin of interest...
    ## apply centre of mass transformation from ubi of interest to the new ubi
    newUbiAtoms['x3'] = newUbiAtoms['x2']+x_COM
    newUbiAtoms['y3'] = newUbiAtoms['y2']+y_COM
    newUbiAtoms['z3'] = newUbiAtoms['z2']+z_COM

    
    
    ### FIFTH    move the glysine on C of new ubiquitin to the HZ3 position...

    ## and then push all the points to where the HZ3 was... from the new coordinate...
    ## FIFTH    TRANSFORM ALL THE COORDINATES to the HZ3 from the C OF THE ROTATED GLYCINE...
    ### get the coordinates of HZ3 on the lysine of interest
    # use adapted as the HZ3 is still there... adapted HZ3 to startubi C...

    # lysine .. HZ3 === this disappears
    ## find the index of the HZ3 on the lysine of the starting ubiquitin... index to be deleted
    ## can still pul HZ3 coordiantes from adaptableUbiAtoms as the x, y, z coordinates havent been changed yet
    x_HZ3, y_HZ3, z_HZ3 = conjugatingAtoms.at[index_LYS_HZ3, 'x'], conjugatingAtoms.at[index_LYS_HZ3, 'y'], conjugatingAtoms.at[index_LYS_HZ3, 'z']

    ### get the coordinates of C on glycine on the rotated carbon... working ubi...
    x_C, y_C, z_C = newUbiAtoms.at[index_GLY_C, 'x3'], newUbiAtoms.at[index_GLY_C, 'y3'], newUbiAtoms.at[index_GLY_C, 'z3']

    ### get the transformation matrix.. 
    x_transform, y_transform, z_transform = x_HZ3 - x_C, y_HZ3 - y_C, z_HZ3 - z_C

    ### apply the transformation matrix... 
    newUbiAtoms['x4'] = newUbiAtoms['x3']+x_transform
    newUbiAtoms['y4'] = newUbiAtoms['y3']+y_transform
    newUbiAtoms['z4'] = newUbiAtoms['z3']+z_transform

    ### make x4, y4, z4 the coordinates of the adaptable protein...
    newUbiAtoms['x'] = newUbiAtoms['x4']
    newUbiAtoms['y'] = newUbiAtoms['y4']
    newUbiAtoms['z'] = newUbiAtoms['z4']

    ## drop everything...
    newUbiAtoms = newUbiAtoms.drop(['x1', 'x2', 'x3', 'x4', 'y1', 'y2', 'y3', 'y4', 'z1', 'z2', 'z3', 'z4', 'np_rotated', 'nparray'], axis=1)
    conjugatingAtoms = conjugatingAtoms.drop(['x1', 'y1', 'z1', 'nparray'], axis=1)
    return newUbiAtoms, full_conjugatingAtoms, newUbiBonds, full_conjugatingUbiBonds




## conjugated lysine number... using K63 for now... 
## ATOMIC CHANGES...
## NZ 
# 1006 - NZ === this will turn into a N.am
# 1014 - HZ3 === this disappears

def numbering_for_bonds_and_atoms(newUbiAtoms, conjugatingAtoms, newUbiBonds, conjugatingUbiBonds, ubiquitinConjugationSite, lastBackBoneSubstID):
    
    # lysine .. NZ === atom_type will turn into a N.am
    ## this can be any of the lysines... 48 or 63 
    workingLysine = conjugatingAtoms[conjugatingAtoms['subst_id']==ubiquitinConjugationSite]
    # lysine .. NZ === atom_type will turn into a N.am
    ## find the index of the HZ3... index to be deleted
    indexNZ = workingLysine[workingLysine['atom_name'] == 'NZ'].index[0]
    atomID_NZ = conjugatingAtoms.at[indexNZ, 'atom_id']
    # lysine .. HZ3 === this disappears
    ## find the index of the HZ3... index to be deleted
    indexHZ3 = workingLysine[workingLysine['atom_name'] == 'HZ3'].index[0]
    atomID_HZ3 = conjugatingAtoms.at[indexHZ3, 'atom_id']

    # 1006 - NZ === this will turn into a N.am
    # 1014 - HZ3 === this disappears
    conjugatingAtoms.at[indexNZ, 'atom_type'] = 'N.am'
    conjugatingAtoms = conjugatingAtoms.drop(index=indexHZ3)


    ## adaptableUbiAtoms is here the new ubi...
    ## bonding changes.. 
    ## change..adaptable bonds... & atoms... add 2000.. this could be anything... 
    ## atom id changes...
    newUbiAtoms['atom_id'] = newUbiAtoms['atom_id'] + 1000000
    ## origin_atom_id and target_atom_id change
    newUbiBonds['origin_atom_id'] = newUbiBonds['origin_atom_id'] + 1000000
    newUbiBonds['target_atom_id'] = newUbiBonds['target_atom_id'] + 1000000
    newUbiBonds['bond_id'] = newUbiBonds['bond_id'] + 1000000


    ## conjugated glycine number... always using G76 for now... 
    ## ATOMIC CHANGES...
    # C
    # 1226 - C === nothing changes
    # 1228 - OXT === 

    # glycine .. OXT === this disappears
    ## this can be any of the lysines...
    ## currently always 76...
    workingCterminus_max = newUbiAtoms[newUbiAtoms['status_bit']=='BACKBONE']['subst_id'].max()
    workingCterminusGlycine = newUbiAtoms[newUbiAtoms['subst_id'] == workingCterminus_max]    ## find the index of the OXT... index to be deleted
    ## find the index of the OXT... index to be deleted
    indexC = workingCterminusGlycine[workingCterminusGlycine['atom_name'] == 'C'].index[0]
    atomID_C = newUbiAtoms.at[indexC, 'atom_id']
    indexOXT = workingCterminusGlycine[workingCterminusGlycine['atom_name'] == 'OXT'].index[0]
    atomID_OXT = newUbiAtoms.at[indexOXT, 'atom_id']

    # 1226 - C === nothing changes
    # 1228 - OXT === 
    ## drop ATOMS
    newUbiAtoms= newUbiAtoms.drop(index=indexOXT)

    ## BONDING CHANGES...
    # on the lysine loss of the NZ - HZ3 bond ==== startUbiBonds
    NZ_HZ3_index = conjugatingUbiBonds[(conjugatingUbiBonds['origin_atom_id'] == atomID_NZ) & (conjugatingUbiBonds['target_atom_id']== atomID_HZ3)].index[0]
    # on the glycine loss of the C - OXT bond ==== adaptableUbiBonds
    C_OXT_index = newUbiBonds[(newUbiBonds['origin_atom_id'] == atomID_C) & (newUbiBonds['target_atom_id']== atomID_OXT)].index[0]

    ## drop the BONDS only the C_OXT BOND AS THE NZ_HZ3 get converted below
    newUbiBonds= newUbiBonds.drop(index=C_OXT_index) 

    ## addition of all the adaptable bonds... see what another peptide bond does...
    ## you can actually just replace the C_OXT bond... replace the index...
    ## its complicated because you need to think of this backwards... so you need to replace the NZ_HZ3_index
    # converting the NZ_HZ3 into a NZ_C bond...
    conjugatingUbiBonds.loc[NZ_HZ3_index,'target_atom_id'] = atomID_C
    conjugatingUbiBonds.loc[NZ_HZ3_index,'bond_type'] = 'am'

    conjugatingUbiBonds.loc[NZ_HZ3_index]
    ## add the C NZ bond...
    ## last bit will be the library of number changes as you add a ubiquitin... 

    ## reorder numbers in subst_id & subst_name
    adaptableAA = newUbiAtoms['subst_name'].str[:3]
    ## ALERTTTTT here the value will need to be consdiered this is tough...
    newUbiAtoms['subst_id'] = newUbiAtoms['subst_id']+lastBackBoneSubstID
    newUbiAtoms.loc[:,'subst_id']= newUbiAtoms.loc[:,'subst_id'].astype(str)
    newUbiAtoms['subst_name'] = adaptableAA + newUbiAtoms['subst_id']
    newUbiAtoms.loc[:,'subst_id']= newUbiAtoms.loc[:,'subst_id'].astype(int)

    ## reorder numbers.. and  create a new dictionary
    ## add a million to bond id... 
    ## reorder bond_id and atom_id and change... atom_id, origin_atom_id, and taarget_atom_id
    ## then map the values.. for bonds and atoms... 
    # https://sparkbyexamples.com/pandas/pandas-remap-values-in-column-with-a-dictionary-dict/

    ## atoms && bonds
    concatUbiAtoms = pd.concat([conjugatingAtoms,newUbiAtoms])
    concatUbiAtoms['new_atom_id'] = range(1,len(concatUbiAtoms)+1)
    atomID_Dict = pd.Series(concatUbiAtoms.new_atom_id.values,index=concatUbiAtoms.atom_id).to_dict()

    #Using exhaustive map() dict
    concatUbiAtoms["atom_id"]=concatUbiAtoms['atom_id'].map(atomID_Dict)
    concatUbiBonds = pd.concat([conjugatingUbiBonds,newUbiBonds])
    concatUbiBonds['bond_id'] = range(1,len(concatUbiBonds)+1)
    concatUbiBonds["origin_atom_id"]=concatUbiBonds['origin_atom_id'].map(atomID_Dict)
    concatUbiBonds["target_atom_id"]=concatUbiBonds['target_atom_id'].map(atomID_Dict)
    concatUbiAtoms=concatUbiAtoms.reset_index()
    concatUbiBonds=concatUbiBonds.reset_index()
    concatUbiAtoms=concatUbiAtoms.drop('index', axis=1)
    concatUbiAtoms=concatUbiAtoms.drop('new_atom_id', axis=1)
    concatUbiBonds=concatUbiBonds.drop('index', axis=1)

    concatUbiAtoms.loc[:,'x'] = concatUbiAtoms.loc[:,'x'].astype(float)
    concatUbiAtoms.loc[:,'y'] = concatUbiAtoms.loc[:,'y'].astype(float)
    concatUbiAtoms.loc[:,'z'] = concatUbiAtoms.loc[:,'z'].astype(float)
    
    concatUbiAtoms['x'] = concatUbiAtoms['x'].apply(lambda x: format(x, '.4f'))
    concatUbiAtoms['y'] = concatUbiAtoms['y'].apply(lambda x: format(x, '.4f'))
    concatUbiAtoms['z'] = concatUbiAtoms['z'].apply(lambda x: format(x, '.4f'))
    return concatUbiAtoms, concatUbiBonds


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

def write_mol2(final_ubiAllAtoms,final_ubiAllBonds, file_name):
    maxAtomID = max(final_ubiAllAtoms['atom_id'])
    maxBondID = max(final_ubiAllBonds['bond_id'])
    atoms_combined = push_togther_atoms(final_ubiAllAtoms)
    bonds_combined = push_togther_bonds(final_ubiAllBonds)

    mol2Header = pd.DataFrame({'C': ['#',
                                    '# Created by the Jeffrey Bode Group of ETH Zurich',
                                    '# Software Developed by Eric Kummelstedt',
                                    '# In collaboration with Leo Seild, Sohei Majima and Toshiki Mikami' 
                                    '#', 
                                    '#', 
                                    '@<TRIPOS>MOLECULE', 
                                    'ubi_multimer.pdb', 
                                    str(maxAtomID) + ' ' + str(maxBondID),
                                    'SMALL',
                                    'NO_CHARGES',
                                    '@<TRIPOS>ATOM'
                                    ]})

    pandas_for_mol = pd.concat([mol2Header, atoms_combined, pd.DataFrame({'C': ['@<TRIPOS>BOND']}), bonds_combined, pd.DataFrame({'C': ['@<TRIPOS>SUBSTRUCTURE', '0 pdb 0 ****']}                                                                                                                                                                                                                                                                                             )],ignore_index=True)
    np.savetxt(str(project_root) + "/data/jsons_to_mol2/" + str(file_name) + ".txt", pandas_for_mol, delimiter="", fmt="%s") 
    os.rename(str(project_root) + "/data/jsons_to_mol2/" + str(file_name) + ".txt", str(project_root) + "/data/jsons_to_mol2/" + str(file_name) + ".mol2")
    
    return str(project_root) + "/data/jsons_to_mol2/" + str(file_name) + ".mol2"
    ## add the number of bonds and number of atoms... 3rd line of the mol2 file...
    ### max of each value....  


def json_to_mol(parent_dictionary):
    """
    Relabel ubiquitin chain numbers, process protein branching and validate input dictionary keys.
    Generate a multimer string name based on a parent dictionary that represents proteins and ubiquitin chains.

    Args:
        parent_dictionary (dict or str): Ubiquitin structure as a dictionary or JSON string.

    Returns:
        dict: Updated polyubiquitin dictionary/JSON with relabeled values.
        context: Upadated dictionary with the information regarding the polyubiquitin. Includes;
            - chain_number_list
            - chain_length_list
            - multimer_string_name
            - max_chain_number
            - ABOC_lysines
            - SMAC_lysines
            - conjugated_lysines

    Raises:
        KeyError: If required keys are missing.
        ValueError: If the input is not a dictionary or valid JSON string.

    """
    
    # Create into separate function
    UBIQUITIN_DF = pd.read_csv(str(project_root) + "/data/mol2_database/ubiquitin.txt", names=['#'])
    HISTAG_UBIQUITIN_DF = pd.read_csv(str(project_root) + "/data/mol2_database/histag-ubiquitin.txt", names=['#'])
    ABOC_DF = pd.read_csv(str(project_root) + "/data/mol2_database/ABOC_LG.txt", names=['#'])
    SMAC_DF = pd.read_csv(str(project_root) + "/data/mol2_database/SMAC_LG.txt", names=['#'])

    ATOMS_UBIQUITIN_DF, BONDS_UBIQUITIN_DF = get_BONDS_ATOMS(UBIQUITIN_DF)
    ATOMS_HISTAG_UBIQUITIN_DF, BONDS_HISTAG_UBIQUITIN_DF = get_BONDS_ATOMS(HISTAG_UBIQUITIN_DF)
    ATOMS_ABOC_DF, BONDS_ABOC_DF = get_BONDS_ATOMS(ABOC_DF)
    ATOMS_SMAC_DF, BONDS_SMAC_DF = get_BONDS_ATOMS(SMAC_DF)

    mol2_dictionary = {
        'ATOMS_UBIQUITIN_DF' : ATOMS_UBIQUITIN_DF,
        'ATOMS_HISTAG_UBIQUITIN_DF' : ATOMS_HISTAG_UBIQUITIN_DF,
        'ATOMS_ABOC_DF' : ATOMS_ABOC_DF,
        'ATOMS_SMAC_DF' : ATOMS_SMAC_DF,
        'BONDS_UBIQUITIN_DF' : BONDS_UBIQUITIN_DF,
        'BONDS_HISTAG_UBIQUITIN_DF' : BONDS_HISTAG_UBIQUITIN_DF,
        'BONDS_ABOC_DF' : BONDS_ABOC_DF,
        'BONDS_SMAC_DF' : BONDS_SMAC_DF
    }

    # Ensure that the parent dictionary is a JSON
    parent_dictionary = convert_json_to_dict(parent_dictionary)

    # Ensure that the parent dictionary has all the valid keys 
    validate_protein_keys(parent_dictionary)
    
    # Initialize context object to track global-like variables
    context = {
        "chain_number_list": [1],
        "chain_length_list": [],
        "multimer_string_name": "",
        "max_chain_number" : "",
        "ABOC_lysines" : [],
        "SMAC_lysines": [],
        "free_lysines": [],
        "conjugated_lysines": []
    }

    adapted_dictionary, context, mol2_dictionary = inner_wrapper_json_to_mol(
        parent_dictionary, context, mol2_dictionary
    )

    file_path = write_mol2(mol2_dictionary['BASE_ATOMS_DF'], mol2_dictionary['BASE_BONDS_DF'], context['multimer_string_name'] )

    return file_path

def inner_wrapper_json_to_mol(input_dictionary, context, mol2_dictionary):
    """
    Recursively process nested dictionaries to relabel and extract ubiquitin chain numbers and process protein branching.
    Recursively process nested dictionaries to construct the multimer string name.
    """
    working_dictionary = copy.deepcopy(input_dictionary)

    # Ensure that the working dictionary is a JSON
    working_dictionary = convert_json_to_dict(working_dictionary)

    # Ensure that the working dictionary has all the valid keys 
    validate_protein_keys(working_dictionary)
    
    ## Process node numbers for the pre-order tree treversal of the protein 
    working_dictionary, context = process_current_protein(working_dictionary, context)

    # Append chain information to multimer string
    # Change to pdbid
    if (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH") & (working_dictionary["chain_number"]==1):
        context["multimer_string_name"] += f"his-GG-{working_dictionary['protein']}-{working_dictionary['chain_number']}-("
    elif (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG") & (working_dictionary["chain_number"]==1):
        context["multimer_string_name"] += f"GG-{working_dictionary['protein']}-{working_dictionary['chain_number']}-("
    else: 
        context["multimer_string_name"] += f"{working_dictionary['protein']}-{working_dictionary['chain_number']}-("

    print(working_dictionary["FASTA_sequence"])

    # TODO: add to function
    # Change to pdbid
    if (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH") & (working_dictionary["chain_number"]==1):
        # TODO 
        mol2_dictionary['BASE_ATOMS_DF'] = mol2_dictionary['ATOMS_HISTAG_UBIQUITIN_DF'].copy()
        mol2_dictionary['BASE_BONDS_DF'] = mol2_dictionary['BONDS_HISTAG_UBIQUITIN_DF'].copy()

    elif (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG") & (working_dictionary["chain_number"]==1):
        # TODO 
        mol2_dictionary['BASE_ATOMS_DF'] = mol2_dictionary['ATOMS_UBIQUITIN_DF'].copy()
        mol2_dictionary['BASE_BONDS_DF'] = mol2_dictionary['BONDS_UBIQUITIN_DF'].copy()
    elif working_dictionary["chain_number"] == 1:
        raise TypeError("The sequence found is not in the database.")
    
    # Log current protein details
    log_protein_details(working_dictionary, context)

    # Ensure that the branching has all the valid keys 
    validate_branching_sites(working_dictionary)

    ## get the branching sites
    working_branching_sites = working_dictionary.get("branching_sites", [])

    ## loop through the branches
    for branch in working_branching_sites:
        # log new branching deatils
        log_branching_details(branch, working_dictionary, context)
        
        ## process branch
        branch, working_dictionary, context, mol2_dictionary = process_branch_json_to_mol(branch, working_dictionary, context, mol2_dictionary)

        # log end of branching details
        log_end_of_branching()

    log_end_of_protein(working_dictionary)

    # End of multimer string editing 
    context["multimer_string_name"] += ")"

    # Find the maximum chain number
    context = add_max_chain_number(context)
    
    return working_dictionary, context, mol2_dictionary

def process_branch_json_to_mol(branch, working_dictionary, context, mol2_dictionary):
    """
    Handle the logic for each branch site in the protein structure.
    - Must deal with the following;  
    """
    chain_length_list = context["chain_length_list"]
    chain_number_list = context["chain_number_list"]
    chain_number_index = working_dictionary["chain_number"] - 1

    # Find branching site
    branching_number = find_branching_site(branch["sequence_id"], working_dictionary["FASTA_sequence"])
    ubiquitin_conjugation_site = int(sum(chain_length_list[:chain_number_index])) + branching_number
    conjugating_atoms_start = int(sum(chain_length_list[:chain_number_index])) + 1
    conjugating_atoms_end = int(sum(chain_length_list[:(chain_number_index+1)]))
    last_backbone_subst_id = int(sum(chain_length_list))

    # Handle protecting groups
    if branch["children"] in ["SMAC"]:
        ## add protecting group
        context["multimer_string_name"] += f"<{branch['site_name']}_SMAC>"
        context["SMAC_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

        #TODO insert function to add SMAC, with ubiquitin_conjugation_site
        protecting_group = 'SMAC'
        new_atoms, new_bonds = mol2_dictionary['ATOMS_SMAC_DF'].copy(), mol2_dictionary['BONDS_SMAC_DF'].copy()
        base_atoms, base_bonds = mol2_dictionary['BASE_ATOMS_DF'].copy(), mol2_dictionary['BASE_BONDS_DF'].copy()
        ## rotate new ubiquitin molecule
        new_atoms, conjugating_atoms, new_bonds, conjugating_bonds = rotation_transformation_of_new_PG(base_atoms, base_bonds, new_atoms, new_bonds, ubiquitin_conjugation_site, conjugating_atoms_start, conjugating_atoms_end, protecting_group)
        ## renumber the bonding for the ubiquitin molecule
        mol2_dictionary['BASE_ATOMS_DF'], mol2_dictionary['BASE_BONDS_DF'] = numbering_for_bonds_and_atoms_PG(new_atoms, conjugating_atoms, new_bonds, conjugating_bonds, ubiquitin_conjugation_site, protecting_group)

    elif branch["children"] in ["ABOC"]:
        ## add protecting group
        context["multimer_string_name"] += f"<{branch['site_name']}_ABOC>"
        context["ABOC_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

        #TODO insert function to add ABOC, with ubiquitin_conjugation_site
        protecting_group = 'ABOC'
        new_atoms, new_bonds = mol2_dictionary['ATOMS_ABOC_DF'].copy(), mol2_dictionary['BONDS_ABOC_DF'].copy()
        base_atoms, base_bonds = mol2_dictionary['BASE_ATOMS_DF'].copy(), mol2_dictionary['BASE_BONDS_DF'].copy()
        ## rotate new ubiquitin molecule
        new_atoms, conjugating_atoms, new_bonds, conjugating_bonds = rotation_transformation_of_new_PG(base_atoms, base_bonds, new_atoms, new_bonds, ubiquitin_conjugation_site, conjugating_atoms_start, conjugating_atoms_end, protecting_group)
        ## renumber the bonding for the ubiquitin molecule
        mol2_dictionary['BASE_ATOMS_DF'], mol2_dictionary['BASE_BONDS_DF'] = numbering_for_bonds_and_atoms_PG(new_atoms, conjugating_atoms, new_bonds, conjugating_bonds, ubiquitin_conjugation_site, protecting_group)
        
    # this can change later
    # Handle free lysines
    elif (branch['children'] == "") & (branch['site_name'] in ['M1','K6','K11','K27','K29','K33']): 
        print('')

    # Handle K48 & K63 lysines
    elif (branch['children'] == "") & (branch['site_name'] in ['K48', 'K63']): 
        context["free_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

    # Handle branches that have proteins bound  
    elif isinstance(branch["children"], dict):
        context["multimer_string_name"] += f"<{branch['site_name']}_"
        context["conjugated_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

        #TODO insert function to add ABOC, with ubiquitin_conjugation_site
        # TODO: for pdbid
        # Check that 
        if (branch["children"]["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"):
            
            #TODO insert function to add ubiquitin, with ubiquitin_conjugation_site
            protecting_group = 'N/A'
            new_atoms, new_bonds = mol2_dictionary['ATOMS_UBIQUITIN_DF'].copy(), mol2_dictionary['BONDS_UBIQUITIN_DF'].copy()
            base_atoms, base_bonds = mol2_dictionary['BASE_ATOMS_DF'].copy(), mol2_dictionary['BASE_BONDS_DF'].copy()
            ## rotate new ubiquitin molecule
            new_atoms, conjugating_atoms, new_bonds, conjugating_bonds = rotation_transformation_of_new_ubiquitin(base_atoms, base_bonds, new_atoms, new_bonds, ubiquitin_conjugation_site, conjugating_atoms_start, conjugating_atoms_end)
            ## renumber the bonding for the ubiquitin molecule
            mol2_dictionary['BASE_ATOMS_DF'], mol2_dictionary['BASE_BONDS_DF'] = numbering_for_bonds_and_atoms(new_atoms, conjugating_atoms, new_bonds, conjugating_bonds, ubiquitin_conjugation_site, last_backbone_subst_id)

        branch["children"], context, mol2_dictionary = inner_wrapper_json_to_mol(
            branch["children"], context, mol2_dictionary
        )
        context["multimer_string_name"] += ">"

    return branch, working_dictionary, context, mol2_dictionary


