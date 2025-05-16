## O should turn into O.2 not O.3

## use assign_backbone to find the backbone...
## and find the extra O.2's and H's on Nitrogen..

### create your own importer of mol2 files$
import os
import pandas as pd
### this should all go before the bonding and atom changes...
## as bonding and atom changes delete atoms..

from scipy.spatial.transform import Rotation as R
import numpy as np
import math
import json
import copy

### create your own importer of mol2 files$
import os
import pandas as pd
### this should all go before the bonding and atom changes...
## as bonding and atom changes delete atoms..

from scipy.spatial.transform import Rotation as R
import numpy as np
import math
import json
import copy
import sys

home_dir = os.path.expanduser('~')
local_path = home_dir + '/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/'

sys.path.insert(0, local_path + 'aa_functions_folder')
sys.path.insert(0, local_path + '.data')

from labeling_backbones import *

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
    


## change input to dictionary for: bonding_base_atom_id, lost_base_atom_id, bonding_new_atom_id, lost_new_atom_id
## add in previous backbone atom_id
## add in other end backbone atom_id



def rotation_transformation_mol2(base_atoms, base_bonds, new_atoms, new_bonds, rot_transform_dict):
    ## find the two bonds that are being changed and then rotate everything and change everything..
    ## that way you modularize everything...
    ## for larger reactions you can remove certain groups...
    ## and then do a bond replacement..
    bonding_base_atom_id = rot_transform_dict['bonding_base_atom_id']
    lost_base_atom_id = rot_transform_dict['lost_base_atom_id']
    previous_base_atom_id = rot_transform_dict['previous_base_atom_id']
    previous2_base_atom_id = rot_transform_dict['previous2_base_atom_id']
    bonding_new_atom_id = rot_transform_dict['bonding_new_atom_id']
    lost_new_atom_id = rot_transform_dict['lost_new_atom_id']
    previous_new_atom_id = rot_transform_dict['previous_new_atom_id']
    previous2_new_atom_id = rot_transform_dict['previous2_new_atom_id']
    nitrogen_on = rot_transform_dict['nitrogen_on']

    ## add in the other values.... 
    multimer_atoms_start = rot_transform_dict['multimer_atoms_start']
    multimer_atoms_end = rot_transform_dict['multimer_atoms_end']
    ubiquitin_conjugation_site = rot_transform_dict['ubiquitin_conjugation_site']
    which_terminusQ = rot_transform_dict['which_terminusQ']
    base_molecule = rot_transform_dict['base_molecule']
    new_molecule = rot_transform_dict['new_molecule']
    last_backbond_substid = rot_transform_dict['last_backbond_substid']

    ## first step is to find the indexes
    bonding_base_index = base_atoms[base_atoms['atom_id'] == bonding_base_atom_id].index[0]
    lost_base_index = base_atoms[base_atoms['atom_id'] == lost_base_atom_id].index[0]
    previous_base_index = base_atoms[base_atoms['atom_id'] == previous_base_atom_id].index[0]
    previous2_base_index = base_atoms[base_atoms['atom_id'] == previous2_base_atom_id].index[0]
    
    bonding_new_index = new_atoms[new_atoms['atom_id'] == bonding_new_atom_id].index[0]
    lost_new_index = new_atoms[new_atoms['atom_id'] == lost_new_atom_id].index[0]
    previous_new_index = new_atoms[new_atoms['atom_id'] == previous_new_atom_id].index[0]
    previous2_new_index = new_atoms[new_atoms['atom_id'] == previous2_new_atom_id].index[0]

    conjugating_atoms = base_atoms[(multimer_atoms_start <= base_atoms['subst_id']) &  (multimer_atoms_end  >= base_atoms['subst_id'])]

    base_atoms['x'] = base_atoms['x'].astype(float)
    base_atoms['y'] = base_atoms['y'].astype(float)
    base_atoms['z'] = base_atoms['z'].astype(float)
    
    # ubi of interest... find the aa numbers that matter...
    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...
    conjugating_atoms['c'] = conjugating_atoms['atom_name'].str[:1]
    atomID_mass_Dict = {'N':'14', 'C': '12', 'O':'16', 'H':'1', 'S': '32'}
    conjugating_atoms['c'] = conjugating_atoms['c'].map(atomID_mass_Dict)
    conjugating_atoms['c'] = conjugating_atoms['c'].astype(int)
    conjugating_atoms['x'] = conjugating_atoms['x'].astype(float)
    conjugating_atoms['y'] = conjugating_atoms['y'].astype(float)
    conjugating_atoms['z'] = conjugating_atoms['z'].astype(float)
    ## finding centre of mass
    conjugating_atoms['cx'] = conjugating_atoms['c']*conjugating_atoms['x']
    conjugating_atoms['cy'] = conjugating_atoms['c']*conjugating_atoms['y']
    conjugating_atoms['cz'] = conjugating_atoms['c']*conjugating_atoms['z']
   
    x_COM=conjugating_atoms['cx'].sum()/conjugating_atoms['c'].sum()
    y_COM=conjugating_atoms['cy'].sum()/conjugating_atoms['c'].sum()
    z_COM=conjugating_atoms['cz'].sum()/conjugating_atoms['c'].sum()
    
    ## this might not be needed.. 
    ## drop the new columns
    conjugating_atoms=conjugating_atoms.drop('c', axis=1)
    conjugating_atoms=conjugating_atoms.drop('cx', axis=1)
    conjugating_atoms=conjugating_atoms.drop('cy', axis=1)
    conjugating_atoms=conjugating_atoms.drop('cz', axis=1)

    ## Centre the Coordinates on the centre of mass (COM)...
    ## probably dont need this for base atoms... 
    conjugating_atoms['x1'] = conjugating_atoms['x']-x_COM
    conjugating_atoms['y1'] = conjugating_atoms['y']-y_COM
    conjugating_atoms['z1'] = conjugating_atoms['z']-z_COM
    conjugating_atoms['nparray'] = conjugating_atoms[['x1', 'y1', 'z1']].values.tolist()

    # new ubiquitin... 
    ## this could also be calculated from the conjugating ubiquitin
    ## there might be a slight difference in centre of mass... as the OXT is missing...
    # find the centre of mass --- use to determine the vector... you can do this yourself...
    # just use centre of mass at the moment... good enough...
    new_atoms['c'] = new_atoms['atom_name'].str[:1]
    atomID_mass_Dict = {'N':'14', 'C': '12', 'O':'16', 'H':'1', 'S': '32', 'L': '17'}
    new_atoms['c'] = new_atoms['c'].map(atomID_mass_Dict)
    new_atoms['c'] = new_atoms['c'].astype(int)
    new_atoms['x'] = new_atoms['x'].astype(float)
    new_atoms['y'] = new_atoms['y'].astype(float)
    new_atoms['z'] = new_atoms['z'].astype(float)
    
    new_atoms['cx'] = new_atoms['c']*new_atoms['x']
    new_atoms['cy'] = new_atoms['c']*new_atoms['y']
    new_atoms['cz'] = new_atoms['c']*new_atoms['z']
    x_NEW_COM=new_atoms['cx'].sum()/new_atoms['c'].sum()
    y_NEW_COM=new_atoms['cy'].sum()/new_atoms['c'].sum()
    z_NEW_COM=new_atoms['cz'].sum()/new_atoms['c'].sum()
    ## drop the new columns
    new_atoms=new_atoms.drop('c', axis=1)
    new_atoms=new_atoms.drop('cx', axis=1)
    new_atoms=new_atoms.drop('cy', axis=1)
    new_atoms=new_atoms.drop('cz', axis=1)

    ## Centre the Coordinates on the centre of mass (COM)...
    new_atoms['x1'] = new_atoms['x']-x_NEW_COM
    new_atoms['y1'] = new_atoms['y']-y_NEW_COM
    new_atoms['z1'] = new_atoms['z']-z_NEW_COM
    new_atoms['nparray'] = new_atoms[['x1', 'y1', 'z1']].values.tolist()
    
    ## find the coordinates of glysine carbon and lysine carbon
    ## need to decide whether this is done by bond vector or done by the centre of mass to central atom...
    x_bonding_base,y_bonding_base,z_bonding_base  = base_atoms.at[bonding_base_index, 'x'],base_atoms.at[bonding_base_index, 'y'],base_atoms.at[bonding_base_index, 'z']
    x_lost_base,y_lost_base,z_lost_base  = base_atoms.at[lost_base_index, 'x'],base_atoms.at[lost_base_index, 'y'],base_atoms.at[lost_base_index, 'z']
    x_previous_base,y_previous_base,z_previous_base  = base_atoms.at[previous_base_index, 'x'],base_atoms.at[previous_base_index, 'y'],base_atoms.at[previous_base_index, 'z']
    x_previous2_base, y_previous2_base, z_previous2_base  = base_atoms.at[previous2_base_index, 'x'],base_atoms.at[previous2_base_index, 'y'],base_atoms.at[previous2_base_index, 'z']

    x_bonding_new,y_bonding_new,z_bonding_new  = new_atoms.at[bonding_new_index, 'x'],new_atoms.at[bonding_new_index, 'y'],new_atoms.at[bonding_new_index, 'z']
    x_lost_new,y_lost_new,z_lost_new  = new_atoms.at[lost_new_index, 'x'],new_atoms.at[lost_new_index, 'y'],new_atoms.at[lost_new_index, 'z']
    x_previous_new,y_previous_new,z_previous_new  = new_atoms.at[previous_new_index, 'x'],new_atoms.at[previous_new_index, 'y'],new_atoms.at[previous_new_index, 'z']
    x_previous2_new,y_previous2_new,z_previous2_new  = new_atoms.at[previous2_new_index, 'x'],new_atoms.at[previous2_new_index, 'y'],new_atoms.at[previous2_new_index, 'z']

    ## SECOND    find your rotation matrix.. 
    ## maybe need to minus the points here... flip the coordinates...
    ### here flip the glysine coordinates as we are looking for a smaller rotation... 
    ## we are rotating around the centree of mass
    ## because we are going to the one 180 position of the final position...


    ## if using the previous2 might not need to use 
    if new_molecule == 'peptide': 
        ## using the backbone bone two back 
        vec1 = np.array([x_previous2_new- x_previous_new, 
                         y_previous2_new- y_previous_new, 
                         z_previous2_new- z_previous_new])
    elif new_molecule == 'protein':
        ## centre of mass change... 
        vec1 = np.array([x_NEW_COM - x_bonding_new, 
                         y_NEW_COM - y_bonding_new, 
                         z_NEW_COM - z_bonding_new])

    if base_molecule == 'peptide': 
        ## using the backbone bone two back 
        vec2 = np.array([x_previous_base-x_previous2_base, 
                         y_previous_base-y_previous2_base, 
                         z_previous_base-z_previous2_base])
    elif base_molecule == 'protein':
        ## centre of mass change... 
        vec2 = np.array([x_bonding_base-x_COM, 
                         y_bonding_base-y_COM, 
                         z_bonding_base-z_COM])

    
    #vec1 = np.array([x_previous2_new- x_previous_new, 
    #                 y_previous2_new- y_previous_new, 
    #                 z_previous2_new- z_previous_new])
    #
    #vec2 = np.array([x_previous_base-x_previous2_base, 
    #                 y_previous_base-y_previous2_base, 
    #                 z_previous_base-z_previous2_base])



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
    new_atoms['np_rotated'] = new_atoms['nparray'].map(lambda x: r.apply(np.array(x)))
    new_atoms['x2'] = new_atoms['np_rotated'].map(lambda x: x[0])
    new_atoms['y2'] = new_atoms['np_rotated'].map(lambda x: x[1])
    new_atoms['z2'] = new_atoms['np_rotated'].map(lambda x: x[2])

    ## FOURTH    return the points to the centre of mass of the ubiquitin of interest...
    ## apply centre of mass transformation from ubi of interest to the new ubi
    new_atoms['x3'] = new_atoms['x2']+x_COM
    new_atoms['y3'] = new_atoms['y2']+y_COM
    new_atoms['z3'] = new_atoms['z2']+z_COM

    ### FIFTH move the glysine on C of new ubiquitin to the HZ3 position...

    ## and then push all the points to where the HZ3 was... from the new coordinate...
    ## FIFTH    TRANSFORM ALL THE COORDINATES to the HZ3 from the C OF THE ROTATED GLYCINE...
    ### get the coordinates of HZ3 on the lysine of interest
    ### 

    ## figure out bonding changes... 
    ## here use the minus of the base previous... 
    x_bonding_base,y_bonding_base,z_bonding_base  = base_atoms.at[bonding_base_index, 'x'],base_atoms.at[bonding_base_index, 'y'],base_atoms.at[bonding_base_index, 'z']
    x_lost_base,y_lost_base,z_lost_base  = base_atoms.at[lost_base_index, 'x'],base_atoms.at[lost_base_index, 'y'],base_atoms.at[lost_base_index, 'z']
    x_previous_base,y_previous_base,z_previous_base  = base_atoms.at[previous_base_index, 'x'],base_atoms.at[previous_base_index, 'y'],base_atoms.at[previous_base_index, 'z']
    x_previous2_base, y_previous2_base, z_previous2_base  = base_atoms.at[previous2_base_index, 'x'],base_atoms.at[previous2_base_index, 'y'],base_atoms.at[previous2_base_index, 'z']

    ## maybe include *1.5...
    x_base_new_position = x_bonding_base + (x_previous_base - x_previous2_base)
    y_base_new_position = y_bonding_base + (y_previous_base - y_previous2_base)
    z_base_new_position = z_bonding_base + (z_previous_base - z_previous2_base)

    # get the coordinates of C on the rotated Nitrogen... working ubi...
    x_bonding_new_rotated, y_bonding_new_rotated, z_bonding_new_rotated = new_atoms.at[bonding_new_index, 'x3'], new_atoms.at[bonding_new_index, 'y3'], new_atoms.at[bonding_new_index, 'z3']
    
    ### get the transformation matrix.. move N to OXT
    x_transform, y_transform, z_transform =  x_base_new_position - x_bonding_new_rotated, y_base_new_position - y_bonding_new_rotated, z_base_new_position - z_bonding_new_rotated

    ### apply the transformation matrix... 
    new_atoms['x4'] = new_atoms['x3']+x_transform
    new_atoms['y4'] = new_atoms['y3']+y_transform
    new_atoms['z4'] = new_atoms['z3']+z_transform

    ### make x4, y4, z4 the coordinates of the adaptable protein...
    new_atoms['x'] = new_atoms['x4']
    new_atoms['y'] = new_atoms['y4']
    new_atoms['z'] = new_atoms['z4']

    ## drop everything...
    new_atoms = new_atoms.drop(['x1', 'x2', 'x3', 'x4', 'y1', 'y2', 'y3', 'y4', 'z1', 'z2', 'z3', 'z4', 'np_rotated', 'nparray'], axis=1)
    # NOT NEEDED base_atoms = base_atoms.drop(['x1', 'y1', 'z1', 'nparray'], axis=1)

    return base_atoms, base_bonds, new_atoms, new_bonds

# nitrogen on new or base
def connectivity_changes_amide_bond(base_atoms, base_bonds, new_atoms, new_bonds, rot_transform_dict):
    
    bonding_base_atom_id = rot_transform_dict['bonding_base_atom_id']
    lost_base_atom_id = rot_transform_dict['lost_base_atom_id']
    previous_base_atom_id = rot_transform_dict['previous_base_atom_id']
    previous2_base_atom_id = rot_transform_dict['previous2_base_atom_id']
    bonding_new_atom_id = rot_transform_dict['bonding_new_atom_id']
    lost_new_atom_id = rot_transform_dict['lost_new_atom_id']
    previous_new_atom_id = rot_transform_dict['previous_new_atom_id']
    previous2_new_atom_id = rot_transform_dict['previous2_new_atom_id']

    ## add in the other values.... 
    multimer_atoms_start = rot_transform_dict['multimer_atoms_start']
    multimer_atoms_end = rot_transform_dict['multimer_atoms_end']
    ubiquitin_conjugation_site = rot_transform_dict['ubiquitin_conjugation_site']
    which_terminusQ = rot_transform_dict['which_terminusQ']
    protecting_group = rot_transform_dict['protecting_group']
    base_molecule = rot_transform_dict['base_molecule']
    new_molecule = rot_transform_dict['new_molecule']
    last_backbond_substid = rot_transform_dict['last_backbond_substid']

    nitrogen_on = rot_transform_dict['nitrogen_on']

    ## atom id changes...
    new_atoms['atom_id'] = new_atoms['atom_id'] + 1000000
    ## origin_atom_id and target_atom_id change
    new_bonds['origin_atom_id'] = new_bonds['origin_atom_id'] + 1000000
    new_bonds['target_atom_id'] = new_bonds['target_atom_id'] + 1000000
    new_bonds['bond_id'] = new_bonds['bond_id'] + 1000000

    lost_new_atom_id = lost_new_atom_id + 1000000
    bonding_new_atom_id = bonding_new_atom_id + 1000000

    ## first step is to find the indexes
    bonding_base_index = base_atoms[base_atoms['atom_id'] == bonding_base_atom_id].index[0]
    lost_base_index = base_atoms[base_atoms['atom_id'] == lost_base_atom_id].index[0]
    bonding_new_index = new_atoms[new_atoms['atom_id'] == bonding_new_atom_id].index[0]
    lost_new_index = new_atoms[new_atoms['atom_id'] == lost_new_atom_id].index[0]
    
    ## change atom_type to 'N.am'
    if nitrogen_on == 'new': 
        new_atoms.at[bonding_new_index, 'atom_type'] = 'N.am'
    elif nitrogen_on == 'base': 
        base_atoms.at[bonding_base_index, 'atom_type'] = 'N.am'

    ## c is the origin atom -- n is the target atom... so the C-O bond turns into a C-N bond and called amide...
    if nitrogen_on == 'new': 
        C_O3_index = base_bonds[((base_bonds['origin_atom_id'] == bonding_base_atom_id) & (base_bonds['target_atom_id'] == lost_base_atom_id))|((base_bonds['origin_atom_id']== lost_base_atom_id) & (base_bonds['target_atom_id'] == bonding_base_atom_id))].index[0]
        new_bonds.loc[C_O3_index,'origin_atom_id'] = bonding_new_atom_id
        new_bonds.loc[C_O3_index,'target_atom_id'] = bonding_base_atom_id
        base_bonds.loc[C_O3_index,'bond_type'] = 'am'
        amide_bond = base_bonds.loc[[C_O3_index]].copy()
        ## for base bonds we need to lose the index....because the bonding new base_if is different
        base_bonds = base_bonds.drop(index= C_O3_index)
        
    elif (nitrogen_on == 'base') & (which_terminusQ == 'protecting_group'): 
        C_LG_index = new_bonds[((new_bonds['origin_atom_id'] == bonding_new_atom_id) & (new_bonds['target_atom_id'] == lost_new_atom_id))|((new_bonds['origin_atom_id']== lost_new_atom_id) & (new_bonds['target_atom_id'] == bonding_new_atom_id))].index[0]
        new_bonds.loc[C_LG_index,'origin_atom_id'] = bonding_new_atom_id
        new_bonds.loc[C_LG_index,'target_atom_id'] = bonding_base_atom_id
        new_bonds.loc[C_LG_index,'bond_type'] = 'am'
        amide_bond = new_bonds.loc[[C_LG_index]].copy()
        new_bonds = new_bonds.drop(index= C_LG_index)

    elif nitrogen_on == 'base': 
        C_O3_index = new_bonds[((new_bonds['origin_atom_id'] == bonding_new_atom_id) & (new_bonds['target_atom_id'] == lost_new_atom_id))|((new_bonds['origin_atom_id']== lost_new_atom_id) & (new_bonds['target_atom_id'] == bonding_new_atom_id))].index[0]
        new_bonds.loc[C_O3_index,'origin_atom_id'] = bonding_new_atom_id
        new_bonds.loc[C_O3_index,'target_atom_id'] = bonding_base_atom_id
        new_bonds.loc[C_O3_index,'bond_type'] = 'am'
        amide_bond = new_bonds.loc[[C_O3_index]].copy()
        new_bonds = new_bonds.drop(index= C_O3_index)
 
    ### for some reason can't find
    new_atoms, new_bonds = remove_atom_from_bonds_atoms(lost_new_atom_id, new_atoms, new_bonds)
    base_atoms, base_bonds = remove_atom_from_bonds_atoms(lost_base_atom_id, base_atoms, base_bonds)

    ## made need some love and care
    if nitrogen_on == 'new': 
        new_bonds = pd.concat([amide_bond, new_bonds])
    elif nitrogen_on == 'base': 
        new_bonds = pd.concat([new_bonds, amide_bond])

    return base_atoms, base_bonds, new_atoms, new_bonds


## do + 1 or -1 and then reorder everything after concatenating the atoms/bonds.. 
## figure out way to work with protecting groups
def new_atom_subst_ids(new_atoms, base_atoms, rot_transform_dict):
    
    ## rot_terminus...
    which_terminusQ = rot_transform_dict['which_terminusQ']
    last_backbond_substid = rot_transform_dict['last_backbond_substid']
    ubiquitin_conjugation_site = rot_transform_dict['ubiquitin_conjugation_site']
    protecting_group = rot_transform_dict['protecting_group']

    ## 
    adaptable_AA = new_atoms['subst_name'].str[:3]
    if which_terminusQ == 'c_terminus':
        base_atoms_max = base_atoms[base_atoms['status_bit']=='BACKBONE']['subst_id'].max()
        base_atoms_min = base_atoms[base_atoms['status_bit']=='BACKBONE']['subst_id'].min()
        new_atoms['subst_id'] = base_atoms_max + 1
        new_atoms['subst_id']= new_atoms['subst_id'].astype(str)
        ## protecting groups have their unique identifier + the amino acid number they are attached t 
        new_atoms['subst_name'] = adaptable_AA + new_atoms['subst_id']
        new_atoms['subst_id']= new_atoms['subst_id'].astype(int)

    elif which_terminusQ == 'n_terminus':
        new_atoms_max = new_atoms[new_atoms['status_bit']=='BACKBONE']['subst_id'].max()
        new_atoms_min = new_atoms[new_atoms['status_bit']=='BACKBONE']['subst_id'].min()
        new_atoms['subst_id'] = new_atoms_min - 1
        new_atoms['subst_id']= new_atoms['subst_id'].astype(str)
        ## protecting groups have their unique identifier + the amino acid number they are attached t 
        new_atoms['subst_name'] = adaptable_AA + new_atoms['subst_id']
        new_atoms['subst_id']= new_atoms['subst_id'].astype(int)

    ## this could be cleaner... 
    elif which_terminusQ == 'lysine':
        new_atoms['subst_id'] = new_atoms['subst_id']+ last_backbond_substid
        new_atoms['subst_id']= new_atoms['subst_id'].astype(str)
        ## protecting groups have their unique identifier + the amino acid number they are attached t 
        new_atoms['subst_name'] = adaptable_AA + new_atoms['subst_id']
        new_atoms['subst_id']= new_atoms['subst_id'].astype(int)

    elif which_terminusQ == 'protecting_group': 
        new_atoms['subst_id'] = ubiquitin_conjugation_site +1000000
        new_atoms['subst_id']= new_atoms['subst_id'].astype(str)
        ## protecting groups have their unique identifier + the amino acid number they are attached to
        new_atoms['subst_name'] = protecting_group + '_' + str(ubiquitin_conjugation_site)
        new_atoms['subst_id']= new_atoms['subst_id'].astype(int)
        
    return new_atoms


def reformat_subst_ids(atoms):
    
    ## add something from protecting groups... 

    ## code for things that arent protecting groups..
    adaptable_AA = atoms['subst_name'].str[:3]

    ## get the subst_ids
    subst_ids = list(atoms['subst_id'].unique())

    ## provide new subst_ids
    new_subst_ids = list(range(1,len(subst_ids)+1))
    
    ## mapping dict
    subst_id_mapping_dict = {}
    for old, new in zip(subst_ids, new_subst_ids):
        subst_id_mapping_dict[old] = new

    ## apply mapping
    atoms["subst_id"]=atoms["subst_id"].map(subst_id_mapping_dict)

    ## new change the subst_ids 
    atoms['subst_id']= atoms['subst_id'].astype(str)
    ## protecting groups have their unique identifier + the amino acid number they are attached t 
    atoms['subst_name'] = adaptable_AA + atoms['subst_id']
    atoms['subst_id']= atoms['subst_id'].astype(int)

    return atoms



## NEED TO DO: 
# make function compatible with c & n terminus addition 
# make function compatible with PG addition 
# make function compatible with ubi addition 

# output from find_peptide_backbone should be atoms, bonds and dictionary 
# bonding_base_atom_id etc etc should be defined at the start... 

## which_terminusQ can any of the lysines... K6 to K63... 
## and maybe some logic around peptide vs protein... 
## inputs :
# multimer_atoms_start
# multimer_atoms_end
# ubiquitin_conjugation_site
# which_terminusQ 
# base_molecule = peptide/protein
# new_molecule = peptide/protein

def terminus_addition(atoms, bonds, new_atoms, new_bonds, json_dictionary):

    ## pull out the variables 
    multimer_atoms_start = json_dictionary['multimer_atoms_start']
    multimer_atoms_end = json_dictionary['multimer_atoms_end']
    ubiquitin_conjugation_site = json_dictionary['ubiquitin_conjugation_site']
    which_terminusQ = json_dictionary['which_terminusQ']
    protecting_group = json_dictionary['protecting_group']
    base_molecule = json_dictionary['base_molecule']
    new_molecule = json_dictionary['new_molecule']
    last_backbond_substid = json_dictionary['last_backbond_substid']

    # define the new ubiquitin... and find the lysine of choice... 
    ## copy so that the value isnt assigned...
    ## this is be a starting place.. 

    ## adaptable atoms is a copy of the ubiquitin of choice...
    ## this will change to be a subset of the original ubiquitin...
    ## adaptable ubi is the ubi unit being added...
    ## this will not be the case...
    base_atoms = atoms.copy()
    base_bonds = bonds.copy()

    ## adaptable atoms is a copy of the ubiquitin of choice...
    ## this will change to be a subset of the original ubiquitin...
    ## adaptable ubi is the ubi unit being added...
    new_atoms = new_atoms.copy()
    new_bonds = new_bonds.copy()

    ## working glycine belongs to the new ubiquitin... never changes...
    ## need to find centre of mass for this molecule...
    ### find the sc terminal glycine..
    base_atoms_max = base_atoms[base_atoms['status_bit']=='BACKBONE']['subst_id'].max()
    base_atoms_min = base_atoms[base_atoms['status_bit']=='BACKBONE']['subst_id'].min()
    new_atoms_max = new_atoms[new_atoms['status_bit']=='BACKBONE']['subst_id'].max()
    new_atoms_min = new_atoms[new_atoms['status_bit']=='BACKBONE']['subst_id'].min()
    
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

    if which_terminusQ == 'c_terminus':
        base_atoms_substid = base_atoms_max
        multimer_atoms_start = base_atoms_min
        multimer_atoms_end = base_atoms_max
        ubiquitin_conjugation_site = ''
    elif which_terminusQ == 'n_terminus':
        base_atoms_substid = base_atoms_min
        multimer_atoms_start = base_atoms_min
        multimer_atoms_end = base_atoms_max
        ubiquitin_conjugation_site = ''
    
    ########################################################################################################
    ########################################################################################################

    if (which_terminusQ == 'c_terminus') | (which_terminusQ == 'n_terminus'):

        ## change the below to a dictionary
        base_atoms, base_bonds, base_result_dictionary =  find_peptide_backbone(base_atoms, base_bonds, aa_position=which_terminusQ, subst_id = base_atoms_substid)    
        
        print(base_result_dictionary)

        base_peptide_backboneQ = base_result_dictionary['peptide_backboneQ']
        base_backbone_atom_id_list = base_result_dictionary['backbone_atom_id_list']
        base_nitrogen_H_atom_id_list = base_result_dictionary['nitrogen_H_atom_id_list']
        base_o3_atom_id_list = base_result_dictionary['o3_atom_id_list']
        base_co2_atom_id_list = base_result_dictionary['co2_atom_id_list']

        ## find the backbone of the 
        N_atom_id_base = base_backbone_atom_id_list[0]
        C3_atom_id_base = base_backbone_atom_id_list[1]
        C2_atom_id_base = base_backbone_atom_id_list[2]
        O2_atom_id_base = base_backbone_atom_id_list[3]
        if which_terminusQ == 'c_terminus':
            O3_atom_id_base = base_o3_atom_id_list[0]

        ## change the below to a dictionary
        new_atoms, new_bonds, new_result_dictionary =  find_peptide_backbone(new_atoms, new_bonds, aa_position='individual', subst_id = new_atoms_min)    
    

        new_peptide_backboneQ = new_result_dictionary['peptide_backboneQ']
        new_backbone_atom_id_list = new_result_dictionary['backbone_atom_id_list']
        new_nitrogen_H_atom_id_list = new_result_dictionary['nitrogen_H_atom_id_list']
        new_o3_atom_id_list = new_result_dictionary['o3_atom_id_list']
        new_co2_atom_id_list = new_result_dictionary['co2_atom_id_list']

        ## find the backbone atom ids 
        N_atom_id_new = new_backbone_atom_id_list[0]
        C3_atom_id_new = new_backbone_atom_id_list[1]
        C2_atom_id_new = new_backbone_atom_id_list[2]
        O2_atom_id_new = new_backbone_atom_id_list[3]
        O3_atom_id_new = new_o3_atom_id_list[0]
    
    ########################################################################################################

    ## this can be a number...
    elif which_terminusQ == 'lysine':

        ## get the specific multimer and lysine...
        conjugating_atoms = base_atoms[(multimer_atoms_start <= base_atoms['subst_id']) &  (multimer_atoms_end  >= base_atoms['subst_id'])]
        working_lysine = conjugating_atoms[conjugating_atoms['subst_id'] == ubiquitin_conjugation_site]

        ## get the value of NZ
        index_lys_NZ = working_lysine[working_lysine['atom_name'] == 'NZ'].index[0]

        # get the atom_ids 
        NZ_atom_id_base = conjugating_atoms.at[index_lys_NZ, 'atom_id']

        # find the atom ids of the atoms next to NZ   
        # find the non-H atom_id ie C3
        nitro_atom_ids, nitro_atom_types = find_bonded_atoms(NZ_atom_id_base, base_atoms, base_bonds)
        C3_previous_list_index = [index for index in range(len(nitro_atom_types)) if nitro_atom_types[index] == 'C.3'] 
        C3_previous_atom_id_base = nitro_atom_ids[C3_previous_list_index[0]]

        # find the atom ids of the atoms next to NZ         
        # find the non-H, non-N atom_id ie C2 previous 
        C3_atom_ids, C3_atom_types = find_bonded_atoms(C3_previous_atom_id_base, base_atoms, base_bonds)
        C3_previous2_list_index = [index for index in range(len(C3_atom_types)) if C3_atom_types[index] == 'C.3'] 
        C3_previous2_atom_id_base = C3_atom_ids[C3_previous2_list_index[0]]

        ## change the below to a dictionary
        new_atoms, new_bonds, new_result_dictionary =  find_peptide_backbone(new_atoms, new_bonds, aa_position='c_terminus', subst_id = new_atoms_max)    
        
        new_peptide_backboneQ = new_result_dictionary['peptide_backboneQ']
        new_backbone_atom_id_list = new_result_dictionary['backbone_atom_id_list']
        new_nitrogen_H_atom_id_list = new_result_dictionary['nitrogen_H_atom_id_list']
        new_o3_atom_id_list = new_result_dictionary['o3_atom_id_list']
        new_co2_atom_id_list = new_result_dictionary['co2_atom_id_list']

        ## find the backbone atom ids 
        N_atom_id_new = new_backbone_atom_id_list[0]
        C3_atom_id_new = new_backbone_atom_id_list[1]
        C2_atom_id_new = new_backbone_atom_id_list[2]
        O2_atom_id_new = new_backbone_atom_id_list[3]
        O3_atom_id_new = new_o3_atom_id_list[0]

    elif which_terminusQ == 'protecting_group':
        ## BASE 
        ## get the specific multimer and lysine...
        conjugating_atoms = base_atoms[(multimer_atoms_start <= base_atoms['subst_id']) &  (multimer_atoms_end  >= base_atoms['subst_id'])]
        working_lysine = conjugating_atoms[conjugating_atoms['subst_id'] == ubiquitin_conjugation_site]
        print('ubiquitin_conjugation_site: ' +str(ubiquitin_conjugation_site))
        print('ubiquitin_conjugation_site: ' +str(working_lysine))


        ## get the value of NZ
        index_lys_NZ = working_lysine[working_lysine['atom_name'] == 'NZ'].index[0]

        # get the atom_ids 
        NZ_atom_id_base = conjugating_atoms.at[index_lys_NZ, 'atom_id']

        # find the atom ids of the atoms next to NZ   
        # find the non-H atom_id ie C3
        nitro_atom_ids, nitro_atom_types = find_bonded_atoms(NZ_atom_id_base, base_atoms, base_bonds)
        C3_previous_list_index = [index for index in range(len(nitro_atom_types)) if nitro_atom_types[index] == 'C.3'] 
        C3_previous_atom_id_base = nitro_atom_ids[C3_previous_list_index[0]]

        # find the atom ids of the atoms next to NZ         
        # find the non-H, non-N atom_id ie C2 previous 
        C3_atom_ids, C3_atom_types = find_bonded_atoms(C3_previous_atom_id_base, base_atoms, base_bonds)
        C3_previous2_list_index = [index for index in range(len(C3_atom_types)) if C3_atom_types[index] == 'C.3'] 
        C3_previous2_atom_id_base = C3_atom_ids[C3_previous2_list_index[0]]

        ## change the below to a dictionary
        index_LG = new_atoms[new_atoms['atom_name'] == 'LG'].index[0]
        LG_atom_id_new = new_atoms.at[index_LG, 'atom_id']
        LG_atom_ids, LG_atom_types = find_bonded_atoms(LG_atom_id_new, new_atoms, new_bonds)
        C2_previous_list_index = [index for index in range(len(LG_atom_types)) if LG_atom_types[index] == 'C.2'] 
        C2_previous_atom_id_new = LG_atom_ids[C2_previous_list_index[0]]
        

    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

    rot_transform_dict = {}

    #### is it c or n terminus
    ## function inputs -- here you could probably change for c vs n terminus... 
    if which_terminusQ == 'c_terminus': 

        nitrogen_on = 'new'        
        ## first you do the coordinate changes then you remove the bonds... 
        ## check that there are three H's in nitrogen_H_atom_id_list
        if len(new_nitrogen_H_atom_id_list) == 3: 
            ## while H1 is 
            Ha_removed_atom_id = new_nitrogen_H_atom_id_list[-3]
            ## H2 is just removed 
            Hb_removed_atom_id = new_nitrogen_H_atom_id_list[-2]
        else: 
            TypeError('aksjdhfkjahsdkhf')
        
        ## define the bonding 
        bonding_base_atom_id  = C2_atom_id_base
        lost_base_atom_id = O3_atom_id_base
        previous_base_atom_id = C3_atom_id_base
        previous2_base_atom_id = N_atom_id_base

        bonding_new_atom_id = N_atom_id_new 
        lost_new_atom_id = Ha_removed_atom_id
        previous_new_atom_id = C3_atom_id_new
        previous2_new_atom_id = C2_atom_id_new

        ## define the add to dictionary 
        rot_transform_dict['bonding_base_atom_id'] = bonding_base_atom_id 
        rot_transform_dict['lost_base_atom_id'] = lost_base_atom_id
        rot_transform_dict['previous_base_atom_id'] = previous_base_atom_id
        rot_transform_dict['previous2_base_atom_id'] = previous2_base_atom_id
        rot_transform_dict['bonding_new_atom_id'] = bonding_new_atom_id
        rot_transform_dict['lost_new_atom_id'] = lost_new_atom_id
        rot_transform_dict['previous_new_atom_id'] = previous_new_atom_id
        rot_transform_dict['previous2_new_atom_id'] = previous2_new_atom_id
        rot_transform_dict['nitrogen_on'] = nitrogen_on
        
        ## add in the other values.... 
        rot_transform_dict['multimer_atoms_start'] = multimer_atoms_start
        rot_transform_dict['multimer_atoms_end'] = multimer_atoms_end
        rot_transform_dict['ubiquitin_conjugation_site'] = ubiquitin_conjugation_site
        rot_transform_dict['which_terminusQ'] = which_terminusQ
        rot_transform_dict['protecting_group'] = protecting_group
        rot_transform_dict['base_molecule'] = base_molecule
        rot_transform_dict['new_molecule'] = new_molecule
        rot_transform_dict['last_backbond_substid'] = last_backbond_substid


        ## SPECIAL CHANGES...    
        # loss of extra H (Hb_removed_atom_id) both in new bonds and new atoms 
        new_atoms, new_bonds = remove_atom_from_bonds_atoms(Hb_removed_atom_id, new_atoms, new_bonds)

        # O.2 is called 'O'
        O2_index_base = base_atoms[base_atoms['atom_id'] == O2_atom_id_base].index[0]
        base_atoms.at[O2_index_base, 'atom_name'] = 'O'

    ########################################################################################################

    elif which_terminusQ == 'n_terminus': 

        nitrogen_on = 'base'

        ## first you do the coordinate changes then you remove the bonds... 
        ## check that there are three H's in nitrogen_H_atom_id_list
        if len(base_nitrogen_H_atom_id_list) == 3: 
            ## while H1 is 
            Ha_removed_atom_id = base_nitrogen_H_atom_id_list[-3]
            ## H2 is just removed 
            Hb_removed_atom_id = base_nitrogen_H_atom_id_list[-2]
        else: 
            TypeError('aksjdhfkjahsdkhf')

        bonding_base_atom_id  = N_atom_id_base
        lost_base_atom_id = Ha_removed_atom_id
        previous_base_atom_id = C3_atom_id_base
        previous2_base_atom_id = C2_atom_id_base
        
        bonding_new_atom_id = C2_atom_id_new 
        lost_new_atom_id = O3_atom_id_new
        previous_new_atom_id = C3_atom_id_new
        previous2_new_atom_id = N_atom_id_new

        ## define the add to dictionary 
        rot_transform_dict['bonding_base_atom_id'] = bonding_base_atom_id 
        rot_transform_dict['lost_base_atom_id'] = lost_base_atom_id
        rot_transform_dict['previous_base_atom_id'] = previous_base_atom_id
        rot_transform_dict['previous2_base_atom_id'] = previous2_base_atom_id
        rot_transform_dict['bonding_new_atom_id'] = bonding_new_atom_id
        rot_transform_dict['lost_new_atom_id'] = lost_new_atom_id
        rot_transform_dict['previous_new_atom_id'] = previous_new_atom_id
        rot_transform_dict['previous2_new_atom_id'] = previous2_new_atom_id
        rot_transform_dict['nitrogen_on'] = nitrogen_on

        ## add in the other values.... 
        rot_transform_dict['multimer_atoms_start'] = multimer_atoms_start
        rot_transform_dict['multimer_atoms_end'] = multimer_atoms_end
        rot_transform_dict['ubiquitin_conjugation_site'] = ubiquitin_conjugation_site
        rot_transform_dict['which_terminusQ'] = which_terminusQ
        rot_transform_dict['protecting_group'] = protecting_group
        rot_transform_dict['base_molecule'] = base_molecule
        rot_transform_dict['new_molecule'] = new_molecule
        rot_transform_dict['last_backbond_substid'] = last_backbond_substid

        ## SPECIAL CHANGES...    
        # loss of extra H (Hb_removed_atom_id) both in new bonds and new atoms 
        base_atoms, base_bonds = remove_atom_from_bonds_atoms(Hb_removed_atom_id, base_atoms, base_bonds)
        
        # check the O2 and O3 dont have any H on it... 
        O2_index_new = new_atoms[new_atoms['atom_id'] == O2_atom_id_new].index[0]
        new_atoms.at[O2_index_new, 'atom_name'] = 'O'

    ########################################################################################################

    elif which_terminusQ == 'lysine':
        
        ## check nitrogen hydrogen list
        nitrogen_on = 'base'
        
        Ha_removed_atom_id, base_atoms, base_bonds = fix_bonding_nitrogen(NZ_atom_id_base, base_atoms, base_bonds)
        
        bonding_base_atom_id  = NZ_atom_id_base
        lost_base_atom_id = Ha_removed_atom_id
        previous_base_atom_id = C3_previous_atom_id_base
        previous2_base_atom_id = C3_previous2_atom_id_base
        
        bonding_new_atom_id = C2_atom_id_new 
        lost_new_atom_id = O3_atom_id_new
        previous_new_atom_id = C3_atom_id_new
        previous2_new_atom_id = N_atom_id_new

        ## define the add to dictionary 
        rot_transform_dict['bonding_base_atom_id'] = bonding_base_atom_id 
        rot_transform_dict['lost_base_atom_id'] = lost_base_atom_id
        rot_transform_dict['previous_base_atom_id'] = previous_base_atom_id
        rot_transform_dict['previous2_base_atom_id'] = previous2_base_atom_id
        rot_transform_dict['bonding_new_atom_id'] = bonding_new_atom_id
        rot_transform_dict['lost_new_atom_id'] = lost_new_atom_id
        rot_transform_dict['previous_new_atom_id'] = previous_new_atom_id
        rot_transform_dict['previous2_new_atom_id'] = previous2_new_atom_id
        rot_transform_dict['nitrogen_on'] = nitrogen_on

        ## add in the other values.... 
        rot_transform_dict['multimer_atoms_start'] = multimer_atoms_start
        rot_transform_dict['multimer_atoms_end'] = multimer_atoms_end
        rot_transform_dict['ubiquitin_conjugation_site'] = ubiquitin_conjugation_site
        rot_transform_dict['protecting_group'] = protecting_group
        rot_transform_dict['which_terminusQ'] = which_terminusQ
        rot_transform_dict['base_molecule'] = base_molecule
        rot_transform_dict['new_molecule'] = new_molecule
        rot_transform_dict['last_backbond_substid'] = last_backbond_substid

        # check the O2 and O3 dont have any H on it... 
        O2_index_new = new_atoms[new_atoms['atom_id'] == O2_atom_id_new].index[0]
        new_atoms.at[O2_index_new, 'atom_name'] = 'O'

    elif which_terminusQ == 'protecting_group':
        ## check nitrogen hydrogen list
        nitrogen_on = 'base'

        Ha_removed_atom_id, base_atoms, base_bonds = fix_bonding_nitrogen(NZ_atom_id_base, base_atoms, base_bonds)
        
        bonding_base_atom_id  = NZ_atom_id_base
        lost_base_atom_id = Ha_removed_atom_id
        previous_base_atom_id = C3_previous_atom_id_base
        previous2_base_atom_id = C3_previous2_atom_id_base

        bonding_new_atom_id = C2_previous_atom_id_new 
        lost_new_atom_id = LG_atom_id_new
        previous_new_atom_id = LG_atom_id_new
        previous2_new_atom_id = C2_previous_atom_id_new

        ## define the add to dictionary 
        rot_transform_dict['bonding_base_atom_id'] = bonding_base_atom_id 
        rot_transform_dict['lost_base_atom_id'] = lost_base_atom_id
        rot_transform_dict['previous_base_atom_id'] = previous_base_atom_id
        rot_transform_dict['previous2_base_atom_id'] = previous2_base_atom_id
        rot_transform_dict['bonding_new_atom_id'] = bonding_new_atom_id
        rot_transform_dict['lost_new_atom_id'] = lost_new_atom_id
        rot_transform_dict['previous_new_atom_id'] = previous_new_atom_id
        rot_transform_dict['previous2_new_atom_id'] = previous2_new_atom_id
        rot_transform_dict['nitrogen_on'] = nitrogen_on

        ## add in the other values.... 
        rot_transform_dict['multimer_atoms_start'] = multimer_atoms_start
        rot_transform_dict['multimer_atoms_end'] = multimer_atoms_end
        rot_transform_dict['ubiquitin_conjugation_site'] = ubiquitin_conjugation_site
        rot_transform_dict['which_terminusQ'] = which_terminusQ
        rot_transform_dict['protecting_group'] = protecting_group
        rot_transform_dict['base_molecule'] = base_molecule
        rot_transform_dict['new_molecule'] = new_molecule
        rot_transform_dict['last_backbond_substid'] = last_backbond_substid


    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

    ## MOL2 changes... 
    ## call the rotation transformation function
    base_atoms, base_bonds, new_atoms, new_bonds = rotation_transformation_mol2(base_atoms, base_bonds, new_atoms, new_bonds, rot_transform_dict)

    ## call the bonding function -- where is the nitrogen on? new or base? this case new...
    ## change N.am in bonds and atoms... 
    base_atoms, base_bonds, new_atoms, new_bonds = connectivity_changes_amide_bond(base_atoms, base_bonds, new_atoms, new_bonds, rot_transform_dict)


    ########################################################################################################
    ########################################################################################################

    ## call the subst_id change function
    new_atoms = new_atom_subst_ids(new_atoms, base_atoms, rot_transform_dict)    

    ## concatenate the dataframes of atoms and bonds
    if which_terminusQ == 'c_terminus': 
        concat_atoms = pd.concat([base_atoms, new_atoms])
        concat_bonds = pd.concat([base_bonds, new_bonds])
        ## need to shift everything with c_terminus & n_terminus
        concat_atoms = reformat_subst_ids(concat_atoms) 
    elif which_terminusQ == 'n_terminus': 
        concat_atoms = pd.concat([new_atoms, base_atoms])
        concat_bonds = pd.concat([new_bonds, base_bonds])
        ## need to shift everything with c_terminus & n_terminus
        concat_atoms = reformat_subst_ids(concat_atoms) 
    elif which_terminusQ == 'lysine': 
        concat_atoms = pd.concat([base_atoms, new_atoms])
        concat_bonds = pd.concat([base_bonds, new_bonds])
    elif which_terminusQ == 'protecting_group': 
        concat_atoms = pd.concat([base_atoms, new_atoms])
        concat_bonds = pd.concat([base_bonds, new_bonds])

    ## call the atom_id change function
    concat_atoms, concat_bonds = reformat_atom_bond_ids(concat_atoms, concat_bonds)

    ## call subst_id reordering function
    return concat_atoms, concat_bonds

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

def write_mol2(final_ubiAllAtoms,final_ubiAllBonds, path):
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
                                    'ubi_multimer.pdb', 
                                    str(maxAtomID) + ' ' + str(maxBondID),
                                    'SMALL',
                                    'NO_CHARGES',
                                    '@<TRIPOS>ATOM'
                                    ]})

    pandas_for_mol = pd.concat([mol2Header, atoms_combined, pd.DataFrame({'C': ['@<TRIPOS>BOND']}), bonds_combined, pd.DataFrame({'C': ['@<TRIPOS>SUBSTRUCTURE', '0 pdb 0 ****']}                                                                                                                                                                                                                                                                                             )],ignore_index=True)
    np.savetxt(path, pandas_for_mol, delimiter="", fmt="%s") 
    os.rename(path, path.replace(".txt", ".mol2"))