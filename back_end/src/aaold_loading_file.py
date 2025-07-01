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
import logging
import os
import pandas as pd
import itertools
import pandas as pd
import xlsxwriter
### this should all go before the bonding and atom changes...
## as bonding and atom changes delete atoms..

from scipy.spatial.transform import Rotation as R
import numpy as np
import math
import json
import copy
import sys
import logging
from back_end.src.aaold_loading_file import * 
import kaleido
## image at the top, followed by reaction & monomer history....
import json
import copy
import ipywidgets as widgets
import plotly.graph_objs as go
import pandas as pd
from plotly.validators.scatter.marker import SymbolValidator
import networkx as nx
import random
import pandas as pd
### for adding ubiquitin image to description...
import plotly.express as px
import pycountry
import numpy as np
from PIL import Image
import sys, os
import plotly.graph_objects as go
import numpy as np
import reportlab
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import Frame, Image

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
pdfmetrics.registerFont(TTFont('Vera', 'Vera.ttf'))
pdfmetrics.registerFont(TTFont('VeraBd', 'VeraBd.ttf'))
pdfmetrics.registerFont(TTFont('VeraIt', 'VeraIt.ttf'))
pdfmetrics.registerFont(TTFont('VeraBI', 'VeraBI.ttf'))
import sys
from pypdf import PdfWriter

# Standard Libraries
import math
import string
from operator import itemgetter

# External Libraries
import six
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D


#### Graphing and PDF functions 

## maybe have the input functions be have everything labelled... 
## and then you add at at a particular chain_number and lysine site.. 
def creating_networkx_nodes_edges(parent_dictionary):
    
    if isinstance(parent_dictionary, str):
        x = parent_dictionary.replace("'", "\"")
        parent_dictionary = json.loads(x)
    
    global chain_number_list
    chain_number_list = [1]
    
    # for node information
    global networkx_NODES
    global networkx_EDGES
    networkx_NODES = []
    networkx_EDGES = []
    
    # for attribute information
    global ubi_nodes
    ubi_nodes = []
    # lysine_NODES = K48_nodes + K63_nodes
    # pg_NODES = SMAC_nodes + ABOC_nodes
    global K48_nodes
    global K11_nodes
    global K63_nodes
    global SMAC_nodes
    global ABOC_nodes
    K48_nodes = []
    K11_nodes = []
    K63_nodes = []
    SMAC_nodes = []
    ABOC_nodes = []

    # for attribute information
    global ubi_lysine_edges
    global lysine_PG_edges
    global lysine_ubi_edges
    ubi_lysine_edges = []
    lysine_PG_edges = []
    lysine_ubi_edges = []

    inner_wrapper_creating_networkx_nodes_edges(parent_dictionary)
    return {"networkx_NODES" : networkx_NODES,
            "networkx_EDGES" : networkx_EDGES, 
            "ubi_nodes" : ubi_nodes, 
            "K48_nodes" : K48_nodes, 
            "K63_nodes" : K63_nodes, 
            "SMAC_nodes" : SMAC_nodes,
            "ABOC_nodes" : ABOC_nodes,
            "ubi_lysine_edges" : ubi_lysine_edges,
            "lysine_PG_edges" : lysine_PG_edges,
            "lysine_ubi_edges" : lysine_ubi_edges}


### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_creating_networkx_nodes_edges(input_dictionary, previous_connection = 'root'):  
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)
    if ('protein' in working_dictionary.keys()):
        logging.info(' ===== START OF PROTEIN ===== ')
        logging.info("Protein: " + working_dictionary['protein'])
        working_dictionary['chain_number'] = chain_number_list[-1]
        logging.info("chain_number_list: " + str(chain_number_list))
        logging.info("UBI NUMBER: " + str(working_dictionary['chain_number']))
        logging.info("LYSINES: " + str(working_dictionary['branching_sites']))
        
        ### add every ubiquitin... and its number... 
        ### concatenate the two... 
        ### if start is root you can ignore... 
        logging.info("===== NODE CHANGES: START  =====")
        #ubinode = str(working_dictionary['protein']) + '_' + str(working_dictionary['chain_number'])
        ubinode = str(working_dictionary['chain_number'])

        ### append to: networkx_NODE_list, ubi_nodes, networkx_EDGE_list
        networkx_NODES.append(ubinode)
        ubi_nodes.append(ubinode)
        if previous_connection != 'root':
            networkx_EDGES.append((previous_connection, ubinode))
            lysine_ubi_edges.append((previous_connection, ubinode))
        logging.info("===== networkx_NODES:" + str(networkx_NODES))
        logging.info("===== ubi_nodes:" + str(ubi_nodes))
        logging.info("===== networkx_EDGES:" + str(networkx_EDGES))
        logging.info("===== lysine_ubi_edges:" + str(lysine_ubi_edges))
        logging.info("===== NODE CHANGES: END  =====")

        ### associate branching_sites for the next for looop... 
        working_branching_sites = working_dictionary['branching_sites']
        ### while still working witht the last ubi you add to the list the next ubi position you will label...
        chain_number_list.append((chain_number_list[-1]+1))
        
        ### adding branching_sites... here can use the chain_number and lysine site...
        ### working_branching_sites = loop_through_immediate_branching_sites(working_branching_sites, x_ubi_num)
        ### working_dictionary['branching_sites'] = working_branching_sites
        for i in working_branching_sites:
            if str(i['site_name']) in ['K48', 'K63']:
                logging.info(' ===== START OF LYSINE SITE =====  ')
                logging.info("UBI NUMBER: " + str(working_dictionary['chain_number']))
                ### which site K48 or K63
                logging.info("LYSINE SITE: " + str(i['site_name']))
                
                logging.info("===== NODE CHANGES: START  =====")
                lysine = str(working_dictionary['chain_number'])+ "_"+ str(i['site_name'])
                ### append to: networkx_NODE_list
                networkx_NODES.append(lysine)
                networkx_EDGES.append((ubinode, lysine))
                ubi_lysine_edges.append((ubinode, lysine))
                logging.info("===== networkx_NODES:" + str(networkx_NODES))
                logging.info("===== networkx_EDGES:" + str(networkx_EDGES))
                logging.info("===== ubi_lysine_edges:" + str(ubi_lysine_edges))

                ### add to the K48 or K63 node list
                if str(i['site_name']) == 'K48':
                    K48_nodes.append(lysine)
                    logging.info("===== K48_nodes:" + str(K48_nodes))
                     
                elif str(i['site_name']) == 'K63':
                    K63_nodes.append(lysine)                
                    logging.info("===== K63_nodes:" + str(K63_nodes)) 

                logging.info("===== NODE CHANGES: END  =====")

            ### if the site has a protecting group or a protein 
                ### if the site is a protecting group
                if (i['children'] =='SMAC' or i['children'] =='ABOC'):
                    logging.info("Protecting Group: " + str(i['children']))

                    logging.info("==== ALERT: START OF OUTPUTS OF INTEREST ===== ")
                    logging.info(str(working_dictionary['chain_number'])+ "_"+ str(i['site_name'])+ "_"+ str(i['children']))
                    logging.info("==== ALERT: END OF OUTPUTS OF INTEREST ===== ")
                    
                    logging.info("===== NODE CHANGES: START  =====")
                    protecting_group = str(working_dictionary['chain_number'])+ "_"+ str(i['site_name'])+ "_"+ str(i['children'])
                    ### append to: networkx_NODE_list
                    networkx_NODES.append(protecting_group)
                    networkx_EDGES.append((lysine, protecting_group))
                    lysine_PG_edges.append((lysine, protecting_group))
                    logging.info("===== networkx_NODES:" + str(networkx_NODES))
                    logging.info("===== networkx_EDGES:" + str(networkx_EDGES))
                    logging.info("===== lysine_PG_edges:" + str(lysine_PG_edges))

                    ### add to the K48 or K63 node list
                    if str(i['children']) == 'SMAC':
                        SMAC_nodes.append(protecting_group)
                        logging.info("===== SMAC_nodes:" + str(SMAC_nodes))

                    elif str(i['children']) == 'ABOC':
                        ABOC_nodes.append(protecting_group)                
                        logging.info("===== K63_nodes:" + str(ABOC_nodes))    
                    

                ### if the site is a protein.... 
                elif i['children'] =='': 
                    ### recursive function that calls itself... 
                    logging.info("===== NODE CHANGES: END  =====")

                ### if the site is a protein.... 
                elif isinstance(i['children'], dict): 
                    ### recursive function that calls itself... 
                    logging.info('NEXT UBIQUITIN: ' + str(i['children']))
                    
                    #### figure out a way to carry the lysine through... perhaps use 'root' as the default
                    ## (str(working_dictionary['protein'])+ str(working_dictionary['chain_number'])+ "_"+ str(i['site_name'])) = the last lysine
                    # so an edge can be created 
                    inner_wrapper_creating_networkx_nodes_edges(i['children'], (str(working_dictionary['chain_number'])+ "_"+ str(i['site_name']))) 
                logging.info(' ===== END OF LYSINE SITE =====  ')
        logging.info(' ===== END OF PROTEIN - UBI NUMBER: ' + str(working_dictionary['chain_number']) + ' =====  ')
    return working_dictionary

### input.. graph, root node, and width...
### output.. 


#### code below provides the x and y coordinates for the tree
def hierarchy_pos(G, root=None, width=0.2, vert_gap = 0.4, vert_loc = 0, xcenter = 0.5):
    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.  
    Licensed under Creative Commons Attribution-Share Alike 
    
    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.
    
    G: the graph (must be a tree)
    
    root: the root node of current branch 
    - if the tree is directed and this is not given, 
      the root will be found and used
    - if the tree is directed and this is given, then 
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given, 
      then a random choice will be used.
    
    width: horizontal space allocated for this branch - avoids overlap with other branches
    
    vert_gap: gap between levels of hierarchy
    
    vert_loc: vertical location of root
    
    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
            logging.info(root)
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=3., vert_gap = 0.3, vert_loc = 0, xcenter = 0.5, pos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''
    
        if pos is None:
            pos = {root:(xcenter,vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)

        ### may need to change this code to include changes when there are more than 2 children... up to 7... 
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            dx = width/len(children) 
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                ### to change whether the graph goes up or down you do vert_loc - vert_gap
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap, 
                                    vert_loc = vert_loc+vert_gap, xcenter=nextx,
                                    pos=pos, parent = root)
        return pos

            
    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)


### defining the x_y coordinates of the nodes
def x_y_coordinates_for_nodes(G, pos):
    edge_x = []
    edge_y = []
    for edge in list(G.edges):
        ## take coordinates from hierachal code above
        ## take the list from the edges... loop through them and pull their coordinates from the hierachal functions output = pos
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        ## gives the coordinates... 
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in list(G.nodes):
        ### here you can add the node positionss.. 
        ## take the list from the nodes... loop through them and pull their coordinates from the hierachal functions output = pos
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    ### separate traces for each species..
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        textposition='middle right',
        marker=dict(
            ### https://plotly.com/python/marker-style/ for marker style
            showscale=False,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            ## no list gives a marker to the whole thing
            opacity=1,
            line=dict(width=2, color="DarkSlateGrey"),
            ## list allow you to change each marker individually 
            symbol=[], 
            color=[],
            size=[]))
    
    #### OUTPUT node_trace && edge_trace
    return (node_trace, edge_trace)

### code to change the marker information of each node, also, symbol etc...

### input.. G the graph... and networkx_noded_JSON
### output.. node trace... 

def add_graph_annotations(G, node_trace, networkx_noded_JSON):
    
    ### 1= Ubiquitin, 2= K48, 3= K63, 4= SMAC, 5= SMAC
    node_color = []
    node_text = []
    node_size = []
    node_symbol = []

    ## LISTS FOR what is what... 
    K48_nodes = networkx_noded_JSON["K48_nodes"]
    K63_nodes = networkx_noded_JSON["K63_nodes"]
    SMAC_nodes = networkx_noded_JSON["SMAC_nodes"]
    ABOC_nodes = networkx_noded_JSON["ABOC_nodes"] 
    ubi_nodes = networkx_noded_JSON["ubi_nodes"]

    ### text for the the labels in the image..
    for node, adjacencies in enumerate(G.adjacency()):
        if adjacencies[0] in ubi_nodes:
            ubi = str(adjacencies[0])
            ## usual ubi node text has been removed...
            #node_text.append(("UBIQUITIN<br><br>Ubiquitin: " + ubi + "<br>Lysine: N/A <br>Protecting Group : N/A"))
            node_text.append((""))
            ### 1= Ubiquitin, 2= K48, 3= K63, 3= K11, 7= SMAC, 8= SMAC
            node_color.append('#551a8b')
            ### ubi node size can be zero to make the ubiquitin the image of ubiquitin == normally = node_size.append(0)
            node_size.append(25)
            node_symbol.append("circle")
        
        
        elif adjacencies[0] in K48_nodes:
            ubi, lys = str(adjacencies[0]).split("_")
            #node_text.append(("K48<br><br>Ubiquitin: " + ubi + "<br>Lysine: " + lys + "<br>Protecting Group : N/A"))
            node_text.append((""))
            ### 1= Ubiquitin, 2= K48, 3= K63, 3= K11, 7= SMAC, 8= SMAC
            node_color.append('#c7e9b4')
            node_size.append(10)
            node_symbol.append("square")
        
        elif adjacencies[0] in K63_nodes:
            ubi, lys = str(adjacencies[0]).split("_")
            #node_text.append(("K63<br><br>Ubiquitin: " + ubi + "<br>Lysine: " + lys + "<br>Protecting Group : N/A"))
            node_text.append((""))
            ### 1= Ubiquitin, 2= K48, 3= K63, 3= K11, 7= SMAC, 8= SMAC
            node_color.append('#7fcdbb')
            node_size.append(10)
            node_symbol.append("square")
        
        elif adjacencies[0] in SMAC_nodes:
            ubi, lys, pg = str(adjacencies[0]).split("_")
            #node_text.append(("SMAC<br><br>Ubiquitin: " + ubi + "<br>Lysine: " + lys + "<br>Protecting Group: " +pg))
            node_text.append((""))
            ### 1= Ubiquitin, 2= K48, 3= K63, 3= K11, 7= SMAC, 8= SMAC
            node_color.append('#1d91c0')
            node_size.append(8)
            node_symbol.append("diamond")
        
        
        elif adjacencies[0] in ABOC_nodes:
            ubi, lys, pg = str(adjacencies[0]).split("_")
            #node_text.append(("ABOC<br><br>Ubiquitin: " + ubi + "<br>Lysine: " + lys + "<br>Protecting Group: " +pg))
            node_text.append((""))
            ### 1= Ubiquitin, 2= K48, 3= K63, 3= K11, 7= SMAC, 8= SMAC
            node_color.append('#0c2c84')
            node_size.append(8)
            node_symbol.append("diamond")

    node_trace.marker.color = node_color
    node_trace.marker.symbol = node_symbol
    node_trace.marker.size = node_size
    node_trace.text = node_text

    return node_trace

def x_y_coordinates_for_legend_in_node_trace(node_trace):
    ## provide coordinates for the legend... 
    ## provide text... 
    ## provide color to match the nodes.... 
    ## connect to the JSON... 
    #### save the changes made by the clicks... 
    list_of_build_action = [' Ubiquitin',' K48',' K63',' SMAC', ' ABOC']

    ### tuples
    x_pos_ubi = [0.85]*len(list_of_build_action)
    x_pos_legend = tuple(x_pos_ubi)

    y_pos_ubi = list(range(0, len(list_of_build_action)))
    ## create a gap between list of ubi and react
    ## that gap is 1
    y_pos_legend = tuple([x / 5 for x in y_pos_ubi])
    symbol_ubi = ['circle','square', 'square','diamond','diamond']

    legend_text = tuple(list_of_build_action)
    legend_symbol = tuple(symbol_ubi)
    
    ## lists 
    color_ubi = ['#551a8b', '#c7e9b4', '#7fcdbb', '#1d91c0', '#0c2c84']
    size_ubi = [25, 10, 10, 8, 8]

    legend_color = color_ubi
    legend_size = size_ubi

    node_trace.marker.color = tuple(legend_color + list(node_trace.marker.color))
    node_trace.marker.symbol = tuple(list(legend_symbol) + list(node_trace.marker.symbol))
    node_trace.marker.size = tuple(legend_size + list(node_trace.marker.size))
    node_trace.text = tuple(list(legend_text) + list(node_trace.text))
    node_trace.x = tuple(list(x_pos_legend) + list(node_trace['x']))
    node_trace.y = tuple(list(y_pos_legend) + list(node_trace['y']))

    return node_trace

### takes ubi_json file and pulls the node and edge trace... 

def defining_node_and_edge_trace(js_data):
    ### save the js file as your working file... 

    ### https://plotly.com/python/network-graphs/
    ### input.. G the graph... and pos
    ### output.. node trace and edge trace
    ### dictionary with everything described... 
    ## networkx the JSON...
    networkx_noded_JSON = creating_networkx_nodes_edges(js_data)
    nodes = networkx_noded_JSON['networkx_NODES']
    edges = networkx_noded_JSON['networkx_EDGES']

    ## graph out of the nodes and edges... 
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    ### input.. G, root node and the width of the graph
    #pos = hierarchy_pos(G, "Ubi1", 0.8)
    pos = hierarchy_pos(G, "1", 0.8)

    #### input.. G & pos 
    ## OUTPUT node_trace && edge_trace 
    node_trace, edge_trace = x_y_coordinates_for_nodes(G, pos)
    node_trace = add_graph_annotations(G, node_trace, networkx_noded_JSON)

    ## return node trace without legend, but also the one with a legend
    return (node_trace, edge_trace)

### maybe add it as a different trace... remove it from the event... (probably best...)

### clickmode... should help get the thing you wanted regarding a 
#f = go.FigureWidget([go.Scatter(x=x, y=y, mode='markers')])
### input.. node trace... edge trace... 
### output.. figure...
def updatable_ubiquitin_figure(node_trace, edge_trace):
        
    f = go.FigureWidget([node_trace,edge_trace])
    f.update_layout(title='<br>Network graph for Ubiquitin created by Bode Group',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest', 
                    margin=dict(b=20,l=5,r=5,t=40),
                    #annotations=[ dict(
                    #          text="",
                    #          showarrow=False,
                    #          xref="paper", yref="paper",
                    #          x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0., 1.0]),#""", range=[0.35, 0.65]"""
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.2, 4.5]))#""", range=[-0.2, 5]"""£
    
    

    scatter_node = f.data[0]
    #edge_trace
    scatter_edge = f.data[1]

    data = f.data
    ###colors = ['#a3a7e4'] * 100
    ###scatter.marker.color = colors
    ###scatter.marker.size = [10] * 100
    return f

def playing_with_ubiquitin(js_data):
    ### here can include it... where the coloring system of the labels is carried forward...
    node_trace, edge_trace = defining_node_and_edge_trace(js_data)
    ## return the figure.. 
    ### node trace...
    node_trace_with_legend = node_trace
    ### default ubi node...
    node_trace_with_legend = x_y_coordinates_for_legend_in_node_trace(node_trace_with_legend)
    ## return the figure.. 
    return updatable_ubiquitin_figure(node_trace_with_legend, edge_trace)




## change colormap to what you want...
## specific for each combination...
### plotting 96 well plates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## SMAC deprotection on deck...
## so first step on deck will be smac deprot...
## if smac deprot then index = index+1
## you should re-input everything...

def build_acceptor_jsons(largest_multimer_size):
    ## getting the acceptors from the dataframe
    starting_acceptor_size= 1
    next_multimer_size= largest_multimer_size
    ## figure out acceptors to use in the next section...
    acceptors_df = pd.read_csv('/Users/ekummelstedt/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/aa_looping_through_builders/.data/data_for_analysis/acceptor_df.csv')
    multimer_summary = pd.read_csv('/Users/ekummelstedt/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/aa_looping_through_builders/.data/reaction_history/' + str(starting_acceptor_size) + 'mer__to_'+ str(next_multimer_size-1)+'multimer_summary.csv')
    unique_acceptors = list(multimer_summary['1'].map(getting_multimer_string_name).unique())
    unique_acceptors_dict = {}
    for i in unique_acceptors: 
        unique_acceptors_dict[i] = len(multimer_summary[multimer_summary['1'].map(getting_multimer_string_name)==i])
        
    acceptor_json_list = []
    for i in list(unique_acceptors_dict.keys()):
        hi = list(acceptors_df[acceptors_df['acceptor_string'] == i]['acceptor_json'].values)
        acceptor_json_list = acceptor_json_list + hi

    ## add in sohei acceptor stuff, but rename it with the new nomenclature...
    acceptor_json_list= list(set(acceptor_json_list))
    ## sohei deprot first so that sohei without deprot overwrites...
    acceptor_dict = {}
    for json_i, string_i in zip(acceptor_json_list, list(unique_acceptors_dict.keys())):
        acceptor_dict[string_i] = json_i
    path_to_file = ".data/jsons/acceptors.json"

    with open(path_to_file, 'w') as file:
        json_string = json.dumps(acceptor_dict, indent=2)
        file.write(json_string)     

    reverse_acceptor_dict = {}
    for json_i, string_i in zip(acceptor_json_list, list(unique_acceptors_dict.keys())):
        reverse_acceptor_dict[json_i] = string_i
    path_to_file = ".data/jsons/acceptors_reverse.json"
    with open(path_to_file, 'w') as file:
        json_string = json.dumps(reverse_acceptor_dict, indent=2)
        file.write(json_string)

    return acceptor_json_list

def simulate_synthesis(largest_multimer, ubi_acceptor_list, ubi_donor_list):
    
    ## 2nd Simulation...
    enzyme_of_deprot_first = 'enzyme'
    starting_acceptor_history_list = []
    for i in ubi_acceptor_list: 
        # acceptor list
        starting_acceptor_history_list.append([i])
    starting_reaction_history_list= [[""]]*len(starting_acceptor_history_list)
    starting_monomer_history_list= [[""]]*len(starting_acceptor_history_list)
    create_synthesis_dataframes(starting_acceptor_history_list, starting_reaction_history_list, starting_monomer_history_list, ubi_donor_list, enzyme_of_deprot_first= enzyme_of_deprot_first)
    starting_size = 1
    export_multimer_pngs(2)
    export_multimer_pdfs(starting_size,2, enzyme_of_deprot_first=enzyme_of_deprot_first)

    starting_size = 1
    for i in range(3,largest_multimer):
        build_next_history(starting_size, ubi_donor_list, i, enzyme_of_deprot_first)
        starting_size = 1
        
    for i in range(3,largest_multimer):
        export_multimer_pngs(i)
        export_multimer_pdfs(starting_size,i, enzyme_of_deprot_first)

    for i in range(3,largest_multimer):
        merge_multimer_pdfs(starting_size,i)

    acceptor_json_list = build_acceptor_jsons(largest_multimer)

    ## 2nd Simulation...
    enzyme_of_deprot_first = 'deprot'
    starting_acceptor_history_list = []
    for i in acceptor_json_list: 
        # acceptor list
        starting_acceptor_history_list.append([i])
    starting_reaction_history_list= [[""]]*len(starting_acceptor_history_list)
    starting_monomer_history_list= [[""]]*len(starting_acceptor_history_list)

    create_synthesis_dataframes(starting_acceptor_history_list, starting_reaction_history_list, starting_monomer_history_list, ubi_donor_list, enzyme_of_deprot_first)
    # starting_size is defined by the input list of acceptors 
    starting_size = 2
    for i in range(4,largest_multimer):
        build_next_history(starting_size, ubi_donor_list, i, enzyme_of_deprot_first)
    for i in range(3,largest_multimer):
        export_multimer_pngs(i)
        export_multimer_pdfs(starting_size,i, enzyme_of_deprot_first)
    for i in range(3,largest_multimer):
        merge_multimer_pdfs(starting_size, i)

    return


def plot_96wells(figure=1, figure_name = 'Test',colorbar_type= 'PuRd', cdata=None, sdata=None, bdata=None, bcolors=None, bmeans=None, **kwargs):
    # from https://github.com/jaumebonet/RosettaSilentToolbox/blob/master/rstoolbox/plot/experimental.py
    """Plot data of a 96 well plate into an equivalent-shaped plot.

    Allows up to three data sources at the same time comming from experiments performed in a
    `96 wells microplate <https://www.wikiwand.com/en/Microtiter_plate>`_. Two of this data
    sources can have to be numerical, and are represented as color and size, while an extra
    boolean dataset can be represented by the color of the border of each well.

    Plot is based on :func:`~matplotlib.pyplot.scatter`; some graphical properties to control
    the visuals (such as ``cmap``), can be provided through this function.

    :param cdata: Data contentainer to be represented by color coding. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **continuos numerical data**.
    :type cdata: :class:`~pandas.DataFrame`
    :param sdata: Data contentainer to be represented by size coding. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **continuos numerical data**.
    :type sdata: :class:`~pandas.DataFrame`
    :param bdata: Data contentainer to be represented by the edge color. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **boolean data**.
    :type bdata: :class:`~pandas.DataFrame`
    :param bcolors: List of the two colors to identify the differences in the border for binary
        data. It has to be a list of 2 colors only. First color represents :data:`True` instances
        while the second color is :data:`False`. Default is ``black`` and ``green``.
    :type bcolors: :func:`list` of :class:`str`
    :param bmeans: List with the meanings of the boolean condition (for the legend). First color
        represents :data:`True` instances while the second color is :data:`False`. Default is
        ``True`` and ``False``.
    :type bmeans: :func:`list` of :class:`str`

    :return: Union[:class:`~matplotlib.figure.Figure`, :class:`~matplotlib.axes.Axes`]

    :raises:
        :ValueError: If input data is not a :class:`~pandas.DataFrame`.
        :ValueError: If :class:`~pandas.DataFrame` do not has the proper shape.
        :ValueError: If ``bcolors`` of ``bmeans`` are provided with sizes different than 2.

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.plot import plot_96wells
           ...: import numpy as np
           ...: import pandas as pd
           ...: import matplotlib.pyplot as plt
           ...: np.random.seed(0)
           ...: df = pd.DataFrame(np.random.randn(8, 12))
           ...: fig, ax = plot_96wells(cdata = df, sdata = -df, bdata = df<0)
           ...: plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

        @savefig plot_96wells_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """

    # Changes in one of this parameters should change others to ensure size fit.
    fig = plt.figure(figure, figsize=(15, 7))
    fig.suptitle(figure_name, fontsize=16)
    ax  = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    top = 2000
    bot = 50
    sizeOfFont = 15
    ticks_font = font_manager.FontProperties(style='normal', size=sizeOfFont, weight='normal')

    # Fixed: THIS CANNOT CHANGE!
    kwargs['x'] = list(range(1, 13)) * 8
    kwargs['y'] = sorted(list(range(1, 9)) * 12)

    # Modified by the input data.
    kwargs.setdefault('s', [top, ] * len(kwargs['y']))
    kwargs.setdefault('c', 'white')
    kwargs.setdefault('edgecolor', ['black', ] * len(kwargs['y']))
    kwargs.setdefault('linewidths', 1.5)
    kwargs.setdefault('cmap', colorbar_type)

    # Color Data
    if cdata is not None:
        if not isinstance(cdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if cdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        kwargs['c'] = cdata.values.flatten()

    # Size Data
    if sdata is not None:
        if not isinstance(sdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if sdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        p = sdata.values.flatten()
        p = ((p - np.min(p)) / np.ptp(p)) * (top - bot) + bot
        kwargs['s'] = p

    # Border Data
    if bdata is not None:
        if not isinstance(bdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if bdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        if not (bdata.dtypes == bool).all():
            raise ValueError('bdata has to be booleans')
        if bcolors is None:
            bcolors = ['black', 'green']
        if bmeans is None:
            bmeans = ['True', 'False']
        if len(bcolors) < 2:
            raise ValueError('At least to border colors need to be provided')
        if len(bmeans) < 2:
            raise ValueError('At least to binary names need to be provided')
        b = bdata.values.flatten()
        b = [bcolors[0] if _ else bcolors[1] for _ in b]
        kwargs['edgecolor'] = b

    # PLOT
    mesh = ax.scatter(**kwargs)

    ## Make Color Bar
    #if cdata is not None:
    #    plt.colorbar(mesh, fraction=0.046, pad=0.04)
#
    ## Make Size Legend
    #slegend = None
    #if sdata is not None:
    #    poslab = 1.35 if cdata is not None else 1
    #    p = sdata.values.flatten()
    #    medv = ((max(p) - min(p)) / 2) + min(p)
    #    topl = "{:.2f}".format(max(sdata.values.flatten()))
    #    botl = "{:.2f}".format(min(sdata.values.flatten()))
    #    medl = "{:.2f}".format(medv)
    #    medv = ((medv - np.min(p)) / np.ptp(p)) * (top - bot) + bot
#
    #    legend_elements = [
    #        Line2D([0], [0], marker='o', color='w', label=topl,
    #               markeredgecolor='black', markersize=math.sqrt(top)),
    #        Line2D([0], [0], marker='o', color='w', label=medl,
    #               markeredgecolor='black', markersize=math.sqrt(medv)),
    #        Line2D([0], [0], marker='o', color='w', label=botl,
    #               markeredgecolor='black', markersize=math.sqrt(bot)),
    #    ]
#
    #    slegend = ax.legend(handles=legend_elements, labelspacing=5.5,
    #                        handletextpad=2, borderpad=2,
    #                        bbox_to_anchor=(poslab, 1.015))
#
    ## Make Border Legend
    #if bdata is not None:
    #    poslab = 1.35 if cdata is not None else 1
    #    legend_elements = [
    #        Line2D([0], [0], marker='o', color='w', label=bmeans[0],
    #               markeredgecolor=bcolors[0], markersize=math.sqrt(top)),
    #        Line2D([0], [0], marker='o', color='w', label=bmeans[1],
    #               markeredgecolor=bcolors[1], markersize=math.sqrt(top))
    #    ]
#
    #    ax.legend(handles=legend_elements, labelspacing=5.5,
    #              handletextpad=2, borderpad=2,
    #              bbox_to_anchor=(poslab, 0.32))
    #    if slegend is not None:
    #        ax.add_artist(slegend)

    # Image aspects
    ax.grid(False)
    ax.set_xticks(range(1, 13))
    ax.xaxis.tick_top()
    ax.set_yticks(range(1, 9))
#    ax.set_yticklabels(string.ascii_uppercase[0:9])
    ax.set_ylim((8.5, 0.48))
    ax.set_aspect(1)
    ax.tick_params(axis='both', which='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    
    ## adding text to 96 well plate
    for text, x, y in zip(cdata.values.flatten(), list(range(1, 13)) * 8, sorted(list(range(1, 9)) * 12)):
        if text == 0: 
            new_text = ''
        else:
            new_text = text 
        if text < 10:
            plt.text(x-0.175, y+0.15, str(new_text), fontsize = 23)
        else:
            plt.text(x-0.275, y+0.15, str(new_text), fontsize = 23)

    return fig, ax

def build_next_history(starting_acceptor_size, ubi_donor_list, next_multimer_size, enzyme_of_deprot_first):
    previous_multimer_size = next_multimer_size-1
    acceptor_history_df = pd.read_csv('.data/core_data/' + str(starting_acceptor_size) + "mer__to_"  + str(previous_multimer_size) + 'mer_acceptor_history.csv')
    reaction_history_df = pd.read_csv('.data/core_data/' + str(starting_acceptor_size) + "mer__to_" + str(previous_multimer_size) + 'mer_final_reaction_history.csv')
    monomer_history_df = pd.read_csv('.data/core_data/' + str(starting_acceptor_size) + "mer__to_" + str(previous_multimer_size) + 'mer_monomer_history.csv')
    
    #column_names = list(acceptor_history_df.columns)
    #for index, name in enumerate(column_names): 
    #    ##create functions that removes reactions if nothing happens with a smac reaction
    #    if index % 2 == 0:
    #        acceptor_history_df[name] = acceptor_history_df.dtypes.apply(lambda x: x.name).to_dict()
    #    ### ### create functions that removes reactions if nothing happens....
    #    elif index % 2 == 1:
    #        acceptor_history_df[name] = acceptor_history_df.dtypes.apply(lambda x: x.name).to_dict()
    #        monomer_history_df[name] = monomer_history_df.dtypes.apply(lambda x: x.name).to_dict()

    acceptor_history_list = acceptor_history_df.values.tolist()
    reaction_history_list = reaction_history_df.values.tolist()
    monomer_history_list = monomer_history_df.values.tolist()
    create_synthesis_dataframes(acceptor_history_list,reaction_history_list,monomer_history_list, ubi_donor_list, enzyme_of_deprot_first)
    return 


# building xlsx files
def building_xlsx_files(monomer_or_enzyme_df, enzyme_mix_monomer_df, acceptor_count_df, count_deprot_df, stock_amounts_dict, final_amounts_dict, file_name, start_row): 
    ## acceptor conc. = 0.05 mmol/L
    E1_info = {
        'enzymes or monomer': ['hUba1'],
        'Molecular_Weight [Da]': [116300],  #  molecular weights
        # 1uM
        'Reaction Concentration [mmol/L]' : [final_amounts_dict['hUba1']/1000]
        }


    ## 
    E2_info = {
        'enzymes or monomer': ['Ube13/Mms2', 'gp78/Ube2g2', 'Ube2K'],
        'Molecular_Weight [Da]': [33800, 27600, 22400],  #  molecular weights
        # 20uM
        'Reaction Concentration [mmol/L]' : [final_amounts_dict['Ube13/Mms2']/1000,
                                            final_amounts_dict['gp78/Ube2g2']/1000,
                                            final_amounts_dict['Ube2K']/1000
                                            ]}

    monomer_info = {
        'enzymes or monomer': ['ubi_ubq_1_K48_ABOC_K63_ABOC', 'ubi_ubq_1_K48_SMAC_K63_ABOC',
        'ubi_ubq_1_K48_ABOC_K63_SMAC', 'ubi_ubq_1_K48_SMAC',
        'ubi_ubq_1_K63_SMAC'],
        'Molecular_Weight [Da]': [8864.13, 8864.13, 8733.00, 8733.00, 8733.00],
        'Reaction Concentration [mmol/L]' : [final_amounts_dict['ubi_ubq_1_K48_ABOC_K63_ABOC']/1000,
                                            final_amounts_dict['ubi_ubq_1_K48_SMAC_K63_ABOC']/1000,
                                            final_amounts_dict['ubi_ubq_1_K48_ABOC_K63_SMAC']/1000,
                                            final_amounts_dict['ubi_ubq_1_K48_SMAC']/1000,
                                            final_amounts_dict['ubi_ubq_1_K63_SMAC']/1000
                                            ]}  #  molecular weights

    E1_info_df = pd.DataFrame(E1_info)
    E2_info_df = pd.DataFrame(E2_info)
    monomer_info_df  = pd.DataFrame(monomer_info)

    ## reaction volume
    monomer_info_df['Reaction Volume [uL]'] = final_amounts_dict['Reaction Volume [uL]']
    E2_info_df['Reaction Volume [uL]'] = final_amounts_dict['Reaction Volume [uL]']
    E1_info_df['Reaction Volume [uL]'] = final_amounts_dict['Reaction Volume [uL]']

    monomer_info_df = monomer_info_df.merge(monomer_or_enzyme_df, on= 'enzymes or monomer')
    E2_info_df = E2_info_df.merge(monomer_or_enzyme_df, on= 'enzymes or monomer')
    E1_info_df['count'] = int(E2_info_df['count'].sum())

    ## get the starting concs.
    all_info_df = pd.concat([E1_info_df, E2_info_df, monomer_info_df], axis = 0).reset_index().drop('index', axis=1)
    ## get the reconstitution volume which is 3 times the final conc. 
    all_info_df['Stock Concentration [mmol/L]'] = 0
    for i in list(all_info_df['enzymes or monomer'].unique()):
        index = all_info_df[all_info_df['enzymes or monomer'] == i].index[0]
        all_info_df.at[index, 'Stock Concentration [mmol/L]'] = stock_amounts_dict[i]/1000

    ## Total values
    all_info_df['Total Mols [mmol]'] = all_info_df['Reaction Concentration [mmol/L]']* (all_info_df['Reaction Volume [uL]'] /1000000) * all_info_df['count']
    all_info_df['Total Mass [mg]'] = all_info_df['Total Mols [mmol]'] * all_info_df['Molecular_Weight [Da]'] 
    all_info_df['Total Stock Volume Needed [uL]'] = (all_info_df['Total Mols [mmol]']/ all_info_df['Stock Concentration [mmol/L]'])*1000000

    ## add in stock conc. 
    ## 
    ## 
    enzyme_mix_monomer_df[['Monomer', 'Enzymes']] = enzyme_mix_monomer_df['enzymes + monomer'].str.split('+', expand=True)
    enzyme_mix_monomer_df['Enzymes'] = enzyme_mix_monomer_df['Enzymes'].str.replace(' ', '')
    enzyme_mix_monomer_df['Monomer'] = enzyme_mix_monomer_df['Monomer'].str.replace(' ', '')



    enzyme_mix_monomer_df['Enzymes Stock Concentration [mmol/L]'] = enzyme_mix_monomer_df['Enzymes'].map(stock_amounts_dict)   
    enzyme_mix_monomer_df['Enzymes Stock Concentration [mmol/L]'] = enzyme_mix_monomer_df['Enzymes Stock Concentration [mmol/L]']/1000
    enzyme_mix_monomer_df['Enzymes Final Concentration [mmol/L]'] = enzyme_mix_monomer_df['Enzymes'].map(final_amounts_dict)   
    enzyme_mix_monomer_df['Enzymes Final Concentration [mmol/L]'] = enzyme_mix_monomer_df['Enzymes Final Concentration [mmol/L]']/1000

    enzyme_mix_monomer_df['Monomer Stock Concentration [mmol/L]'] = enzyme_mix_monomer_df['Monomer'].map(stock_amounts_dict)   
    enzyme_mix_monomer_df['Monomer Stock Concentration [mmol/L]'] = enzyme_mix_monomer_df['Monomer Stock Concentration [mmol/L]']/1000
    enzyme_mix_monomer_df['Monomer Final Concentration [mmol/L]'] = enzyme_mix_monomer_df['Monomer'].map(final_amounts_dict)   
    enzyme_mix_monomer_df['Monomer Final Concentration [mmol/L]'] = enzyme_mix_monomer_df['Monomer Final Concentration [mmol/L]']/1000

    enzyme_mix_monomer_df['hUba1 Stock Concentration [mmol/L]'] = stock_amounts_dict['hUba1']
    enzyme_mix_monomer_df['hUba1 Final Concentration [mmol/L]'] = final_amounts_dict['hUba1'] 


    enzyme_mix_monomer_df['Reaction Volume [uL]'] = final_amounts_dict['Reaction Volume [uL]']
    enzyme_mix_monomer_df['Total Volume [uL]'] = enzyme_mix_monomer_df['count']*enzyme_mix_monomer_df['Reaction Volume [uL]']

    enzyme_mix_monomer_df['Enzymes Stock Volume [uL]'] = (enzyme_mix_monomer_df['Total Volume [uL]'] * enzyme_mix_monomer_df['Enzymes Final Concentration [mmol/L]'])/ enzyme_mix_monomer_df['Enzymes Stock Concentration [mmol/L]']
    enzyme_mix_monomer_df['Monomer Stock Volume [uL]'] = (enzyme_mix_monomer_df['Total Volume [uL]'] * enzyme_mix_monomer_df['Monomer Final Concentration [mmol/L]'])/ enzyme_mix_monomer_df['Monomer Stock Concentration [mmol/L]']
    enzyme_mix_monomer_df['hUba1 Stock Volume [uL]'] = (enzyme_mix_monomer_df['Total Volume [uL]'] * enzyme_mix_monomer_df['hUba1 Final Concentration [mmol/L]'])/ enzyme_mix_monomer_df['hUba1 Stock Concentration [mmol/L]']

    enzyme_mix_monomer_df['Buffer Volume [uL]'] = enzyme_mix_monomer_df['Total Volume [uL]'] - enzyme_mix_monomer_df['Enzymes Stock Volume [uL]'] - enzyme_mix_monomer_df['Monomer Stock Volume [uL]']- enzyme_mix_monomer_df['hUba1 Stock Volume [uL]']

    #enzyme_mix_monomer_df['Enzymes Volume [uL]']
    #enzyme_mix_monomer_df['Monomer Volume [uL]']
    #enzyme_mix_monomer_df['hUba1 Volume [uL]']

    enzyme_mix_monomer_df = enzyme_mix_monomer_df[["enzymes + monomer", 
                                                "E + M_code",
                                                "count",
                                                "Reaction Volume [uL]",
                                                "Total Volume [uL]",
                                                "Monomer",
                                                "Monomer Final Concentration [mmol/L]",
                                                "Monomer Stock Concentration [mmol/L]",
                                                "Monomer Stock Volume [uL]","Enzymes",
                                                "Enzymes Final Concentration [mmol/L]",
                                                "Enzymes Stock Concentration [mmol/L]",
                                                "Enzymes Stock Volume [uL]","hUba1 Final Concentration [mmol/L]",
                                                "hUba1 Stock Concentration [mmol/L]","hUba1 Stock Volume [uL]",
                                                "Buffer Volume [uL]"]]

    acceptor_file_path = pd.read_csv('/Users/ekummelstedt/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/aa_looping_through_builders/.data/data_for_analysis/acceptor_df.csv')
    acceptor_info_df = acceptor_count_df.merge(acceptor_file_path, on= 'acceptor_string')

    acceptor_info_df['Final Concentration [mmol/L]'] = final_amounts_dict['Acceptor Concentration']/1000 #  volumes in µL
    acceptor_info_df['Molecular_Weight [Da]'] = 18522.07
    acceptor_info_df['Total Volume [uL]'] = acceptor_info_df['count']*final_amounts_dict['Reaction Volume [uL]']
    acceptor_info_df['Total Mols [mmol]'] = acceptor_info_df['Final Concentration [mmol/L]']* (acceptor_info_df['Total Volume [uL]'] /1000000)
    acceptor_info_df['Total Mass [mg]'] = acceptor_info_df['Total Mols [mmol]'] * acceptor_info_df['Molecular_Weight [Da]']

    acceptor_info_df = acceptor_info_df[['acceptor_string', 
                                         'acceptor_number', 
                                         'starting_acceptor', 
                                         'enzyme', 
                                         'monomer',
                                         'count',
                                         'Final Concentration [mmol/L]', 
                                         'Molecular_Weight [Da]',
                                         'Total Volume [uL]', 
                                         'Total Mols [mmol]', 
                                         'Total Mass [mg]']]

    count_deprot_df['Total Volume [mL]'] = 0
    count_deprot_df['Reaction Volume [uL]'] = 200

    smac_index = count_deprot_df[count_deprot_df['deprot'] == 'SMAC_deprot'].index[0]
    fake_index = count_deprot_df[count_deprot_df['deprot'] == 'Fake_Wash'].index[0]

    count_deprot_df.at[smac_index, 'Total Volume [mL]'] =  (count_deprot_df.loc[smac_index, 'count'] *2) * count_deprot_df.loc[smac_index, 'Reaction Volume [uL]']/1000
    count_deprot_df.at[fake_index, 'Total Volume [mL]'] = ((count_deprot_df.loc[fake_index, 'count'] *2) +96) * count_deprot_df.loc[fake_index, 'Reaction Volume [uL]']/1000

    count_deprot_df['HEPES Conc. [mmol/L]'] = 0
    count_deprot_df['NaCl Conc. [mmol/L]'] = 0
    count_deprot_df['Imidazole Conc. [mmol/L]'] = 0
    count_deprot_df['PLP Conc. [mmol/L]'] = 0

    count_deprot_df.at[smac_index, 'HEPES Conc. [mmol/L]'] = 150
    count_deprot_df.at[fake_index, 'HEPES Conc. [mmol/L]'] = 100
    count_deprot_df.at[smac_index, 'NaCl Conc. [mmol/L]'] = 150
    count_deprot_df.at[fake_index, 'NaCl Conc. [mmol/L]'] = 500
    count_deprot_df.at[smac_index, 'Imidazole Conc. [mmol/L]'] = 0
    count_deprot_df.at[fake_index, 'Imidazole Conc. [mmol/L]'] = 30
    count_deprot_df.at[smac_index, 'PLP Conc. [mmol/L]'] = 20
    count_deprot_df.at[fake_index, 'PLP Conc. [mmol/L]'] = 0

    count_deprot_df['HEPES MW [g/mol]'] = 238.3012
    count_deprot_df['NaCl MW [g/mol]'] = 58.44
    count_deprot_df['Imidazole MW [g/mol]'] = 68.077
    count_deprot_df['PLP MW [g/mol]'] = 247.142

    count_deprot_df['HEPES mass [mg]'] = count_deprot_df['HEPES Conc. [mmol/L]']*(count_deprot_df['Total Volume [mL]']/1000)*count_deprot_df['HEPES MW [g/mol]']
    count_deprot_df['NaCl mass [mg]'] = count_deprot_df['NaCl Conc. [mmol/L]']*(count_deprot_df['Total Volume [mL]']/1000)*count_deprot_df['NaCl MW [g/mol]']
    count_deprot_df['Imidazole mass [mg]'] = count_deprot_df['Imidazole Conc. [mmol/L]']*(count_deprot_df['Total Volume [mL]']/1000)*count_deprot_df['Imidazole MW [g/mol]']
    count_deprot_df['PLP mass [mg]'] = count_deprot_df['PLP Conc. [mmol/L]']*(count_deprot_df['Total Volume [mL]']/1000)*count_deprot_df['PLP MW [g/mol]']

    count_deprot_df.at[smac_index, 'Total Volume [mL]'] =  (count_deprot_df.loc[smac_index, 'count'] *2) * count_deprot_df.loc[smac_index, 'Reaction Volume [uL]']/1000

    count_deprot_df = count_deprot_df[['deprot', 
                                       'deprot_code',
                                        'count', 
                                        'Total Volume [mL]', 
                                        'Reaction Volume [uL]',
                                        'HEPES Conc. [mmol/L]', 
                                        'HEPES MW [g/mol]',
                                        'HEPES mass [mg]', 
                                        'NaCl Conc. [mmol/L]',
                                        'NaCl MW [g/mol]', 
                                        'NaCl mass [mg]', 
                                        'Imidazole Conc. [mmol/L]', 
                                        'Imidazole MW [g/mol]', 
                                        'Imidazole mass [mg]',
                                        'PLP Conc. [mmol/L]', 
                                        'PLP MW [g/mol]',
                                        'PLP mass [mg]']]
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(".data/excel_cal_concs/" + file_name + str(start_row +1) + " to " + str(start_row+ 14) + ".xlsx", engine="xlsxwriter")

    ## sheet 1 is for just final cal amounts...
    ## sheet 2 provides all the  calculations... maybe add in later...
    ## have both sheets and heat maps solved by end of day...

    ## what do we want
    # how much buffer we need...
    ## SHEET 1 = "Full_Calculations"
    #1 dimer amounts
    acceptor_info_df.to_excel(writer, sheet_name="Full_Calculations",startrow=2 , startcol=0)
    #2 individual amounts of everything
    all_info_df.to_excel(writer, sheet_name="Full_Calculations",startrow=(len(acceptor_info_df)+ 5), startcol=0)
    #3 mix amounts
    enzyme_mix_monomer_df.to_excel(writer, sheet_name="Full_Calculations",startrow=(len(acceptor_info_df)+ 5 + len(all_info_df) + 4) , startcol=0)
    #4 smac or deprot amounts
    count_deprot_df.to_excel(writer, sheet_name="Full_Calculations",startrow=(len(enzyme_mix_monomer_df) + 5 + len(acceptor_info_df)+ 4 + len(all_info_df) + 4) , startcol=0)

    full_calc_sheets = writer.sheets['Full_Calculations']
    full_calc_sheets.write(0,0,(file_name + "  " + str(start_row +1) + " to " + str(start_row+ 14)))
    full_calc_sheets.write(1,0,'Acceptor Calculations')
    full_calc_sheets.write((len(acceptor_info_df)+ 4),0,'Monomer and Enyzme Calculations')
    full_calc_sheets.write((len(acceptor_info_df)+ 5 + len(all_info_df) + 3),0,'Monomer/Enzyme Mix Calculations')
    full_calc_sheets.write((len(enzyme_mix_monomer_df) + 5 + len(acceptor_info_df)+ 4 + len(all_info_df) + 3),0,'Smac/Fake deprot Calculations')

    # Close the Pandas Excel writer and output the Excel file
    writer.close()


def experimental_output_design(file_name, enzyme_of_deprot_first, start_row):
    
    file_path = '/Users/ekummelstedt/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/aa_looping_through_builders/.data/reaction_history/' + file_name + '.csv'
    ## fix the numbering at the end
    df = pd.read_csv(file_path)
    columns = df.columns[4:]
    ## remove first 4 columns
    df_reaction_pathway = df[columns]
    ## get the odd number of columns from columns
    column_zero = []
    column_odds = []
    column_evens = []
    for i in columns.map(int):
        if i == 0: 
            column_zero.append(str(i))
        elif i % 2 == 1:
            column_odds.append(str(i))
        elif i% 2 ==0:
            column_evens.append(str(i))

    ## just pull out the columns with enzyme reactions
    df_starting_acceptor = df_reaction_pathway[column_zero]
    if enzyme_of_deprot_first == 'enzyme': 
        enyzme_columns = column_odds
        deprot_columns = column_evens
    if enzyme_of_deprot_first == 'deprot': 
        enyzme_columns = column_evens
        deprot_columns = column_odds

    df_enzymes = df_reaction_pathway[enyzme_columns]
    df_deprot = df_reaction_pathway[deprot_columns]

    ## take 14 reactions 
    start_row_edit = start_row*2
    start_row_14 = start_row_edit+(14*2)



    df_enzymes = df_enzymes.astype(str)
    df_enzymes_14 = df_enzymes[start_row_edit:start_row_14]
    df_starting_acceptor = df_starting_acceptor.astype(str)
    df_starting_acceptor_14 = df_starting_acceptor[start_row_edit:start_row_14]
    df_deprot = df_deprot.astype(str)
    df_deprot_14 = df_deprot[start_row_edit:start_row_14]

    ##############################################################################
    ############################ ENZYMES + MONOMERS ##############################
    #############################################################################

    ## combined every two rows... 
    ## add a column with the same number every two columns
    # creating groupby list
    b1 = list(range(len(df_enzymes_14)))
    groupby_list = [int(x / 2) for x in b1]
    df_enzymes_14['combining_column'] = groupby_list
    ## use groupby 
    ## no other way to deal with the below...
    if enzyme_of_deprot_first == 'enzyme': 
        if len(enyzme_columns) == 1:
            result = df_enzymes_14.groupby('combining_column').agg({'1': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 2:
            result = df_enzymes_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 3:
            result = df_enzymes_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 4:
            result = df_enzymes_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join,
                                                                    '7': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 5:
            result = df_enzymes_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join,
                                                                    '7': ' + '.join,
                                                                    '9': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 6:
            result = df_enzymes_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join,
                                                                    '7': ' + '.join,
                                                                    '9': ' + '.join,
                                                                    '11': ' + '.join}).reset_index().drop('combining_column', axis=1)

    elif enzyme_of_deprot_first == 'deprot': 
        if len(enyzme_columns) == 1:
            result = df_enzymes_14.groupby('combining_column').agg({'2': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 2:
            result = df_enzymes_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 3:
            result = df_enzymes_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join,
                                                                    '6': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 4:
            result = df_enzymes_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join,
                                                                    '6': ' + '.join,
                                                                    '8': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(enyzme_columns) == 5:
            result = df_enzymes_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join,
                                                                    '6': ' + '.join,
                                                                    '8': ' + '.join,
                                                                    '10': ' + '.join}).reset_index().drop('combining_column', axis=1)


    ## start of library building..
    # have a library already built
    list_of_donors = ['ubi_ubq_1_K48_SMAC', 
                    'ubi_ubq_1_K63_SMAC', 
                    'ubi_ubq_1_K48_SMAC_K63_ABOC', 
                    'ubi_ubq_1_K48_ABOC_K63_SMAC', 
                    'ubi_ubq_1_K48_ABOC_K63_ABOC']

    list_of_enzymes = ['gp78/Ube2g2', 
                        'Ube2K', 
                        'Ube13/Mms2']

    
    

    ## build encoded library.. all combinations of donors + enzymes
    encoded_dictionary= {}
    count = 1 
    for i in list_of_enzymes: 
        for j in list_of_donors:
            ## these are the self reactions
            if ((j == 'ubi_ubq_1_K48_SMAC') & (i=='Ube13/Mms2')) | ((j == 'ubi_ubq_1_K63_SMAC') & (i=='Ube2K'))|((j == 'ubi_ubq_1_K63_SMAC') & (i=='gp78/Ube2g2')):
                count = count
            else:
                encoded_dictionary[j + ' + ' + i] = count 
                count = count+1 
    
    
    ### end of library building... 
    # listify the enzymes 
    list_of_enzymes = []
    for i in enyzme_columns: 
        column_as_list = result[i].tolist()
        list_of_enzymes = list_of_enzymes + column_as_list

    ## change Ube13/Mms2_branching to Ube13/Mms2
    list_of_enzymes_cleaned = list(map(lambda x: x.replace('Ube13/Mms2_branching', 'Ube13/Mms2'), list_of_enzymes))

    ## 
    count_dictionary = {}
    enzymes_monomers = encoded_dictionary.keys()
    for i in enzymes_monomers:
        count_dictionary[i] = list_of_enzymes_cleaned.count(i)
    ## number of times everything appears...
    well_plate_96 = pd.DataFrame()
    for index, i in enumerate(enyzme_columns):
        if index == 0: 
            well_plate_96_i = result[[i]].iloc[0:8]
            well_plate_96_ii = result[[i]].iloc[8:14].reset_index().drop(['index'], axis=1)
            well_plate_96 = pd.concat([well_plate_96_i, well_plate_96_ii], axis=1, ignore_index=True) 
        else:
            well_plate_96_i = result[[i]].iloc[0:8]
            well_plate_96_ii = result[[i]].iloc[8:14].reset_index().drop(['index'], axis=1)
            well_plate_96 = pd.concat([well_plate_96, well_plate_96_i, well_plate_96_ii], axis=1, ignore_index=True) 

    for i in well_plate_96.columns:
        well_plate_96[i] = well_plate_96[i].astype(str)
        ## make sure branching is Ube13/Mms2
        well_plate_96[i] = well_plate_96[i].map(lambda x: x.replace('Ube13/Mms2_branching', 'Ube13/Mms2'))
        well_plate_96[i] = well_plate_96[i].map(encoded_dictionary)

    ## any well on the 96 well plate assign 0
    fillna_value = min(encoded_dictionary.values())- 1
    well_plate_96 = well_plate_96.fillna(fillna_value).astype(int)

    for i in range(12):
        if (i in well_plate_96.columns)==False:
            well_plate_96[i] = fillna_value

    figure_name1 = (file_name + ' E2 + monomer ' + str(start_row +1) + ' to ' + str(start_row+14))
    donor_enzyme_fig, donor_enzyme_ax = plot_96wells(cdata = well_plate_96, figure= 1, figure_name = figure_name1, colorbar_type= 'PuRd')
    donor_enzyme_fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

    encoded_df = pd.DataFrame(encoded_dictionary.items(), columns=['enzymes + monomer', 'E + M_code'])
    count_df = pd.DataFrame(count_dictionary.items(), columns=['enzymes + monomer', 'count'])

    enzyme_mix_monomer_df = encoded_df.merge(count_df)

    reordered_df = enzyme_mix_monomer_df.copy()
    reordered_df[['monomer', 'enzymes']] = reordered_df['enzymes + monomer'].str.split('+', expand=True)
    reordered_df['enzymes'] = reordered_df['enzymes'].str.replace(' ', '', regex=True)
    reordered_df['monomer'] = reordered_df['monomer'].str.replace(' ', '', regex=True)
    reordered_df = reordered_df[['enzymes', 'monomer', 'E + M_code', 'count']]
    enzyme_list = reordered_df['enzymes'].unique().tolist()
    monomer_list = reordered_df['monomer'].unique().tolist()
    
    enzyme_monomer_count_dict = {}
    for i in enzyme_list:
        enzyme_monomer_count_dict[i] = int(reordered_df[reordered_df['enzymes']==i]['count'].sum())
    for i in monomer_list:
        enzyme_monomer_count_dict[i] = int(reordered_df[reordered_df['monomer']==i]['count'].sum())
    monomer_or_enzyme_df = pd.DataFrame(enzyme_monomer_count_dict.items(), columns=['enzymes or monomer', 'count'])
    
    ## table and 96 well map of dimers 
    df_starting_acceptor_14['combining_column'] = groupby_list
    acceptors_df = df_starting_acceptor_14.groupby('combining_column').agg({'0': ' + '.join}).reset_index().drop('combining_column', axis=1)
    acceptors_df['0'] = acceptors_df['0'].str.replace(' + nan', '')
    
    ##############################################################################
    ################################## ACCEPTORS #################################
    ##############################################################################
    acceptor_file_path = '/Users/ekummelstedt/le_code_base/ubiquitin_syn/Ubiquitin_Scripts/aa_looping_through_builders/.data/data_for_analysis/acceptor_df.csv'
    acceptors_map_df = pd.read_csv(acceptor_file_path)
    acceptor_mapping_dict = dict(zip(acceptors_map_df['acceptor_string'],acceptors_map_df['acceptor_number']))
    numbered_acceptors_df = pd.DataFrame()
    numbered_acceptors_df['0'] = acceptors_df['0'].map(acceptor_mapping_dict)
    ### END OF CHANGE 
    
    numbered_acceptors_df = numbered_acceptors_df.astype(int)
    acceptor_well_plate_96_i = numbered_acceptors_df[['0']].iloc[0:8]
    acceptor_well_plate_96_ii = numbered_acceptors_df[['0']].iloc[8:14].reset_index().drop(['index'], axis=1)
    acceptor_well_plate_96 = pd.concat([acceptor_well_plate_96_i, acceptor_well_plate_96_ii], axis=1, ignore_index=True) 

    ## any well on the 96 well plate assign 0
    fillna_value = 0
    acceptor_well_plate_96 = acceptor_well_plate_96.fillna(fillna_value).astype(int)

    for i in range(12):
        if (i in acceptor_well_plate_96.columns)==False:
            acceptor_well_plate_96[i] = fillna_value
    
    figure_name2 = (file_name + ' acceptors ' + str(start_row +1) + ' to ' + str(start_row+14))
    acceptor_fig, acceptor_ax = plot_96wells(figure=2, figure_name = figure_name2, colorbar_type= 'Blues', cdata = acceptor_well_plate_96)
    acceptor_fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

    acceptor_count_dictionary = {'acceptor_string': [], 'count': []}
    for i in list(acceptors_df['0'].unique()):
        acceptor_count_dictionary['acceptor_string'].append(i)
        acceptor_count_dictionary['count'].append(acceptors_df['0'].tolist().count(i))
    acceptor_count_df = pd.DataFrame.from_dict(acceptor_count_dictionary)
    
    #############################################################################
    ############################ DEPROTECTIONS ##################################
    #############################################################################

    ## combined every two rows... 
    ## add a column with the same number every two columns
    # creating groupby list
    b1 = list(range(len(df_deprot_14)))
    groupby_list = [int(x / 2) for x in b1]
    df_deprot_14['combining_column'] = groupby_list
    ## use groupby 
    ## no other way to deal with the below...
    if enzyme_of_deprot_first == 'deprot': 
        if len(deprot_columns) == 1:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'1': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 2:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 3:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 4:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join,
                                                                    '7': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 5:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                    '3': ' + '.join,
                                                                    '5': ' + '.join,
                                                                    '7': ' + '.join,
                                                                    '9': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 6:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'1': ' + '.join,
                                                                            '3': ' + '.join,
                                                                            '5': ' + '.join,
                                                                            '7': ' + '.join,
                                                                            '9': ' + '.join,
                                                                            '11': ' + '.join}).reset_index().drop('combining_column', axis=1)

    elif enzyme_of_deprot_first == 'enzyme': 
        if len(deprot_columns) == 1:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'2': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 2:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 3:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join,
                                                                    '6': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 4:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join,
                                                                    '6': ' + '.join,
                                                                    '8': ' + '.join}).reset_index().drop('combining_column', axis=1)
        elif len(deprot_columns) == 5:
            result_deprot = df_deprot_14.groupby('combining_column').agg({'2': ' + '.join,
                                                                    '4': ' + '.join,
                                                                    '6': ' + '.join,
                                                                    '8': ' + '.join,
                                                                    '10': ' + '.join}).reset_index().drop('combining_column', axis=1)


    for i in result_deprot.columns:
        result_deprot[i] = result_deprot[i].str.replace('nan + ', '')

    ## build encoded library.. all combinations of donors + enzymes
    encoded_deprot_dictionary= {'SMAC_deprot' : 1, 'Fake_Wash': 2}
    
    ### end of library building... 
    # listify the enzymes 
    list_of_deprots = []
    for i in deprot_columns: 
        column_as_list = result_deprot[i].tolist()
        list_of_deprots = list_of_deprots + column_as_list

    ## 
    count_deprot_dictionary = {}
    deprots = encoded_deprot_dictionary.keys()
    for i in deprots:
        count_deprot_dictionary[i] = list_of_deprots.count(i)
    ## number of times everything appears...
    deprot_well_plate_96 = pd.DataFrame()
    for index, i in enumerate(deprot_columns):
        if index == 0: 
            deprot_well_plate_96_i = result_deprot[[i]].iloc[0:8]
            deprot_well_plate_96_ii = result_deprot[[i]].iloc[8:14].reset_index().drop(['index'], axis=1)
            deprot_well_plate_96 = pd.concat([deprot_well_plate_96_i, deprot_well_plate_96_ii], axis=1, ignore_index=True) 
        else:
            deprot_well_plate_96_i = result_deprot[[i]].iloc[0:8]
            deprot_well_plate_96_ii = result_deprot[[i]].iloc[8:14].reset_index().drop(['index'], axis=1)
            deprot_well_plate_96 = pd.concat([deprot_well_plate_96, deprot_well_plate_96_i, deprot_well_plate_96_ii], axis=1, ignore_index=True) 

    for i in deprot_well_plate_96.columns:
        deprot_well_plate_96[i] = deprot_well_plate_96[i].astype(str)
        ## make sure branching is Ube13/Mms2
        deprot_well_plate_96[i] = deprot_well_plate_96[i].map(encoded_deprot_dictionary)

    ## any well on the 96 well plate assign 0
    fillna_value = 0
    deprot_well_plate_96 = deprot_well_plate_96.fillna(fillna_value).astype(int)

    for i in range(12):
        if (i in deprot_well_plate_96.columns)==False:
            deprot_well_plate_96[i] = fillna_value
    
    figure_name3 = (file_name + ' deprotections ' + str(start_row +1) + ' to ' + str(start_row+14))
    deprot_fig, deprot_ax = plot_96wells(cdata = deprot_well_plate_96, figure= 3, figure_name = figure_name3, colorbar_type= 'BuGn')
    deprot_fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

    final_count_deprot_dictionary = {'deprot': list(count_deprot_dictionary.keys()), 'count': list(count_deprot_dictionary.values()), 'deprot_code': list(encoded_deprot_dictionary.values())}
    count_deprot_df = pd.DataFrame.from_dict(final_count_deprot_dictionary)

    donor_enzyme_fig.savefig(('.data/plate_map_images/' + figure_name1 + '.png'))
    acceptor_fig.savefig(('.data/plate_map_images/' + figure_name2 + '.png'))
    deprot_fig.savefig(('.data/plate_map_images/' + figure_name3 + '.png'))

    ## add smac/fake deprot...
    return monomer_or_enzyme_df, enzyme_mix_monomer_df, acceptor_count_df, count_deprot_df, donor_enzyme_fig, acceptor_fig, deprot_fig