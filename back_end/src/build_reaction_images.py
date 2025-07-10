from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import re
import subprocess
import sys
import plotly.io as pio
import plotly

# Dynamically get the backend path relative to this script location
script_path = Path(__file__).resolve()
project_root = script_path.parents[2]  # Go up to project root (adjust if needed)
src_path = project_root / 'back_end' / 'src'
sys.path.insert(0, str(src_path))

from utils.utils import *
from utils.logging_utils import *
from main import *
"""
Functions for formatting Plotly traces and extracting NetworkX graphs from nested ubiquitin JSON.
"""

import networkx as nx
import plotly.graph_objects as go

# ============================================================

multimer_size = 4

# download CSV files
def download_data_dict(multimer_size):
    input_dir = project_root / 'back_end' / 'data' / 'filtered_reaction_database' / f'multimer_size_{multimer_size}'

    combined_database = pd.read_csv(input_dir / 'combined_database.csv', index_col=0)
    context_history = pd.read_csv(input_dir / 'context_history.csv', index_col=0)
    donor_history = pd.read_csv(input_dir / 'donor_history.csv', index_col=0)
    reaction_history = pd.read_csv(input_dir / 'reaction_history.csv', index_col=0)
    ubiquitin_history = pd.read_csv(input_dir / 'ubiquitin_history.csv', index_col=0)

    return {
        'combined_database': combined_database,
        'context_history': context_history,
        'donor_history': donor_history,
        'reaction_history': reaction_history,
        'ubiquitin_history': ubiquitin_history
    }

# Create plate dataframes for the selected indexes
data_dict = download_data_dict(multimer_size)
combined_database = data_dict['combined_database']
context_history = data_dict['context_history']
donor_history = data_dict['donor_history']
reaction_history = data_dict['reaction_history']
ubiquitin_history = data_dict['ubiquitin_history']

# Create the plate dataframes for the selected indexes
# Select the indexes of the multimer_ids that are used in synthesis
selected_ids = ["Ub4_5","Ub4_10","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_8", "Ub4_11","Ub4_7"]
indexes = list()

# Get the indexes of the selected multimer_ids that are used in synthesis
for id in selected_ids:
    new_index = int(combined_database[(combined_database['multimer_id'] == id) & (combined_database['used_in_synthesis'] == 1)]['index'].unique()[0])
    indexes.append(new_index)

# ============================================================


idx = indexes[0]
test_line = ubiquitin_history[ubiquitin_history['index'] == idx]


# ===================== Utility Functions =====================

def format_node_trace_legend(node_trace, legendgroup):
    """
    Format the node trace for consistent legend placement in Plotly.

    Args:
        node_trace (go.Scatter): The Plotly node trace to be formatted.
        legendgroup (str): The legend group label for this trace.

    Returns:
        go.Scatter: The formatted node trace with legendgroup and showlegend set.
    """
    node_trace.legendgroup = legendgroup
    node_trace.showlegend = True
    return node_trace


def annotate_graph_nodes(fig, pos, labels, font_size=12):
    """
    Add graph node labels as annotations to a Plotly figure.

    Args:
        fig (go.Figure): The Plotly figure to annotate.
        pos (dict): Mapping of node to (x, y) positions.
        labels (dict): Mapping of node to label text.
        font_size (int, optional): Font size for annotations.

    Returns:
        go.Figure: The updated Plotly figure with node annotations added.
    """
    for node, (x, y) in pos.items():
        fig.add_annotation(
            x=x, y=y,
            text=labels.get(node, str(node)),
            showarrow=False,
            font=dict(size=font_size)
        )
    return fig


def build_plotly_traces_from_graph(G, pos, node_color='blue', edge_color='gray'):
    """
    Construct Plotly node and edge scatter traces from a NetworkX graph and layout positions.

    Args:
        G (networkx.Graph): The input NetworkX graph.
        pos (dict): Mapping of node to (x, y) positions.
        node_color (str, optional): Color for nodes.
        edge_color (str, optional): Color for edges.

    Returns:
        tuple: (node_trace, edge_trace), both as go.Scatter objects.
    """
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1, color=edge_color),
        hoverinfo='none',
        mode='lines'
    )

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=False,
            color=node_color,
            size=10,
            line_width=2
        )
    )
    return node_trace, edge_trace


def extract_graph_from_ubiquitin_json(ub_json, node_list=None, edge_list=None, parent=None, edge_labels=None):
    """
    Recursively parse a nested ubiquitin JSON object into lists of nodes and edges for NetworkX.
    Also returns a list of edge labels (K48/K63 or None).
    """
    if node_list is None:
        node_list = []
    if edge_list is None:
        edge_list = []
    if edge_labels is None:
        edge_labels = []
    node_id = ub_json.get('id')
    node_list.append(node_id)
    for child in ub_json.get('children', []):
        child_id = child.get('id')
        if parent is not None:
            edge_list.append((parent, node_id))
            edge_labels.append(None)
        # If this child is a K48/K63, add an edge from this node to the child and label it
        if child_id and (child_id.endswith('K48') or child_id.endswith('K63')):
            edge_list.append((node_id, child_id))
            edge_labels.append(child_id[-3:])
        extract_graph_from_ubiquitin_json(child, node_list, edge_list, node_id, edge_labels)
    return node_list, edge_list, edge_labels


# ===================== Visualization Function =====================

def visualize_ubiquitin_tree(ubiquitin_json, title="Ubiquitin Tree Visualization"):
    import networkx as nx
    import plotly.graph_objects as go

    # Parse nodes, edges, and edge labels
    node_list, edge_list, edge_labels = extract_graph_from_ubiquitin_json(ubiquitin_json)

    # Create NetworkX graph
    G = nx.DiGraph()
    G.add_nodes_from(node_list)
    G.add_edges_from(edge_list)

    # Hierarchical layout
    def hierarchy_pos(G, root=None, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
        if root is None:
            root = list(nx.topological_sort(G))[0]
        def _hierarchy_pos(G, root, leftmost, width, vert_gap, vert_loc, xcenter, pos=None, parent=None):
            if pos is None:
                pos = {root: (xcenter, vert_loc)}
            else:
                pos[root] = (xcenter, vert_loc)
            children = list(G.successors(root))
            if len(children) != 0:
                dx = width / len(children)
                nextx = xcenter - width/2 - dx/2
                for child in children:
                    nextx += dx
                    pos = _hierarchy_pos(G, child, leftmost, width=dx, vert_gap=vert_gap,
                                         vert_loc=vert_loc+vert_gap, xcenter=nextx, pos=pos, parent=root)
            return pos
        return _hierarchy_pos(G, root, 0, width, vert_gap, vert_loc, xcenter)

    roots = [n for n, d in G.in_degree() if d == 0]
    root = roots[0] if roots else None
    pos = hierarchy_pos(G, root=root, width=1.0, vert_gap=0.2, vert_loc=0, xcenter=0.5)

    # Build edge traces, highlighting K48/K63 edges
    edge_x = []
    edge_y = []
    edge_colors = []
    for i, edge in enumerate(edge_list):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        color = 'red' if edge_labels[i] in ('K48', 'K63') else 'gray'
        edge_colors.append(color)
    # Plotly does not support per-segment color, so plot highlighted and normal edges separately
    edge_trace_highlight = go.Scatter(
        x=[edge_x[j] for j in range(0, len(edge_x), 3) if edge_colors[j//3]=='red' for k in range(3)],
        y=[edge_y[j] for j in range(0, len(edge_y), 3) if edge_colors[j//3]=='red' for k in range(3)],
        line=dict(width=2, color='red'),
        hoverinfo='none',
        mode='lines',
        showlegend=True,
        name='K48/K63 Edge'
    )
    edge_trace_normal = go.Scatter(
        x=[edge_x[j] for j in range(0, len(edge_x), 3) if edge_colors[j//3]!='red' for k in range(3)],
        y=[edge_y[j] for j in range(0, len(edge_y), 3) if edge_colors[j//3]!='red' for k in range(3)],
        line=dict(width=1, color='gray'),
        hoverinfo='none',
        mode='lines',
        showlegend=True,
        name='Other Edge'
    )

    # Build node trace (all nodes, no K48/K63 special nodes)
    node_x = [pos[node][0] for node in G.nodes() if node in pos]
    node_y = [pos[node][1] for node in G.nodes() if node in pos]
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=False,
            color='lightgray',
            size=10,
            line_width=2
        ),
        showlegend=False
    )

    # Create figure and add traces
    fig = go.Figure()
    fig.add_trace(edge_trace_normal)
    fig.add_trace(edge_trace_highlight)
    fig.add_trace(node_trace)

    # Annotate nodes (only those with positions)
    fig = annotate_graph_nodes(fig, pos, labels={node: node for node in node_list if node in pos})

    # Update layout
    fig.update_layout(
        title=title,
        title_x=0.5,
        showlegend=True,
        margin=dict(l=20, r=20, t=40, b=20),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='black',
        paper_bgcolor='black',
    )

    return fig


if __name__ == "__main__":
    # Example k48_dimer_ubiquitin structure
    k48_dimer_ubiquitin = {
        "protein": "1ubq-histag",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": {"protein": "1ubq",
                                                                        "chain_number": 2,
                                                                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                                                        "chain_length": 76,
                                                                        "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                                                                            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                                                                            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                                                                            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                                                                            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                                                                            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                                                                            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                                                                            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}]}},
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}
        ]
    }

    # Convert the branching_sites structure to the expected nested JSON format for visualization
    def convert_branching_sites_to_json(protein_dict, parent_id=""):
        # Recursively convert only K48 and K63 branching_sites to a nested children structure with unique node ids
        node_id = f"{protein_dict['protein']}_chain{protein_dict['chain_number']}" if not parent_id else f"{parent_id}/{protein_dict['protein']}_chain{protein_dict['chain_number']}"
        node = {
            'id': node_id,
            'children': []
        }
        for i, site in enumerate(protein_dict.get('branching_sites', [])):
            if site['site_name'] in ('K48', 'K63'):
                # Make a unique id for each site by including the parent node id and index
                site_id = f"{node_id}/{site['site_name']}_{i}"
                if isinstance(site.get('children'), dict):
                    child = convert_branching_sites_to_json(site['children'], parent_id=site_id)
                    child['id'] = site_id
                    node['children'].append(child)
                else:
                    node['children'].append({'id': site_id, 'children': []})
        return node

    ubiquitin_json = convert_branching_sites_to_json(k48_dimer_ubiquitin)
    fig = visualize_ubiquitin_tree(ubiquitin_json, title="K48 Dimer Ubiquitin Tree")
    # Save and open the image
    image_path = "/Users/ekummelstedt/Desktop/k48_dimer_ubiquitin_tree.png"
    fig.write_image(image_path, format="png")
    print(f"Saved ubiquitin tree image to: {image_path}")
    subprocess.run(["open", image_path])  # For macOS

