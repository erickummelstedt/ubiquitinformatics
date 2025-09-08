# Ubiquitinformatics

A Python and React application for ubiquitin chain simulation and visualization.

### 🐳 Option 1: Docker Setup (Recommended)
**No Python or Node.js installation needed!**
1. **Git** - For downloading the project
   - Check: `git --version` in Terminal
   - Install: [https://git-scm.com/download](https://git-scm.com/download)

2. **Docker Desktop** - Handles everything else automatically
   - Download: [https://www.docker.com/get-started](https://www.docker.com/get-started)
   - Install and start Docker Desktop
   - Ensure Docker is in your PATH: Run `docker --version` in Terminal to verify

```bash
# Step 1: Download the project
git clone https://github.com/erickummelstedt/ubiquitinformatics.git
cd ubiquitinformatics

# Step 2: Run everything with Docker (takes 2-3 minutes)
# Mac/Linux:
./start.sh

# Windows, working on developing a login page for windows
```

Access at: http://localhost:5173

### 🛠️ Option 2: Manual Setup (If you prefer local installation)
**- requires installing Python and Node.js:**

```bash
# Download the project
git clone https://github.com/erickummelstedt/ubiquitinformatics.git
cd ubiquitinformatics

# Run everything locally
# Mac/Linux:
./manual-start.sh
```

This will automatically set up Python environment, install dependencies, run simulation, and start the web interface.

Access at: http://localhost:5173

---

## 🧩 How Simulation Works: System Outline

Below is the structure and logic used for ubiquitin chain simulation and visualization in Ubiquitinformatics. We describe how the ubiG (ubiquitin graph) is defined as a mutable JSON object, which forms the basis for all simulation and visualization steps. The main functions are outlined here; all other functions in the codebase serve as utility functions supporting this core system.

## 📦 Ubiquitin Graph (ubiG) JSON Structure

Ubiquitinformatics represents ubiquitin chains as nested JSON objects ("ubiG"). Each ubiG describes a monomer or a branched chain, and can be arbitrarily nested to represent complex topologies. The JSON ubiG is mathematically a graph G(V, E), specifically a rooted tree, where each node (vertex) is a protein and each edge represents a conjugation (branching) event. These trees are characterized by a single origin point, unidirectional flow, and the capacity for controlled, multi-site branching. The JSON is converted into a Python dictionary for processing and simulation within the code.

The G(V,E) & JSON is numbered using preorder tree traversal, reflecting biological chain directionality and indexing each ubiquitin unit in the order: K63, K48, K33, K29, K27, K11, K6, M1. Therefore, every branch is fully numbered/traversed—regardless of its size—before moving to the next lysine in the current chain.

The base monomers—5 ubiquitin donors and 3 acceptors—are outlined in `building_blocks.py` and serve as the starting templates for chain construction and simulation.

**Basic Structure:**
```json
{
  "protein": "1ubq",                // Protein identifier
  "chain_number": 1,                  // Position in the chain/multimer
  "FASTA_sequence": "...",           // Amino acid sequence
  "chain_length": 76,                 // Number of residues
  "branching_sites": [                // List of possible branching sites
    {
      "site_name": "K63",            // Lysine or N-terminal site
      "sequence_id": "NIQ(K)EST",    // Sequence context
      "children": ""                  // Empty (no branch), label, or nested ubiG
    },
    // ... more sites ...
  ]
}
```

**Branching Example:**
A site can contain another ubiG as its `children`, allowing for recursive branching. For example, protein 1 with a K48-conjugated branch to protein 2:
```json
{
  "protein": "1ubq",
  "chain_number": 1,
  "FASTA_sequence": "...",
  "chain_length": 76,
  "branching_sites": [
    {
      "site_name": "K63",
      "sequence_id": "NIQ(K)EST",
      "children": ""
    },
    {
      "site_name": "K48",
      "sequence_id": "FAG(K)QLE",
      "children": {
        "protein": "1ubq",
        "chain_number": 2,
        "FASTA_sequence": "...",
        "chain_length": 76,
        "branching_sites": [
          {
            "site_name": "K63",
            "sequence_id": "NIQ(K)EST",
            "children": ""
          },
          {
            "site_name": "K48",
            "sequence_id": "FAG(K)QLE",
            "children": ""
          },
          // ... other lysine sites for protein 2 ...
        ]
      }
    },
    // ... other lysine sites for protein 1 ...
  ]
}
```

This structure enables representation of both simple and highly branched ubiquitin chains. All simulation and visualization in Ubiquitinformatics is based on this format.

Note: Protein 1 can also be a root protein other than ubiquitin, such as GFP or a His-tagged ubiquitin. This allows the system to represent fusion proteins or tagged constructs as the starting point of the chain.

---

We now introduce the principal functions underpinning the simulation framework, all implemented in `main.py`:
- iterate_through_ubiquitin
- ubiquitin_simulation

The output of `iterate_through_ubiquitin` is twofold:
1. An updated ubiquitin graph (dictionary/JSON) with chain numbers assigned by preorder traversal.
2. A context dictionary that encapsulates key properties of the ubiquitin graph, including:
   - `chain_number_list`: All chain numbers in the graph
   - `chain_length_list`: Lengths of each monomer
   - `multimer_string_name`: String representation of the chain architecture
   - `max_chain_number`: Total number of ubiquitin units
   - `ABOC_lysines`, `SMAC_lysines`: Lysines with protecting groups
   - `free_lysines`: Lysines available for further conjugation
   - `conjugated_lysines`: Lysines already conjugated to another ubiquitin

This context file summarizes the key properties of the ubiquitin graph and is used for analysis of the graph.

The `ubiquitin_simulation` function models the addition, removal, or modification of ubiquitin units and protecting groups within the graph. It applies these changes recursively, updating both the graph structure and its context to reflect the outcome of each simulated reaction.

**Example: E2 K63 Reaction**
Given the following ubiG JSON:
```json
{
  "protein": "1ubq",
  "chain_number": 1,
  "FASTA_sequence": "...",
  "chain_length": 76,
  "branching_sites": [
    {
      "site_name": "K63",
      "sequence_id": "NIQ(K)EST",
      "children": ""
    },
    {
      "site_name": "K48",
      "sequence_id": "FAG(K)QLE",
      "children": {
        "protein": "1ubq",
        "chain_number": 2,
        "FASTA_sequence": "...",
        "chain_length": 76,
        "branching_sites": [
          {
            "site_name": "K63",
            "sequence_id": "NIQ(K)EST",
            "children": ""
          },
          {
            "site_name": "K48",
            "sequence_id": "FAG(K)QLE",
            "children": ""
          }
          // ... other lysine sites for protein 2 ...
        ]
      }
    }
    // ... other lysine sites for protein 1 ...
  ]
}
```
Calling `ubiquitin_simulation(ubiG, new_ubiquitin, 'K63')` will add a new ubiquitin unit at the K63 site of chain 1, updating the `children` field of the K63 branching site.

Note: For E2 reactions, the code uses the `sequence_id` (e.g., "NIQ(K)EST") rather than just the site name (e.g., "K63"). This ensures biological relevance when identifying conjugation sites and accounts for differences in ubiquitin sequence length.

**Example: SMAC Deprotection**
If the K63 site has a protecting group:
```json
{
  "site_name": "K63",
  "sequence_id": "NIQ(K)EST",
  "children": "SMAC"
}
```
Calling `ubiquitin_simulation(ubiG, '', 'SMAC_deprot')` will remove the SMAC group, setting `children` to "" and making the site available for conjugation.

---

## 🧬 Simulation Workflow

The simulation proceeds by sequentially applying key functions within `create_reaction_histories`. At each step, `simulate_E2_steps` and `simulate_deprot_steps` are used to expand a growing set of ubiG graphs. For every addition, the type of reaction is noted and the appropriate enzyme is assigned using `assign_correct_E2_enzyme`. Unproductive pathways—those that do not result in a reaction or that exceed the allowed number of ubiquitin units—are systematically excluded from further analysis.

## 🔄 Breadth-First Search (BFS) and Pathway Filtering

The core simulation forms a breadth-first search (BFS) over possible reaction pathways, systematically exploring all valid chain topologies. This BFS approach ensures that all productive pathways are considered, while undesired or unproductive pathways are filtered out. The BFS logic is implemented in `run_file.py`, which also utilizes utility functions from `data_cleaning.py` for preprocessing and validation. 

### ABOC Filtering and SMAC Limiting
Once the initial database of chain architectures is created, additional filtering is applied in `run_file.py` to restrict entries based on the presence of ABOC protecting groups and the number of SMAC groups. Instead of excluding pathways during simulation, the filtering process analyzes the final products in the reaction database and assesses how many ABOC or SMAC groups each product contains.

## 📚 Comprehensive Linkage Database Generation

In addition to the main simulation, `run_file.py` applies functions from `all_linkages.py` to perform a second, comprehensive BFS. This process generates a database of every possible dimer to pentamer configuration, cataloging all valid linkage architectures. This database supports downstream analysis and visualization of ubiquitin chain diversity.

## 💾 Database Output

Both the filtered database (with ABOC/SMAC constraints) and the comprehensive linkage database are saved as CSV files in `back_end/Data`. These CSVs provide a structured record of all simulated chain architectures and are used for further analysis, visualization, and experimental design. In addition, a database of JSONs representing individual chain architectures is also stored in this directory, supporting more detailed inspection and downstream computational workflows.

For further details, see the relevant sections in `run_file.py`, `all_linkages.py`, and `data_cleaning.py`.

## 🏷️ Graph Isomorphism and Nomenclature Systems

Lastly, `all_linkages.py` handles graph isomorphism by directly inserting the JSON chain architecture into NetworkX, a Python package supporting advanced graph theory applications. This enables robust comparison and classification of chain topologies.

Additionally, `nomenclature.py` translates the logic from the JSON file into five defined nomenclature systems:
1. Preorder numbering where the nodes are labeled (e.g., A = 1, B = 2)
2. Preorder numbering where the edges are labeled (e.g., A = K63, B = K48)
3. Graph-based nomenclature with preorder numbering (Bode/Majima/Kummelstedt)
4. Chemistry-style all node labeled nomenclature (hover for info) (Bode)
5. Graph-based nomenclature without preorder (Strieter/Shestoperova/Ivanov)

These nomenclature systems provide multiple perspectives for describing and analyzing ubiquitin chain architectures, supporting both chemical and graph-theoretical interpretations.




