# Backend Component Structure

## Overview
The ubiquitinformatics backend is a Python-based scientific computing platform designed for simulating, analyzing, and synthesizing complex ubiquitin chain structures. Built around a modular architecture, it provides comprehensive functionality for protein biochemistry research, synthesis planning, and nomenclature standardization.

## Core Architecture

### 1. **main.py** - Central Orchestrator
**Purpose**: Core engine for ubiquitin chain processing and validation
**Key Functions**:
- `iterate_through_ubiquitin()`: Primary processing function that traverses ubiquitin structures using pre-order tree traversal
- `ubiquitin_simulation()`: Simulates chemical reactions (K48/K63 conjugation, deprotection)
- `ubiquitin_building()`: Constructs new ubiquitin chains by adding monomers to specific lysine sites
- `validate_protein_keys()`: Ensures structural integrity of ubiquitin dictionaries
- `validate_branching_sites()`: Validates all 8 required lysine sites (M1, K6, K11, K27, K29, K33, K48, K63)

**Technical Implementation**:
- Recursive processing with context tracking for chain numbering and nomenclature generation
- Comprehensive error handling with detailed assertion messages
- Support for protecting groups (SMAC/ABOC) and free lysine identification
- JSON/dictionary interconversion with validation

**Dependencies**: utils, logging_utils, building_blocks
**Data Flow**: Receives ubiquitin structures → processes through tree traversal → outputs numbered chains with contexts

### 2. **fast_api.py** - API Gateway
**Purpose**: FastAPI-based REST interface providing endpoints for frontend integration
**Key Endpoints**:
- `/api/submit-selection`: Processes tetramers/pentamers selection for synthesis planning
- `/api/submit-ubxy`: Handles UbX_Y nomenclature queries and conversions
- `/api/submit-json-output`: Processes direct JSON structure submissions
- `/api/analyze-subgraphs`: Performs containment analysis between multimer sets
- `/api/analyze-subgraphs-stream`: Streaming version with real-time progress updates
- `/api/reaction-path-statistics`: Generates statistical analysis of reaction pathways
- `/api/submit_nomenclature_request`: Comprehensive nomenclature conversion service

**Technical Implementation**:
- CORS middleware for local development (localhost:5173)
- Base64 encoding for file transfers (Excel, Python scripts, PNGs, JSON)
- Streaming responses for long-running analyses
- Dynamic path resolution for cross-platform compatibility
- Comprehensive error handling with JSON responses

**Dependencies**: plotting, main, nomenclature, all_linkages
**Data Flow**: HTTP requests → data processing → base64-encoded responses with multiple file formats

### 3. **simulation.py** - Reaction Engine
**Purpose**: Biochemical reaction simulation for ubiquitin synthesis pathways
**Key Functions**:
- `simulate_E2_steps()`: Simulates E2 ligase reactions at K48/K63 sites
- `simulate_deprot_steps()`: Models SMAC deprotection and buffer wash reactions
- `assign_correct_E2_enzyme()`: Determines appropriate E2 enzyme based on topology
- `create_reaction_histories()`: Builds complete synthesis pathways for multimer generation

**Technical Implementation**:
- Reaction validation preventing self-conjugation (e.g., K48_SMAC + K48_reaction)
- Enzyme assignment logic: gp78/Ube2g2 vs Ube2K for K48 reactions based on context
- Chain number tracking to ensure single ubiquitin addition per reaction
- Support for elongation vs branching determination

**Dependencies**: main, utils, building_blocks
**Data Flow**: Acceptor + donor molecules → reaction simulation → product structures with enzyme assignments

### 4. **nomenclature.py** - Scientific Notation System
**Purpose**: Comprehensive nomenclature system supporting 7 different notation formats
**Key Systems**:
- **UbID Format**: A48B-B63C (preorder node labeling with lysine numbers)
- **Preorder ABC**: 1A2-2B3 (numbered nodes with lysine letters)
- **Graph-based**: Ub1,63(Ub2,48(Ub3)) (nested structural representation)
- **Chemical All-Node**: Complex position-based system for mass spectrometry
- **Conjugated Lysines**: [[1,'K63',2],...] (explicit edge representation)

**Key Functions**:
- `parse_compact_edges()`: Converts compact notation to explicit edge lists
- `build_polyubiquitin_from_edges()`: Constructs structures from edge specifications
- `format_nomenclature_preorder_1A2()`: Generates preorder nomenclature
- `conjugated_lysines_to_chemical_all_node_nomenclature()`: Complex chemical notation
- `build_mass_spec_dictionary()`: Extracts FASTA sequences for mass spectrometry

**Technical Implementation**:
- Bidirectional conversion between all nomenclature formats
- Preorder traversal preservation across format conversions
- Lysine mapping dictionaries for systematic notation
- Mass spectrometry FASTA sequence generation

**Dependencies**: main, simulation, utils, building_blocks
**Data Flow**: Structure representations → nomenclature conversion → standardized scientific notation

### 5. **plotting.py** - Visualization and Analysis Engine
**Purpose**: 96-well plate visualization, synthesis planning, and statistical analysis
**Key Functions**:
- `plot_96wells()`: Creates 96-well plate heatmaps with custom colormaps
- `create_xlsx_bytes()`: Generates Excel files with reagent calculations
- `create_opentrons_file_bytes()`: Produces Python scripts for laboratory automation
- `reaction_path_statistics()`: Analyzes linkage patterns and reaction pathways
- `build_reaction_dictionaries_for_UI()`: Prepares synthesis data for frontend

**Technical Implementation**:
- Matplotlib-based visualization with scientific color schemes
- Base64 encoding for web transfer
- Statistical analysis of multimer databases
- Laboratory automation script generation
- Interactive plotting with customizable parameters

**Dependencies**: matplotlib, pandas, numpy, main
**Data Flow**: Synthesis selections → plate layouts → visualization files + automation scripts

### 6. **all_linkages.py** - Extended Linkage Support
**Purpose**: Support for all 7 lysine linkage types beyond K48/K63
**Key Functions**:
- `process_branch_all()`: Handles all lysine types (M1, K6, K11, K27, K29, K33, K48, K63)
- `iterate_through_ubiquitin_all()`: Extended processing for comprehensive linkage support
- `load_multimer_contexts()`: Loads multimer data for analysis
- `analyze_subgraph_containment()`: Performs subgraph analysis between multimer sets

**Technical Implementation**:
- Extended branch processing logic for all 8 sites
- Subgraph containment algorithms for structural analysis
- Progress callback support for streaming operations
- Statistical analysis of multimer relationships

**Dependencies**: main, utils, building_blocks
**Data Flow**: Complex structures → all-linkage processing → comprehensive analysis results

### 7. **building_blocks.py** - Molecular Definitions
**Purpose**: Defines fundamental ubiquitin building blocks and molecular structures
**Key Components**:
- `ubi_ubq_1`: Standard ubiquitin monomer structure
- `histag_ubi_ubq_1`: His-tagged ubiquitin for purification
- `ubi_ubq_1_K48_SMAC`: K48-protected ubiquitin building block
- `ubi_ubq_1_K63_SMAC`: K63-protected ubiquitin building block

**Technical Implementation**:
- JSON-formatted molecular structures with complete FASTA sequences
- Branching site definitions for all 8 lysine positions
- Protecting group specifications (SMAC/ABOC)
- Chain length and sequence validation data

**Dependencies**: None (foundational module)
**Data Flow**: Provides molecular building blocks → used throughout system for construction

### 8. **data_cleaning.py** - Data Pipeline
**Purpose**: Data preprocessing and filtering for multimer databases
**Key Functions**:
- Database filtering for synthesis-compatible structures
- Data validation and consistency checking
- CSV processing for reaction histories
- Statistical data preparation

**Technical Implementation**:
- Pandas-based data manipulation
- Filtering algorithms for synthesis planning
- Data consistency validation
- Export preparation for analysis tools

**Dependencies**: pandas, main
**Data Flow**: Raw databases → filtered datasets → analysis-ready data

## Utility Modules

### **utils/utils.py** - Core Utilities
- JSON/dictionary conversion functions
- String matching and validation utilities
- Error handling helpers
- Data type conversion tools

### **utils/logging_utils.py** - Logging Framework
- Structured logging for debugging
- Protein processing trace logging
- Branch processing documentation
- Error tracking and reporting

## Data Architecture

### **Input Data Formats**:
- JSON ubiquitin structures with complete molecular definitions
- CSV databases (reaction_history, context_history, ubiquitin_history, donor_history)
- Compact nomenclature strings (A63B-B48C format)
- UbX_Y identifiers (Ub4_1, Ub5_847 format)

### **Output Data Formats**:
- Base64-encoded visualizations (PNG heatmaps)
- Excel files with synthesis calculations
- Python automation scripts for Opentrons
- JSON reaction sequences for frontend
- CSV statistical analysis results
- Mass spectrometry FASTA files

## API Integration Patterns

### **Request Processing Flow**:
1. FastAPI receives HTTP request with selection data
2. Data validation and format conversion
3. Database loading and filtering
4. Processing through appropriate simulation/analysis modules
5. Visualization generation and file creation
6. Base64 encoding for web transfer
7. JSON response with embedded files

### **Streaming Operations**:
- Real-time progress updates for long-running analyses
- Server-sent events for subgraph analysis
- Timing estimation and completion forecasting
- Error handling with graceful degradation

## Scientific Applications

### **Synthesis Planning**:
- Automated reaction sequence generation
- Enzyme selection optimization
- Protecting group strategy planning
- 96-well plate layout optimization

### **Structural Analysis**:
- Subgraph containment analysis
- Multimer relationship mapping
- Statistical linkage pattern analysis
- Nomenclature standardization

### **Laboratory Integration**:
- Opentrons automation script generation
- Reagent calculation spreadsheets
- Mass spectrometry sample preparation
- Quality control validation

## Performance Considerations

### **Optimization Strategies**:
- Recursive processing with context preservation
- Streaming for memory-efficient large dataset analysis
- Base64 encoding for efficient file transfer
- Caching of computed structures and contexts

### **Scalability Features**:
- Modular architecture enabling independent scaling
- Progress tracking for user experience
- Error isolation preventing cascade failures
- Resource monitoring and estimation

## Extensibility Guidelines

### **Adding New Linkage Types**:
1. Update `building_blocks.py` with new molecular definitions
2. Extend validation in `main.py` for new sites
3. Add nomenclature support in `nomenclature.py`
4. Update visualization in `plotting.py`

### **New Analysis Methods**:
1. Create new module following existing patterns
2. Add FastAPI endpoint in `fast_api.py`
3. Implement streaming if analysis is long-running
4. Add appropriate data export formats

### **Frontend Integration**:
1. Define clear API contracts
2. Implement comprehensive error handling
3. Use Base64 encoding for file transfers
4. Provide progress updates for long operations

This backend architecture provides a robust foundation for ubiquitin research, combining rigorous scientific accuracy with modern software engineering practices to enable complex biochemical analysis and synthesis planning.
