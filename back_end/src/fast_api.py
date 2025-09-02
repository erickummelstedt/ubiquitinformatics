from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
import base64
import io
import json
import asyncio
import logging

# Configure logging to write to a file instead of the terminal
log_file_path = "application.log"
file_handler = logging.FileHandler(log_file_path)
file_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)

# Clear existing handlers and add the file handler
logging.getLogger().handlers.clear()
logging.getLogger().addHandler(file_handler)
logger = logging.getLogger(__name__)

app = FastAPI()

# Allow CORS for local frontend development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],  # Use your frontend URL for development
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/api/submit-selection")
async def submit_selection(request: Request):
    data = await request.json()
    selected_ids = data.get("labels", [])
    page = data.get("page", "")

    import sys
    from pathlib import Path
    import pandas as pd

    # Dynamically get the backend path relative to this file's location
    current_path = Path(__file__).resolve()
    project_root = current_path.parents[2]  # Go up to project root
    sys.path.insert(0, str(project_root))
    local_path = project_root / 'back_end'
    sys.path.insert(0, str(local_path))

    # Import project modules
    import src.utils.utils
    import src.utils.logging_utils
    import src.main as main
    import src.plotting as plotting

    # Determine multimer size from page
    if page == 'tetramers':
        multimer_size = 4
    elif page == 'pentamers':
        multimer_size = 5
    else:
        return {"error": "Invalid page value", "status": "error"}

    # Download CSV files
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

    # Download the data dictionary
    data_dict = download_data_dict(multimer_size)
    combined_database = data_dict['combined_database']
    context_history = data_dict['context_history']
    donor_history = data_dict['donor_history']
    reaction_history = data_dict['reaction_history']
    ubiquitin_history = data_dict['ubiquitin_history']

    # Convert selected ids to a list of indexes
    indexes = []
    for id in selected_ids:
        try:
            new_index = int(combined_database[(combined_database['multimer_id'] == id) & (combined_database['used_in_synthesis'] == 1)]['index'].unique()[0])
            indexes.append(new_index)
        except Exception as e:
            return {"error": f"ID {id} not found or invalid: {str(e)}", "status": "error"}

    # Create plate dataframes
    output_dict = plotting.inner_create_plate_dfs(data_dict, indexes, multimer_size)

    # Extract the plate dataframes
    enzymes_donors_96 = output_dict['enzymes_donors_96']
    deprots_96 = output_dict['deprots_96']
    acceptors_96 = output_dict['dimer_acceptors_96']

    # Generate plots
    fig1, ax1 = plotting.plot_96wells(cdata=enzymes_donors_96, figure=1, figure_name='Plate map: Enzyme + Donor Mixes', colorbar_type='enzymes_donors')
    fig2, ax2 = plotting.plot_deprotection_cycles(cdata=deprots_96, figure=2, figure_name='Deprotection Cycles', colorbar_type='deprots')
    fig3, ax3 = plotting.plot_96wells(cdata=acceptors_96, figure=3, figure_name='Plate map: Acceptors', colorbar_type='acceptors')

    # Helper to convert figure to base64
    def fig_to_base64(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        img_bytes = buf.read()
        buf.close()
        return base64.b64encode(img_bytes).decode('utf-8')

    fig1_b64 = fig_to_base64(fig1)
    fig2_b64 = fig_to_base64(fig2)
    fig3_b64 = fig_to_base64(fig3)

    # Generate Excel file as base64
    excel_bytes = plotting.create_xlsx_bytes(output_dict)
    excel_bytes.seek(0)
    excel_b64 = base64.b64encode(excel_bytes.read()).decode('utf-8')

    # Generate python opentrons file as base64
    opentrons_bytes = plotting.create_opentrons_file_bytes(output_dict)
    opentrons_bytes.seek(0)
    opentrons_b64 = base64.b64encode(opentrons_bytes.read()).decode('utf-8')

    # Convert reaction_sequences_dicts to bytes
    reaction_sequences_dicts = plotting.build_reaction_dictionaries_for_UI(data_dict, indexes, multimer_size)
    reaction_sequences_bytes = io.BytesIO()
    reaction_sequences_bytes.write(json.dumps(reaction_sequences_dicts).encode('utf-8'))
    reaction_sequences_bytes.seek(0)
    reaction_sequences_b64 = base64.b64encode(reaction_sequences_bytes.read()).decode('utf-8')

    # Return as JSON with base64-encoded PNGs, Excel, Opentrons Python file, and reaction sequences
    return JSONResponse(content={
        "received_labels": indexes,
        "page": page,
        "status": "ok",
        "figures": {
            "enzymes_donors_96": fig1_b64,
            "deprots_96": fig2_b64,
            "acceptors_96": fig3_b64,
            "reagent_calculations.xlsx": excel_b64,
            "opentrons.py": opentrons_b64,
            "reaction_sequences.json": reaction_sequences_b64
        }
    })

# New endpoint to handle UbX_Y submission
@app.post("/api/submit-ubxy")
async def submit_ubxy(request: Request):

    import sys
    from pathlib import Path
    import pandas as pd
    import copy

    # Dynamically get the backend path relative to this file's location
    current_path = Path(__file__).resolve()
    project_root = current_path.parents[2]  # Go up to project root
    sys.path.insert(0, str(project_root))
    local_path = project_root / 'back_end'
    sys.path.insert(0, str(local_path))

    # Import project modules
    import src.utils.utils
    import src.utils.logging_utils
    import src.main as main
    import src.plotting as plotting
    import src.nomenclature as nomenclature
    
    try:
        data = await request.json()
        ubxy_value = data.get("ubxy", None)
        if not ubxy_value:
            return JSONResponse(content={"status": "error", "message": "No UbX_Y value provided"}, status_code=400)

        # Validate that ubxy_value begins with 'U' or 'A'
        if not (ubxy_value.startswith('U') or ubxy_value.startswith('1')):
            return JSONResponse(content={"status": "error", "message": "UbX_Y value must begin with 'U' or '1'"}, status_code=400)

        logger.info(f"Received UbX_Y value: {ubxy_value}")

        # If it starts with 'A', fix the nomenclature
        if ubxy_value.startswith('1'):
            # Original UbX_Y processing for values starting with 'U'
            multimer_size = int(nomenclature.multimer_length_from_nomenclature(ubxy_value))  # Rough estimate based on length
        else:
            # Original UbX_Y processing for values starting with 'U'
            multimer_size = int(ubxy_value.replace("Ub", "").split('_')[0])
        
        # Function to load data
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

        # Load the data
        data_dict = download_data_dict(multimer_size)
        combined_database = data_dict['combined_database']
        context_history = data_dict['context_history']
        donor_history = data_dict['donor_history']
        reaction_history = data_dict['reaction_history']
        ubiquitin_history = data_dict['ubiquitin_history']

        # If it starts with 'A', convert from nomenclature format to UbX_Y format 
        if ubxy_value.startswith('1'):
            # Convert 'A1B2C3' style to 'U4_1' style
            try:
                nomenclature_value = ubxy_value  # e.g., 'A1B2C3'
                
                parsed_edges = nomenclature.parse_compact_edges(ubxy_value)
                print("parsed_edges: ", parsed_edges)   
                output_ubiG_json = nomenclature.build_polyubiquitin_from_edges(parsed_edges)
                print("output_ubiG_json: ", output_ubiG_json )
                # Load multimers JSON file and convert to dictionary
                file_path1 = project_root / 'front_end' / 'src' / 'data' / f'multimer_id_to_json{multimer_size}.json'
                with open(file_path1, 'r') as f:
                    multimers_dict = json.load(f)

                # Find the key in multimers_dict that corresponds to output_ubiG_json
                print("output_ubiG_json: ", output_ubiG_json )
                print("ubxy_value found:", ubxy_value)

                ubxy_value = None
                for key, value in multimers_dict.items():
                    if value == str(output_ubiG_json):
                        ubxy_value = key
                        break
                if ubxy_value is None:
                    return JSONResponse(content={"status": "error", "message": "Generated structure not found in multimers database"}, status_code=404)
                    
                print("ubxy_value found:", ubxy_value)

                logger.info(f"Converted {nomenclature_value} to {ubxy_value}")  # Debug print
            except Exception as e:
                return JSONResponse(content={"status": "error", "message": f"Error in nomenclature conversion: {str(e)}"}, status_code=400)

        # Get unique indices for the specified multimer_id and ensure all values are integers
        indexes = ubiquitin_history.loc[ubiquitin_history["multimer_id"] == ubxy_value, "index"].dropna().unique()
        indexes = [int(i) for i in indexes]

        if len(indexes) > 0:
            final_multimer = ubiquitin_history[ubiquitin_history["index"] == indexes[0]]["final_multimer"].iloc[0]
            ubxy_value = ubiquitin_history[ubiquitin_history["index"] == indexes[0]]["multimer_id"].iloc[0]
        else:
            logger.info(f"No indexes found for multimer_id: {ubxy_value}")

        # Pull the final multimer json and context
        output_json, output_context = main.iterate_through_ubiquitin(final_multimer)
        
        # Pull formatted edges
        edges = output_context['conjugated_lysines']
        formatted_edges = ', '.join([f"{src} -> {site} -> {dst}" for src, site, dst in edges])
        
        # New version - pre-order with A, B, C
        nomenclature_preorder_ABC = nomenclature.format_nomenclature_preorder_ABC(edges)

        # New version - Jeffs pre-order with K63, K48, K33, K29, K27, K11, K6 and A, B ,C as nodes
        nomenclature_preorder_jeff = nomenclature.format_nomenclature_preorder_jeff(edges)
        # ========== 

        # Correct nomenclature assignments
        strieter_nomenclature_wo_preorder = output_context['nomenclature_wo_preorder']
        kummelstedt_nomenclature_w_preorder = output_context['nomenclature_w_preorder']
        
        # Old version - Jeffs without pre-order - NEED TO FIX AGAIN
        # jeff_k48_k63_nomenclature = nomenclature.conjugated_lysines_to_nomenclature(edges)
        jeff_K48_K63_nomenclature = nomenclature.conjugated_lysines_to_jeff_K48_K63_nomenclature(edges)
        # jeff_full_nomenclature_ABC = ....
        jeff_all_lysines_nomenclature = nomenclature.conjugated_lysines_to_jeff_all_lysines_nomenclature(edges)
        # jeff_multiple_symbols = ....
        jeff_multiple_symbols = nomenclature.conjugated_lysines_to_jeffs_multiple_symbols(edges)
        # jeff_full_nomenclature_numbered = ....
        jeff_numerical_nomenclature = nomenclature.tree_nomenclature_to_numerical_system(jeff_all_lysines_nomenclature)
        # jeff_multiple_symbols_eric_numbering
        jeff_multiple_symbols_eric_numbering = nomenclature.conjugated_lysines_to_jeffs_multiple_symbols_eric_numbering(edges)

        # Convert reaction_sequences_dicts to bytes
        reaction_sequences_dicts = plotting.build_reaction_dictionaries_for_UI(data_dict, indexes, multimer_size)
        reaction_sequences_bytes = io.BytesIO()
        reaction_sequences_bytes.write(json.dumps(reaction_sequences_dicts).encode('utf-8'))
        reaction_sequences_bytes.seek(0)
        reaction_sequences_b64 = base64.b64encode(reaction_sequences_bytes.read()).decode('utf-8')

        # Return the entered UbX_Y value
        return JSONResponse(content={
            "status": "ok", "ubxy": ubxy_value, 
            "reaction_sequences_b64": reaction_sequences_b64,
            "formatted_edges": formatted_edges,
            "nomenclature_preorder_ABC": nomenclature_preorder_ABC,
            "nomenclature_preorder_jeff": nomenclature_preorder_jeff,
            "strieter_nomenclature_wo_preorder": strieter_nomenclature_wo_preorder,
            "kummelstedt_nomenclature_w_preorder": kummelstedt_nomenclature_w_preorder,
            "jeff_K48_K63_nomenclature": jeff_K48_K63_nomenclature,
            "jeff_all_lysines_nomenclature": jeff_all_lysines_nomenclature,
            "jeff_multiple_symbols": jeff_multiple_symbols,
            "jeff_numerical_nomenclature": jeff_numerical_nomenclature,
            "jeff_multiple_symbols_eric_numbering": jeff_multiple_symbols_eric_numbering
            })
    except Exception as e:
        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=500)

# New endpoint to handle jsonOutput submission
@app.post("/api/submit-json-output")
async def submit_json_output(request: Request):
    import sys
    from pathlib import Path
    import pandas as pd

    # Dynamically get the backend path relative to this file's location
    current_path = Path(__file__).resolve()
    project_root = current_path.parents[2]  # Go up to project root
    sys.path.insert(0, str(project_root))
    local_path = project_root / 'back_end'
    sys.path.insert(0, str(local_path))

    # Import project modules
    import src.utils.utils
    import src.utils.logging_utils
    import src.main as main
    import src.plotting as plotting
    import src.nomenclature as nomenclature

    try:
        data = await request.json()
        json_output = data.get("jsonOutput", None)
        if not json_output:
            return JSONResponse(content={"status": "error", "message": "No jsonOutput provided"}, status_code=400)

        # Assuming jsonOutput is a list of dictionaries, determine multimer_size
        multimer_size = len(json_output) + 1

        # Function to load data
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

        # Load the data
        data_dict = download_data_dict(multimer_size)
        combined_database = data_dict['combined_database']
        context_history = data_dict['context_history']
        donor_history = data_dict['donor_history']
        reaction_history = data_dict['reaction_history']
        ubiquitin_history = data_dict['ubiquitin_history']

        # Get indexes for the final multimer from jsonOutput
        indexes = plotting.get_indexes_for_final_multimer(json_output, ubiquitin_history)

        if len(indexes) > 0:
            final_multimer = ubiquitin_history[ubiquitin_history["index"] == indexes[0]]["final_multimer"].iloc[0]
            ubxy_value = ubiquitin_history[ubiquitin_history["index"] == indexes[0]]["multimer_id"].iloc[0]
        else:
            logger.info(f"No indexes found for multimer_id: {ubxy_value}")

        # Pull the final multimer json and context
        output_json, output_context = main.iterate_through_ubiquitin(final_multimer)
        
        # Pull formatted edges
        edges = output_context['conjugated_lysines']
        formatted_edges = ', '.join([f"{src} -> {site} -> {dst}" for src, site, dst in edges])
        
        # New version - pre-order with A, B, C
        nomenclature_preorder_ABC = nomenclature.format_nomenclature_preorder_ABC(edges)

        # New version - Jeffs pre-order with K63, K48, K33, K29, K27, K11, K6 and A, B ,C as nodes
        nomenclature_preorder_jeff = nomenclature.format_nomenclature_preorder_jeff(edges)
        
        # Correct nomenclature assignments
        strieter_nomenclature_wo_preorder = output_context['nomenclature_wo_preorder']
        kummelstedt_nomenclature_w_preorder = output_context['nomenclature_w_preorder']
        
        # Old version - Jeffs without pre-order - NEED TO FIX AGAIN
        # jeff_k48_k63_nomenclature = nomenclature.conjugated_lysines_to_nomenclature(edges)
        jeff_K48_K63_nomenclature = nomenclature.conjugated_lysines_to_jeff_K48_K63_nomenclature(edges)
        # jeff_full_nomenclature_ABC = ....
        jeff_all_lysines_nomenclature = nomenclature.conjugated_lysines_to_jeff_all_lysines_nomenclature(edges)
        # jeff_multiple_symbols = ....
        jeff_multiple_symbols = nomenclature.conjugated_lysines_to_jeffs_multiple_symbols(edges)
        # jeff_full_nomenclature_numbered = ....
        jeff_numerical_nomenclature = nomenclature.tree_nomenclature_to_numerical_system(jeff_all_lysines_nomenclature)
        # jeff_multiple_symbols_eric_numbering
        jeff_multiple_symbols_eric_numbering = nomenclature.conjugated_lysines_to_jeffs_multiple_symbols_eric_numbering(edges)

        # Convert reaction_sequences_dicts to bytes
        reaction_sequences_dicts = plotting.build_reaction_dictionaries_for_UI(data_dict, indexes, multimer_size)
        reaction_sequences_bytes = io.BytesIO()
        reaction_sequences_bytes.write(json.dumps(reaction_sequences_dicts).encode('utf-8'))
        reaction_sequences_bytes.seek(0)
        reaction_sequences_b64 = base64.b64encode(reaction_sequences_bytes.read()).decode('utf-8')

        # Return the entered UbX_Y value
        return JSONResponse(content={
            "status": "ok", "ubxy": ubxy_value, 
            "reaction_sequences_b64": reaction_sequences_b64,
            "formatted_edges": formatted_edges,
            "nomenclature_preorder_ABC": nomenclature_preorder_ABC,
            "nomenclature_preorder_jeff": nomenclature_preorder_jeff,
            "strieter_nomenclature_wo_preorder": strieter_nomenclature_wo_preorder,
            "kummelstedt_nomenclature_w_preorder": kummelstedt_nomenclature_w_preorder,
            "jeff_K48_K63_nomenclature": jeff_K48_K63_nomenclature,
            "jeff_all_lysines_nomenclature": jeff_all_lysines_nomenclature,
            "jeff_multiple_symbols": jeff_multiple_symbols,
            "jeff_numerical_nomenclature": jeff_numerical_nomenclature,
            "jeff_multiple_symbols_eric_numbering": jeff_multiple_symbols_eric_numbering
            })
    except Exception as e:
        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=500)


# New endpoint for subgraph analysis with streaming
@app.post("/api/analyze-subgraphs-stream")
async def analyze_subgraphs_stream(request: Request):
    """
    Streaming subgraph containment analysis with real-time timing information
    """
    import sys
    from pathlib import Path
    
    # Dynamically get the backend path relative to this file's location
    current_path = Path(__file__).resolve()
    project_root = current_path.parents[2]  # Go up to project root
    sys.path.insert(0, str(project_root))
    local_path = project_root / 'back_end'
    sys.path.insert(0, str(local_path))

    # Import project modules
    import src.main as main
    import src.all_linkages as linkages
    
    try:
        import pandas as pd
        import base64
        import time
        
        data = await request.json()
        higher_level_size = data.get("higher_level_size", 5)  # Default to pentamers
        n_level_size = data.get("n_level_size", 4)  # Default to tetramers
        higher_level_lysine_ids = set(data.get("higher_level_lysine_ids", ["K48", "K63"]))  # Default to K48/K63
        n_level_lysine_ids = set(data.get("n_level_lysine_ids", ["K48", "K63"]))  # Default to K48/K63

        # Load the data
        higher_level_data = linkages.load_multimer_contexts(project_root, higher_level_size)
        n_level_data = linkages.load_multimer_contexts(project_root, n_level_size)

        # Filter by lysine types
        higher_level_dict = linkages.get_multimer_edges_by_lysines(higher_level_data, higher_level_lysine_ids)
        n_level_dict = linkages.get_multimer_edges_by_lysines(n_level_data, n_level_lysine_ids)

        async def stream_analysis():
            """Generator function for streaming analysis results"""
            results = None
            
            # Send initial status
            yield f"data: {json.dumps({'type': 'status', 'message': 'Starting analysis...'})}\n\n"
            
            # Run the analysis with progress callback
            analysis_start_time = time.time()
            
            # We need to run the analysis in a way that allows async streaming
            import threading
            import queue
            
            result_queue = queue.Queue()
            update_queue = queue.Queue()
            
            def run_analysis():
                def streaming_callback(progress_data):
                    # Put all progress updates into the queue for streaming
                    update_queue.put(progress_data)
                
                try:
                    results = linkages.analyze_subgraph_containment(higher_level_dict, n_level_dict, streaming_callback)
                    result_queue.put(("success", results))
                except Exception as e:
                    result_queue.put(("error", str(e)))
            
            # Start analysis in separate thread
            analysis_thread = threading.Thread(target=run_analysis)
            analysis_thread.start()
            
            # Stream updates while analysis runs
            while analysis_thread.is_alive():
                try:
                    # Check for any updates (timing or progress)
                    while not update_queue.empty():
                        progress_data = update_queue.get_nowait()
                        
                        if progress_data.get("type") == "timing_analysis":
                            # This is the key timing update after 10 iterations!
                            timing_event = {
                                "type": "timing_update",
                                "data": {
                                    "completed_iterations": progress_data.get("completed_iterations"),
                                    "elapsed_time": progress_data.get("elapsed_time"),
                                    "avg_time_per_iteration": progress_data.get("avg_time_per_iteration"),
                                    "estimated_total_time": progress_data.get("estimated_total_time"),
                                    "estimated_remaining_time": progress_data.get("estimated_remaining_time"),
                                    "estimated_total_seconds": progress_data.get("estimated_total_seconds"),
                                    "estimated_remaining_seconds": progress_data.get("estimated_remaining_seconds")
                                }
                            }
                            # Ensure proper JSON encoding
                            try:
                                json_str = json.dumps(timing_event, ensure_ascii=True)
                                yield f"data: {json_str}\n\n"
                            except Exception as json_error:
                                print(f"JSON encoding error: {json_error}")
                                yield f"data: {json.dumps({'type': 'error', 'message': 'JSON encoding error'})}\n\n"
                        
                        elif progress_data.get("type") == "progress":
                            # Regular progress updates
                            progress_event = {
                                "type": "progress_update",
                                "data": {
                                    "current": progress_data.get("current"),
                                    "total": progress_data.get("total"),
                                    "message": progress_data.get("message", "")
                                }
                            }
                            # Ensure proper JSON encoding
                            try:
                                json_str = json.dumps(progress_event, ensure_ascii=True)
                                yield f"data: {json_str}\n\n"
                            except Exception as json_error:
                                print(f"JSON encoding error: {json_error}")
                                yield f"data: {json.dumps({'type': 'error', 'message': 'JSON encoding error'})}\n\n"
                    
                    await asyncio.sleep(0.1)  # Small delay to prevent busy waiting
                except queue.Empty:
                    pass
            
            # Get final results
            try:
                result_type, result_data = result_queue.get_nowait()
                if result_type == "error":
                    yield f"data: {json.dumps({'type': 'error', 'message': result_data})}\n\n"
                    return
                else:
                    results = result_data
            except queue.Empty:
                yield f"data: {json.dumps({'type': 'error', 'message': 'Analysis completed but no results received'})}\n\n"
                return
            
            total_analysis_time = time.time() - analysis_start_time

            # Convert results to CSV
            df = pd.DataFrame.from_dict(results, orient='index').fillna(0).astype(int)
            csv_bytes = df.to_csv().encode('utf-8')
            csv_b64 = base64.b64encode(csv_bytes).decode('utf-8')

            # Send final results
            final_result = {
                "type": "final_results",
                "data": {
                    "status": "ok",
                    "results": results,
                    "csv_b64": csv_b64,
                    "analysis_metadata": {
                        "total_analysis_time": total_analysis_time,
                        "higher_level_count": len(higher_level_dict),
                        "n_level_count": len(n_level_dict),
                        "total_comparisons": len(higher_level_dict) * len(n_level_dict)
                    }
                }
            }

            # Ensure proper JSON encoding for final results
            try:
                json_str = json.dumps(final_result, ensure_ascii=True)
                yield f"data: {json_str}\n\n"
            except Exception as json_error:
                print(f"JSON encoding error for final results: {json_error}")
                yield f"data: {json.dumps({'type': 'error', 'message': 'JSON encoding error for final results'})}\n\n"

        return StreamingResponse(
            stream_analysis(),
            media_type="text/plain",
            headers={
                "Cache-Control": "no-cache",
                "Connection": "keep-alive",
                "Content-Type": "text/event-stream",
            }
        )
        
    except Exception as e:
        # For streaming errors, we need to return a regular JSON response
        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=500)

# New endpoint for subgraph analysis (simplified)
@app.post("/api/analyze-subgraphs")
async def analyze_subgraphs(request: Request):
    """
    Simple subgraph containment analysis with timing information
    """
    import sys
    from pathlib import Path
    
    # Dynamically get the backend path relative to this file's location
    current_path = Path(__file__).resolve()
    project_root = current_path.parents[2]  # Go up to project root
    sys.path.insert(0, str(project_root))
    local_path = project_root / 'back_end'
    sys.path.insert(0, str(local_path))

    # Import project modules
    import src.main as main
    import src.all_linkages as linkages

    
    try:
        import pandas as pd
        import base64
        import time
        
        data = await request.json()
        higher_level_size = data.get("higher_level_size", 5)  # Default to pentamers
        n_level_size = data.get("n_level_size", 4)  # Default to tetramers
        higher_level_lysine_ids = set(data.get("higher_level_lysine_ids", ["K48", "K63"]))  # Default to K48/K63
        n_level_lysine_ids = set(data.get("n_level_lysine_ids", ["K48", "K63"]))  # Default to K48/K63

        # Load the data
        higher_level_data = linkages.load_multimer_contexts(project_root, higher_level_size)
        n_level_data = linkages.load_multimer_contexts(project_root, n_level_size)

        # Filter by lysine types
        higher_level_dict = linkages.get_multimer_edges_by_lysines(higher_level_data, higher_level_lysine_ids)
        n_level_dict = linkages.get_multimer_edges_by_lysines(n_level_data, n_level_lysine_ids)

        # Store timing information for response
        timing_info = {}
        
        def progress_callback(progress_data):
            """Callback to capture timing information"""
            if progress_data.get("type") == "timing_analysis":
                nonlocal timing_info
                timing_info = {
                    "completed_iterations": progress_data.get("completed_iterations"),
                    "elapsed_time": progress_data.get("elapsed_time"),
                    "avg_time_per_iteration": progress_data.get("avg_time_per_iteration"),
                    "estimated_total_time": progress_data.get("estimated_total_time"),
                    "estimated_remaining_time": progress_data.get("estimated_remaining_time"),
                    "estimated_total_seconds": progress_data.get("estimated_total_seconds"),
                    "estimated_remaining_seconds": progress_data.get("estimated_remaining_seconds")
                }

        # Run the analysis with progress callback
        analysis_start_time = time.time()
        results = linkages.analyze_subgraph_containment(higher_level_dict, n_level_dict, progress_callback)
        total_analysis_time = time.time() - analysis_start_time

        # Convert results to CSV bytes
        df = pd.DataFrame.from_dict(results, orient='index').fillna(0).astype(int)
        csv_bytes = df.to_csv().encode('utf-8')
        csv_b64 = base64.b64encode(csv_bytes).decode('utf-8')

        response_content = {
            "status": "ok",
            "results": results,
            "csv_b64": csv_b64,
            "analysis_metadata": {
                "total_analysis_time": total_analysis_time,
                "higher_level_count": len(higher_level_dict),
                "n_level_count": len(n_level_dict),
                "total_comparisons": len(higher_level_dict) * len(n_level_dict)
            }
        }
        
        # Add timing information if available
        if timing_info:
            response_content["timing_analysis"] = timing_info

        return JSONResponse(content=response_content)
        
    except Exception as e:
        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=500)


@app.post("/api/reaction-path-statistics")
async def reaction_path_statistics_endpoint(request: Request):
    """
    FastAPI endpoint that uses reaction_path_statistics from plotting module.
    Analyzes reaction pathways and linkage patterns for multimers.
    """
    try:
        data = await request.json()
        multimer_size = data.get("multimer_size", 5)
        pathway_type = data.get("pathway_type", "all")

        import sys
        from pathlib import Path
        import pandas as pd

        # Dynamically get the backend path relative to this file's location
        current_path = Path(__file__).resolve()
        project_root = current_path.parents[2]  # Go up to project root
        sys.path.insert(0, str(project_root))
        local_path = project_root / 'back_end'
        sys.path.insert(0, str(local_path))

        # Import project modules (with error handling)
        try:
            import src.main as main
            import src.plotting as plotting
            import src.all_linkages as linkages
            import src.data_cleaning as data_cleaning
        except ImportError as e:
            return JSONResponse(content={"status": "error", "message": f"Import error: {str(e)}"}, status_code=500)

        if pathway_type == 'aboc':# Function to load filtered data
            
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

            # Load the filtered data
            data_dict = download_data_dict(multimer_size)
            context_history = data_dict['context_history']
            ubiquitin_history = data_dict['ubiquitin_history']
        
        elif pathway_type == 'all':
            # Function to load all data (unfiltered)
            def download_all_data_dict(multimer_size):
                input_dir = project_root / 'back_end' / 'data' / 'reaction_database' / f'multimer_size_{multimer_size}'
                context_history = pd.read_csv(input_dir / 'context_history.csv', index_col=0)
                donor_history = pd.read_csv(input_dir / 'donor_history.csv', index_col=0)
                reaction_history = pd.read_csv(input_dir / 'reaction_history.csv', index_col=0)
                ubiquitin_history = pd.read_csv(input_dir / 'ubiquitin_history.csv', index_col=0)
                return {
                    'context_history': context_history,
                    'donor_history': donor_history,
                    'reaction_history': reaction_history,
                    'ubiquitin_history': ubiquitin_history
                }

            # Load the all data
            data_dict = download_all_data_dict(multimer_size)
            context_history = data_dict['context_history']
            ubiquitin_history = data_dict['ubiquitin_history']

        # Load multimers JSON file
        file_path = project_root / 'front_end' / 'src' / 'data' / f'multimer_id_to_json{multimer_size}.json'
        with open(file_path, 'r') as f:
            multimers = json.load(f)

        # Call the reaction_path_statistics function
        json_with_reaction_information = plotting.reaction_path_statistics(ubiquitin_history, context_history, multimers, project_root, multimer_size)
        
        return JSONResponse(content={
            "status": "ok",
            "data": json_with_reaction_information,
            "multimer_size": multimer_size,
            "message": f"Successfully analyzed reaction path statistics for {len(json_with_reaction_information)} multimers of size {multimer_size} ({'Aboc-saturated' if pathway_type == 'aboc' else 'all'} pathways)"
        })
        
    except Exception as e:
        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=500)



# New endpoint to handle UbX_Y submission
@app.post("/api/submit_nomenclature_request")
async def submit_nomenclature_request(request: Request):

    import sys
    from pathlib import Path
    import pandas as pd
    import copy

    # Dynamically get the backend path relative to this file's location
    current_path = Path(__file__).resolve()
    project_root = current_path.parents[2]  # Go up to project root
    sys.path.insert(0, str(project_root))
    local_path = project_root / 'back_end'
    sys.path.insert(0, str(local_path))

    # Import project modules
    import src.utils.utils
    import src.utils.logging_utils
    import src.main as main
    import src.plotting as plotting
    import src.nomenclature as nomenclature
    
    try:
        data = await request.json()
        ubxy_value = data.get("ubxy", None)
        if not ubxy_value:
            return JSONResponse(content={"status": "error", "message": "No UbX_Y value provided"}, status_code=400)

        # Validate that ubxy_value begins with 'U' or 'A'
        if not (ubxy_value.startswith('U') or ubxy_value.startswith('1')):
            return JSONResponse(content={"status": "error", "message": "UbX_Y value must begin with 'U' or '1'"}, status_code=400)

        logger.info(f"Received UbX_Y value: {ubxy_value}")

        # Original UbX_Y processing for values starting with 'U'
        # Split UbX_Y value into X and Y components
        ubxy_parts = ubxy_value.replace("Ub", "").split('_')
        X = int(ubxy_parts[0])  # Multimer size (4 or 5)
        Y = int(ubxy_parts[1])  # Index number
        multimer_size = X
        
        # Function to load JSON data
        def download_jsons(multimer_size):
            input_dir = project_root / 'back_end' / 'data' / 'all_jsons'
            
            # Load multimers JSON files
            with open(input_dir / f'{multimer_size}_multimers_jsons.json', 'r') as f:
                multimer_jsons = json.load(f)
            
            # Load contexts JSON files  
            with open(input_dir / f'{multimer_size}_multimers_contexts.json', 'r') as f:
                multimer_contexts = json.load(f)
                
            return {
                'multimer_jsons': multimer_jsons,
                'multimer_contexts': multimer_contexts
            }

        # Load the data
        data_dict = download_jsons(multimer_size)
        multimer_jsons = data_dict['multimer_jsons']
        multimer_contexts = data_dict['multimer_contexts']

        # Find the multimer in the JSON data
        if str(Y) in multimer_jsons:
            final_multimer = multimer_jsons[str(Y)]
        else:
            return JSONResponse(content={"status": "error", "message": f"Multimer {ubxy_value} not found in database"}, status_code=404)
        
        # Pull the final multimer json and context
        output_json, output_context = main.iterate_through_ubiquitin(final_multimer)
        
        # Convert output_json dictionary to JSON string for saving
        output_json_string = json.dumps(output_json, indent=2)

        # Pull formatted edges
        edges = output_context['conjugated_lysines']
        formatted_edges = ', '.join([f"{src} -> {site} -> {dst}" for src, site, dst in edges])

        # New version - pre-order with A, B, C
        nomenclature_preorder_ABC = nomenclature.format_nomenclature_preorder_ABC(edges)

        # New version - Jeffs pre-order with K63, K48, K33, K29, K27, K11, K6 and A, B ,C as nodes
        nomenclature_preorder_jeff = nomenclature.format_nomenclature_preorder_jeff(edges)
        
        # Correct nomenclature assignments
        strieter_nomenclature_wo_preorder = output_context['nomenclature_wo_preorder']
        kummelstedt_nomenclature_w_preorder = output_context['nomenclature_w_preorder']

        # ======= Strieter-style FASTA generation =======
        
        # Add mass spec dictionary to output_json
        FASTA_sequences = nomenclature.build_mass_spec_dictionary(output_json)

        # Build the text file content
        txt_file_content = ">Ppt"
        
        # Add >Ppt section with FASTA sequences
        for key, fasta_sequence in FASTA_sequences.items():
            txt_file_content += f"\n{key};{fasta_sequence}\n"
        
        # Add >Isf section with nomenclature
        txt_file_content += f">Isf\n{ubxy_value};{strieter_nomenclature_wo_preorder}\n"
        
        # ======== End of Strieter-style FASTA generation =======
        
        # Old version - Jeffs without pre-order 
        # jeff_full_nomenclature_ABC = ....
        jeff_all_lysines_nomenclature = nomenclature.conjugated_lysines_to_jeff_all_lysines_nomenclature(edges)
        # jeff_full_nomenclature_numbered = ....
        jeff_numerical_nomenclature = nomenclature.tree_nomenclature_to_numerical_system(jeff_all_lysines_nomenclature)
        # jeff_multiple_symbols = ....
        jeff_multiple_symbols = nomenclature.conjugated_lysines_to_jeffs_multiple_symbols(edges)
        # jeff_multiple_symbols_eric_numbering
        jeff_multiple_symbols_eric_numbering = nomenclature.conjugated_lysines_to_jeffs_multiple_symbols_eric_numbering(edges)

        # Return the entered UbX_Y value
        return JSONResponse(content={
            "status": "ok", "ubxy": ubxy_value, 
            "output_json": output_json_string,
            "txt_file_content": txt_file_content,
            "formatted_edges": formatted_edges,
            "nomenclature_preorder_ABC": nomenclature_preorder_ABC,
            "nomenclature_preorder_jeff": nomenclature_preorder_jeff,
            "strieter_nomenclature_wo_preorder": strieter_nomenclature_wo_preorder,
            "kummelstedt_nomenclature_w_preorder": kummelstedt_nomenclature_w_preorder,
            "jeff_all_lysines_nomenclature": jeff_all_lysines_nomenclature,
            "jeff_numerical_nomenclature": jeff_numerical_nomenclature,
            "jeff_multiple_symbols": jeff_multiple_symbols,
            "jeff_multiple_symbols_eric_numbering": jeff_multiple_symbols_eric_numbering
            })
    except Exception as e:
        return JSONResponse(content={"status": "error", "message": str(e)}, status_code=500)
