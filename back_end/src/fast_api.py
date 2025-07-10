from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import base64
import io

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
    fig1, ax1 = plotting.plot_96wells(cdata=enzymes_donors_96, figure=1, figure_name='Plate map: Enzyme + Donor Mixes', colorbar_type='PuRd')
    fig2, ax2 = plotting.plot_96wells(cdata=deprots_96, figure=2, figure_name='Plate map: Deprotections', colorbar_type='BuGn')
    fig3, ax3 = plotting.plot_96wells(cdata=acceptors_96, figure=3, figure_name='Plate map: Acceptors', colorbar_type='Blues')

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

    # Return as JSON with base64-encoded PNGs, Excel, and Opentrons Python file
    return JSONResponse(content={
        "received_labels": indexes,
        "page": page,
        "status": "ok",
        "figures": {
            "enzymes_donors_96": fig1_b64,
            "deprots_96": fig2_b64,
            "acceptors_96": fig3_b64,
            "reagent_calculations.xlsx": excel_b64,
            "opentrons.py": opentrons_b64
        }
    })


