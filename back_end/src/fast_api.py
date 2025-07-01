from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware

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
    labels = data.get("labels", [])
    page = data.get("page", "")
    # TODO: Integrate your Python logic here
    # For now, just echo back the received data
    return {"received_labels": labels, "page": page, "status": "ok"}

