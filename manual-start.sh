#!/bin/bash

echo "🧬 Starting Ubiquitinformatics - Complete Setup"
echo "=============================================="
echo ""

# Kill any existing processes on ports first
echo "🔄 Checking for existing servers..."
for port in 5173 5174 5175 5176 5177 5178 5179 5180 8000; do
  PID=$(lsof -ti :$port)
  if [ ! -z "$PID" ]; then
    echo "Killing process $PID on port $port"
    kill -9 $PID 2>/dev/null
  fi
done

# Also kill any uvicorn processes
pkill -f "uvicorn.*fast_api" 2>/dev/null || true

# Step 1: Check if virtual environment exists
if [ ! -d ".venv" ]; then
  echo "📦 Setting up Python environment..."
  python3 -m venv .venv
fi

# Step 2: Activate virtual environment
echo "🔧 Activating Python environment..."
source .venv/bin/activate

# Step 3: Install Python dependencies if needed
if ! python -c "import uvicorn" 2>/dev/null; then
  echo "📥 Installing Python dependencies..."
  pip install --upgrade pip
  pip install -r requirements.txt
fi

# Step 4: Run the simulation ONCE
echo ""
echo "🔬 Running simulation (this may take a few minutes)..."
cd back_end/src
python run_file.py
cd ../..

echo ""
echo "✅ Simulation completed!"
echo ""

# Step 5: Install Node.js and npm if missing
if ! command -v npm &> /dev/null; then
  echo "❌ 'npm' is not installed. Please install Node.js from https://nodejs.org/."
  exit 1
fi

# Install concurrently globally if missing
if ! command -v concurrently &> /dev/null; then
  echo "Installing 'concurrently' globally..."
  npm install -g concurrently
fi

# Step 6: Install Node.js dependencies if needed
if [ -d "front_end" ]; then
  echo "📥 Installing frontend dependencies..."
  cd front_end
  npm install
  cd ..
fi

# Step 7: Start the web interface
echo "🌐 Starting web interface..."
echo ""
echo "This will open your browser automatically at http://localhost:5173"
echo ""
echo "Press Ctrl+C to stop the servers when you're done"
echo ""

# Set default port if not provided
PORT=${1:-5173}

echo "Starting development servers..."
cd front_end
npm run dev:all -- --port $PORT &
SERVER_PID=$!

# Wait a few seconds for the server to start, then open the browser
echo "Waiting for servers to start..."
sleep 5

# Check if servers are running
if lsof -i :$PORT > /dev/null && lsof -i :8000 > /dev/null; then
  echo "✅ Both servers started successfully!"
  echo "Frontend: http://localhost:$PORT"
  echo "Backend API: http://localhost:8000"
  open http://localhost:$PORT/
else
  echo "❌ Failed to start servers. Check for errors above."
fi

echo "Press Ctrl+C to stop all servers"
wait $SERVER_PID
