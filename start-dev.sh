#!/bin/sh
# Kill any process using port 5173-5180 (Vite default range) and port 8000 (FastAPI)
echo "Checking for existing servers..."
for port in 5173 5174 5175 5176 5177 5178 5179 5180 8000; do
  PID=$(lsof -ti :$port)
  if [ ! -z "$PID" ]; then
    echo "Killing process $PID on port $port"
    kill -9 $PID 2>/dev/null
  fi
done

# Also kill any uvicorn processes
pkill -f "uvicorn.*fast_api" 2>/dev/null || true

# Check if virtual environment exists, create if it doesn't
if [ ! -d ".venv" ]; then
  echo "Virtual environment not found. Creating .venv..."
  python3 -m venv .venv
fi

echo "Activating virtual environment..."
source .venv/bin/activate

# Install Node.js and npm if missing
if ! command -v npm &> /dev/null; then
  echo "❌ 'npm' is not installed. Please install Node.js from https://nodejs.org/."
  exit 1
fi

# Install concurrently globally if missing
if ! command -v concurrently &> /dev/null; then
  echo "Installing 'concurrently' globally..."
  npm install -g concurrently
fi

# Install frontend dependencies
if [ -d "front_end" ]; then
  echo "Installing frontend dependencies..."
  cd front_end
  npm install
  cd ..
fi

# Check backend dependencies (venv should already be activated)
echo "Checking backend dependencies..."
if ! python -c "import uvicorn" 2>/dev/null; then
  echo "Installing backend dependencies..."
  pip install --upgrade pip
  pip install -r requirements_pip.txt
fi

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
