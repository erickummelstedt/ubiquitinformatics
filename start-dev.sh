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

echo "Setting up Python environment..."
cd back_end
if [ ! -d ".venv" ]; then
  python3 -m venv .venv
fi
. .venv/bin/activate
pip install --upgrade pip
pip install -r ../requirements_pip.txt
cd ..

echo "Starting development servers..."
cd front_end
npm run dev:all -- --port 5173 &
SERVER_PID=$!

# Wait a few seconds for the server to start, then open the browser
echo "Waiting for servers to start..."
sleep 5

# Check if servers are running
if lsof -i :5173 > /dev/null && lsof -i :8000 > /dev/null; then
  echo "✅ Both servers started successfully!"
  echo "Frontend: http://localhost:5173"
  echo "Backend API: http://localhost:8000"
  open http://localhost:5173/
else
  echo "❌ Failed to start servers. Check for errors above."
fi

echo "Press Ctrl+C to stop all servers"
wait $SERVER_PID
