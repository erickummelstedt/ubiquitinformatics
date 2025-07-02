#!/bin/sh
# Kill any process using port 5173-5180 (Vite default range)
for port in 5173 5174 5175 5176 5177 5178 5179 5180; do
  PID=$(lsof -ti :$port)
  if [ ! -z "$PID" ]; then
    kill -9 $PID 2>/dev/null
  fi
done

# Set up Python virtual environment and install requirements
cd back_end
if [ ! -d ".venv" ]; then
  python3 -m venv .venv
fi
. .venv/bin/activate
pip install --upgrade pip
pip install -r ../requirements_pip.txt
cd ..

cd front_end
npm run dev:all -- --port 5173 &
# Wait a few seconds for the server to start, then open the browser
sleep 3
open http://localhost:5173/
