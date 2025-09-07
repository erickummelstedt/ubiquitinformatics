#!/bin/bash

echo "🧬 Starting Ubiquitinformatics with Docker"
echo "=========================================="
echo ""

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed!"
    echo ""
    echo "Please install Docker Desktop:"
    echo "🔗 https://www.docker.com/get-started"
    echo ""
    echo "After installing, restart this script."
    exit 1
fi

# Check if Docker is running, start it if not
if ! docker info &> /dev/null; then
    echo "🔧 Docker is not running. Starting Docker Desktop..."
    echo ""
    open -a Docker
    echo "⏳ Waiting for Docker to start..."
    
    # Wait up to 60 seconds for Docker to start
    for i in {1..12}; do
        sleep 5
        if docker info &> /dev/null; then
            echo "✅ Docker is now running!"
            break
        fi
        if [ $i -eq 12 ]; then
            echo "❌ Docker failed to start within 60 seconds."
            echo "Please start Docker Desktop manually and try again."
            exit 1
        fi
        echo "   Still waiting... (${i}0 seconds elapsed)"
    done
fi

echo "✅ Docker is ready!"
echo ""
echo "🔧 Building application (this may take a few minutes the first time)..."

# Build and start with docker compose
if command -v "docker compose" &> /dev/null; then
    docker compose up --build -d
elif command -v docker-compose &> /dev/null; then
    docker-compose up --build -d
else
    echo "❌ Docker Compose not found. Please make sure Docker Desktop is properly installed."
    exit 1
fi

echo ""
echo "⏳ Waiting for simulation to complete and web interface to start..."
echo "   This may take 2-3 minutes the first time..."

# Wait for the services to be ready
sleep 5
echo "   🔬 Simulation running..."

# Check if the frontend is responding (wait up to 3 minutes)
for i in {1..36}; do
    if curl -s http://localhost:5173 > /dev/null 2>&1; then
        echo "   ✅ Web interface is ready!"
        break
    fi
    sleep 5
    if [ $((i % 6)) -eq 0 ]; then
        echo "   ⏳ Still waiting... (${i}0 seconds elapsed)"
    fi
done

echo ""
echo "🎉 SUCCESS! Ubiquitinformatics is now running!"
echo ""
echo "👉 Opening your web browser..."
echo "   Frontend: http://localhost:5173"
echo "   Backend API: http://localhost:8000"
echo ""

# Try to open the browser automatically
if command -v open &> /dev/null; then
    open http://localhost:5173
elif command -v xdg-open &> /dev/null; then
    xdg-open http://localhost:5173
else
    echo "Please manually open your browser and go to: http://localhost:5173"
fi

echo ""
echo "To stop the application later, run: ./stop.sh"
echo "To view logs, run: docker compose logs -f"