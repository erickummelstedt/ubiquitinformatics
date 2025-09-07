#!/bin/bash

echo "🛑 Stopping Ubiquitinformatics..."

# Stop with docker compose
if command -v "docker compose" &> /dev/null; then
    docker compose down
elif command -v docker-compose &> /dev/null; then
    docker-compose down
else
    echo "❌ Docker Compose not found."
    exit 1
fi

echo ""
echo "✅ Application stopped!"
echo ""
echo "Your simulation data is preserved for next time."
echo "To start again, run: ./start.sh"