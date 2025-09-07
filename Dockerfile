# Simple all-in-one Dockerfile for Ubiquitinformatics
FROM python:3.11-slim

# Install Node.js and other dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire project
COPY . .

# Install frontend dependencies
WORKDIR /app/front_end
RUN npm install

# Go back to main directory
WORKDIR /app

# Create a startup script that runs the simulation once, then starts the web interface
RUN echo '#!/bin/bash\n\
echo "🔬 Running simulation (this may take a few minutes)..."\n\
cd /app/back_end/src\n\
python run_file.py\n\
echo "✅ Simulation completed!"\n\
echo "🌐 Starting web interface..."\n\
cd /app/front_end\n\
npx concurrently "npm run dev -- --host 0.0.0.0 --port 5173" "cd ../back_end && uvicorn src.fast_api:app --reload --host 0.0.0.0 --port 8000"' > /app/start_app.sh && chmod +x /app/start_app.sh

# Expose the ports
EXPOSE 5173
EXPOSE 8000

# Run the startup script
CMD ["/app/start_app.sh"]