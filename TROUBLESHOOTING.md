# Simple Setup Troubleshooting

## If the application doesn't start:

### 1. Make sure Docker is running
- **Windows/Mac**: Look for the Docker whale icon in your system tray/menu bar
- **Linux**: Run `sudo systemctl start docker`

### 2. Check if port 8000 is already in use
```bash
# Find what's using port 8000
lsof -i :8000    # Mac/Linux
netstat -ano | findstr :8000    # Windows

# If something is using it, stop it or change the port
```

### 3. Still having problems?

**Reset everything and try again:**

```bash
# Stop any running containers
./stop.sh    # or stop.bat on Windows

# Remove old containers and start fresh
docker system prune -f
./start.sh   # or start.bat on Windows
```

**Check the logs:**
```bash
docker-compose logs
```

## Common Issues:

### "Docker is not installed"
- Download Docker Desktop from [docker.com](https://www.docker.com/get-started)
- Follow the installation instructions for your operating system

### "Docker is not running"
- **Windows/Mac**: Start Docker Desktop from your applications
- **Linux**: Run `sudo systemctl start docker`

### "Port already in use"
- Another application is using port 8000
- Close any other development servers or web applications
- Or change the port in `docker-compose.yml` from `8000:8000` to `8001:8000`

### Application loads but shows errors
- Wait a few more minutes - it can take time to download everything the first time
- Check the logs: `docker-compose logs`
- Try restarting: `./stop.sh` then `./start.sh`

### Need help?
Create an issue on GitHub with:
1. Your operating system (Windows/Mac/Linux)
2. The error message you're seeing
3. The output of `docker version`
