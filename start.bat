@echo off
setlocal enabledelayedexpansion

echo 🧬 Starting Ubiquitinformatics with Docker
echo ==========================================
echo.

REM Check if Docker is installed
docker --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Docker is not installed!
    echo.
    echo Please install Docker Desktop:
    echo 🔗 https://www.docker.com/get-started
    echo.
    echo After installing, restart this script.
    pause
    exit /b 1
)

REM Check if Docker is running, start it if not
docker info >nul 2>&1
if errorlevel 1 (
    echo 🔧 Docker is not running. Starting Docker Desktop...
    echo.
    start "" "Docker Desktop"
    echo ⏳ Waiting for Docker to start...
    
    REM Wait up to 60 seconds for Docker to start
    for /l %%i in (1,1,12) do (
        timeout /t 5 /nobreak >nul
        docker info >nul 2>&1
        if not errorlevel 1 (
            echo ✅ Docker is now running!
            goto docker_ready
        )
        if %%i==12 (
            echo ❌ Docker failed to start within 60 seconds.
            echo Please start Docker Desktop manually and try again.
            pause
            exit /b 1
        )
        set /a seconds=%%i*5
        echo    Still waiting... (!seconds! seconds elapsed^)
    )
)

:docker_ready
echo ✅ Docker is ready!
echo.
echo 🔧 Building application (this may take a few minutes the first time...)

REM Build and start with docker compose
docker compose version >nul 2>&1
if not errorlevel 1 (
    docker compose up --build -d
) else (
    docker-compose version >nul 2>&1
    if not errorlevel 1 (
        docker-compose up --build -d
    ) else (
        echo ❌ Docker Compose not found. Please make sure Docker Desktop is properly installed.
        pause
        exit /b 1
    )
)

echo.
echo ⏳ Waiting for simulation to complete and web interface to start...
echo    This may take 2-3 minutes the first time...

REM Wait for the services to be ready
timeout /t 5 /nobreak >nul
echo    🔬 Simulation running...

REM Check if the frontend is responding (wait up to 3 minutes)
for /l %%i in (1,1,36) do (
    curl -s http://localhost:5173 >nul 2>&1
    if not errorlevel 1 (
        echo    ✅ Web interface is ready!
        goto services_ready
    )
    timeout /t 5 /nobreak >nul
    set /a remainder=%%i %% 6
    if !remainder!==0 (
        set /a seconds=%%i*5
        echo    ⏳ Still waiting... (!seconds! seconds elapsed^)
    )
)

:services_ready
echo.
echo 🎉 SUCCESS! Ubiquitinformatics is now running!
echo.
echo 👉 Opening your web browser...
echo    Frontend: http://localhost:5173
echo    Backend API: http://localhost:8000
echo.

REM Open browser automatically
start http://localhost:5173

echo.
echo To stop the application later, run: stop.bat
echo To view logs, run: docker compose logs -f
echo.
echo Press any key to continue...
pause >nul