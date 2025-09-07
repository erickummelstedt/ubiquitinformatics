@echo off

echo 🛑 Stopping Ubiquitinformatics...
echo.

REM Stop and remove containers
docker compose down >nul 2>&1
if errorlevel 1 (
    docker-compose down >nul 2>&1
    if errorlevel 1 (
        echo ❌ Failed to stop containers. Make sure Docker is running.
        pause
        exit /b 1
    )
)

echo ✅ Application stopped!
echo.
echo Your simulation data is preserved for next time.
echo To start again, run: start.bat
echo.
echo Press any key to continue...
pause >nul
