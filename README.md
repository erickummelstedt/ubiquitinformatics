# Ubiquitinformatics

A Python and React application for ubiquitin chain simulation and visualization.

### 🐳 Option 1: Docker Setup (Recommended)
**No Python or Node.js installation needed!**
1. **Git** - For downloading the project
   - Check: `git --version` in Terminal
   - Install: [https://git-scm.com/download](https://git-scm.com/download)

2. **Docker Desktop** - Handles everything else automatically
   - Download: [https://www.docker.com/get-started](https://www.docker.com/get-started)
   - Install and start Docker Desktop
   - Ensure Docker is in your PATH: Run `docker --version` in Terminal to verify

```bash
# Step 1: Download the project
git clone https://github.com/erickummelstedt/ubiquitinformatics.git
cd ubiquitinformatics

# Step 2: Run everything with Docker (takes 2-3 minutes)
# Mac/Linux:
./start.sh

# Windows:
start.bat
```

Access at: http://localhost:5173

### 🛠️ Option 2: Manual Setup (If you prefer local installation)
**- requires installing Python and Node.js:**
```bash
# Download the project
git clone https://github.com/erickummelstedt/ubiquitinformatics.git
cd ubiquitinformatics

# Run everything locally
# Mac/Linux:
./manual-start.sh
```

This will set up Python environment, install dependencies, run simulation, and start the web interface.

---

---

## 📋 What You Need First


### Option 2: Manual Setup (Advanced Users)
1. **Git** + **Python 3.11+** + **Node.js** (see manual setup section below)

---

## ⚠️ Manual Setup (If You Need It)
pip install -r requirements.txt
```

**Option B: Using Conda**
```bash
conda create -n myenv python
conda activate myenv
conda config --append channels conda-forge
conda install --file requirements.txt
```

### Step 3: Run the Application
You need to run these two commands (in this order):

1. Start the Python backend:
   ```bash
   python3 back_end/src/run_file.py
   ```

2. In a **new Terminal window/tab**, start the web interface:
   ```bash
   ./start-dev.sh
   ```

### Step 4: Use the Application
- After running both commands, your web browser should automatically open
- If it doesn't, please [report an issue](https://github.com/erickummelstedt/ubiquitinformatics/issues) with your operating system details