# Ubiquitinformatics

A Python and React application for ubiquitin chain simulation and visualization.

**⚠️ Important: This application currently only works on Mac computers.**
## What You Need Before Starting

Before you can download and run this project, you need to have these programs installed on your Mac:

### Required Software
1. **Git** - For downloading the project from GitHub
   - Check if you have it: Open Terminal and type `git --version`
   - If you don't have it, download from: [https://git-scm.com/download/mac](https://git-scm.com/download/mac)

2. **Python 3.11.3** - The programming language this project uses
   - We recommend using `pyenv` to manage Python versions
   - Install pyenv: [https://github.com/pyenv/pyenv#installation](https://github.com/pyenv/pyenv#installation)

3. **Node.js** - For the web interface
   - Download from: [https://nodejs.org/](https://nodejs.org/) (choose the LTS version)

### Alternative: Conda (Optional)
Instead of the above Python setup, you can use Conda:
- Download Anaconda: [https://www.anaconda.com/download](https://www.anaconda.com/download)
- Or Miniconda: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

## How to Download and Run the Project

### Step 1: Download the Project
1. Open Terminal (Applications → Utilities → Terminal)
2. Navigate to where you want to put the project:
   ```bash
   cd Desktop  # This puts it on your Desktop, or choose another location
   ```
3. Download the project:
   ```bash
   git clone https://github.com/erickummelstedt/ubiquitinformatics.git
   ```
4. Go into the project folder:
   ```bash
   cd ubiquitinformatics
   ```

### Step 2: Set Up the Python Environment

**Option A: Using Python and pip (Recommended)**
```bash
pyenv local 3.11.3
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements_pip.txt
```

**Option B: Using Conda**
```bash
conda create -n myenv python
conda activate myenv
conda config --append channels conda-forge
conda install --file requirements_conda.txt
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
- If it doesn't, go to: [http://localhost:3000](http://localhost:3000)