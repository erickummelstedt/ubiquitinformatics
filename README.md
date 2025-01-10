# ubiquitinformatics

# License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

```Bash
pyenv local 3.11.3
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements_pip.txt
```

```Conda 
conda create -n myenv python
conda activate myenv
conda config --append channels conda-forge
conda install --file requirements_conda.txt


``` NOTES
openbabel is easily install in conda
opentrons is easily installed in pip/venv