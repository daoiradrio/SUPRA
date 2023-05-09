# SUPRA

***Installation***
1. Clone SUPRA repository
2. Install miniforge: https://github.com/conda-forge/miniforge
3. ```conda env create -f supra.yml```, if environment isn't automatically activated: ```conda activate supra```
5. ```ln -s /path/to/python3 /usr/bin/supra-python3``` (```/path/to/python3``` from ```which python3``` in activated supra environment)
6. ```ln -s /path/to/SUPRA/supra-conformer.py ~/bin/supra-conformer```
7. ```ln -s /path/to/SUPRA/supra-refine.py ~/bin/supra-refine```
8. Add ```export PYTHONPATH=/path/to/SUPRA``` to ```~/.bashrc``` 
9. Get ```TURBOMOLE```, add ```export TURBODIR=$HOME/bin/TURBOMOLE``` and ```source $TURBODIR/Config_turbo_env``` to ```~/.bashrc```
10. Put```xtb``` binary in ```~/bin```
