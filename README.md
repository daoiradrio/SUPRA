# SUPRA

***Manual installation***
1. Clone SUPRA repository
2. Install miniforge: https://github.com/conda-forge/miniforge
3. ```conda env create -f supra.yml```, if environment isn't automatically activated: ```conda activate supra```
5. ```ln -s /path/to/python3 /usr/bin/supra-python3``` (```/path/to/python3``` from ```which python3``` in activated supra environment)
6. ```ln -s /path/to/SUPRA/supra-conformer.py ~/bin/supra-conformer```
7. ```ln -s /path/to/SUPRA/supra-refine.py ~/bin/supra-refine```
8. Add ```export PYTHONPATH=/path/to/SUPRA``` to ```~/.bashrc``` 
9. Put ```TURBOMOLE``` in ```~/bin```, add ```export TURBODIR=$HOME/bin/TURBOMOLE``` and ```source $TURBODIR/Config_turbo_env``` to ```~/.bashrc```
10. Put```xtb``` binary in ```~/bin```
11. ```source ~/.bashrc```
12. ```conda activate supra```

***(Partly) automatic installation***
1. Clone SUPRA repository
2. Install miniforge: https://github.com/conda-forge/miniforge
3. ```./install.sh```
4. Put ```TURBOMOLE``` in ```~/bin```
5. Put ```xtb``` binary in ```~/bin```
6. ```source ~/.bashrc```
7. ```conda activate supra```

***Usage***
SUPRA offers interfaces with different functionality. These are explained in the following.
1. The ```supra-conformer``` module expects the path of a .xyz input file containing an optimized molecule structure, e.g. ```supra-conformer -path SUPRA/inputfiles/Alanin.xyz```. The generation of conformer structures starts after a mode (ignoring bonds to terminal bonds or not) and angle increment has been chosen. The generated conformer structures are placed in a directory ```SUPRA_Ouput``` that is placed in the same directory where ```supra-conformer``` has been called.
2. 
