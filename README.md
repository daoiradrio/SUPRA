# SUPRA

## Manual installation
1. Clone SUPRA repository
2. Install miniforge: https://github.com/conda-forge/miniforge
3. ```conda env create -f supra.yml```, if environment isn't automatically activated: ```conda activate supra```
5. ```ln -s /path/to/python3 /usr/bin/supra-python3``` (```/path/to/python3``` from ```which python3``` in activated supra environment)
6. ```ln -s /path/to/SUPRA/supra-conformer.py ~/bin/supra-conformer```
7. ```ln -s /path/to/SUPRA/supra-doubles.py ~/bin/supra-doubles```
8. ```ln -s /path/to/SUPRA/supra-refine.py ~/bin/supra-refine```
9. Add ```export PYTHONPATH=/path/to/SUPRA``` to ```~/.bashrc``` 
10. Put```xtb``` binary in ```~/bin```
11. ```source ~/.bashrc```
12. ```conda activate supra```

## (Partly) automatic installation
1. Clone SUPRA repository
2. Install miniforge: https://github.com/conda-forge/miniforge
3. ```./install.sh```
4. Put ```xtb``` binary in ```~/bin```
5. ```source ~/.bashrc```
6. ```conda activate supra```

## Usage
SUPRA offers interfaces with different functionality. These are explained in the following.
### 1. ```supra-conformer```
The ```supra-conformer``` module expects the path of a .xyz input file containing an optimized molecule structure, e.g. ```supra-conformer -path SUPRA/inputfiles/Alanin.xyz```. The generation of conformer structures starts after a mode (ignoring bonds to terminal bonds or not) and angle increment has been chosen. All structures are optimized with the Universal Force Field and if necessary converged with GFN2-FF. Possible duplicate structures are removed and the final conformer structures are placed in a directory ```SUPRA_Ouput``` that is placed in the same directory where ```supra-conformer``` has been called.
### 2. ```supra-refine```
The ```supra-refine``` module expects the path of a directory containing conformer structures, e.g. ```supra-refine -path SUPRA_Output/```. These can be generated with SUPRA, but don't have to be. It reoptimizes all structures using xTB and removes possible duplicate structures. The final reoptimized conformer structures are placed in the input directory.
### 3. ```supra-doubles```
The ```supra-doubles``` module expects the path of a directory containing conformer structures, e.g. ```supra-doubles -path SUPRA_Output/```. Again, these can be generated with SUPRA, but don't have to be. The module is used to remove duplicate structures without any optimization. In contrast to the other modules this module offers more specific options for deciding whether two conformer structures are equal or not: The optional ```-rmsd``` keyword can be used to adjust the RMSD threshold under which two structures are considered to be equal, e.g. ```supra-refine -path SUPRA_Output/ -rmsd -0.08```. The default value is 0.1. The optional ```-ignore``` keyword determines whether terminal groups are considered in the RMSD-based matching or not by typing ```supra-doubles -path SUPRA_Output/ -ignore all```. If only methyl alike groups (i.e. methyl groups with halides substituted for hydrogen) should be ignored the command is ```supra-doubles -path SUPRA_Output/ -ignore methyl```. Of course, both commands can be combined: ```supra-refine -path SUPRA_Output/ -rmsd -0.08 -ignore methyl```.
