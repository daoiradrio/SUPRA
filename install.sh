#!/bin/bash

supra_path=$(pwd)

conda env create -f supra.yml
sudo ln -s $HOME/miniforge3/envs/supra/bin/python3 /usr/bin/supra-python3
ln -s $supra_path/"supra-conformer.py" ~/bin/supra-conformer
ln -s $supra_path/"supra-refine.py" ~/bin/supra-refine
ln -s $supra_path/"supra-doubles.py" ~/bin/supra-doubles
echo "export PYTHONPATH=$supra_path" >> ~/.bashrc
echo "export TURBODIR=$HOME/bin/TURBOMOLE" >> ~/.bashrc
echo "source $TURBODIR/Config_turbo_env" >> ~/.bashrc
