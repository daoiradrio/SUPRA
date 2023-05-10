#!/bin/bash

conda env create -f supra.yml
conda activate supra

python_path=$(which python3)
supra_path=$(pwd)

sudo ln -s $python_path /usr/bin/supra-python3
ln -s $supra_path/"supra-conformer.py" ~/bin/supra-conformer
ln -s $supra_path/"supra-refine.py" ~/bin/supra-refine
ln -s $supra_path/"supra-doubles.py" ~/bin/supra-doubles
echo "export PYTHONPATH=$supra_path" >> ~/.bashrc
