#!/bin/bash

DIR="~/miniconda3"
if [ -d "$DIR" ]; then
  ### Take action if $DIR exists ###
  echo "Virtual Machine is ready for use."
else
  cd ~/
  wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
  chmod +x Miniconda3-py39_4.12.0-Linux-x86_64.sh
  ./Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -u
  ~/miniconda3/bin/conda init
  cat >>~/.bashrc <<EOL
wget https://raw.githubusercontent.com/CosiMichele/workshop-materials/main/TE_symposium/setup_2.sh ~/setup_2.sh
chmod +x ~/setup_2.sh
~/setup_2.sh
rm ~/setup_2.sh
EOL
  exec bash
  exit 1
fi
