#!/bin/bash

conda_env_path=$1

if [ ! $conda_env_path ]; then
	sed 's|-I/path_to_conda_mlpack_env/include||g' src/Makefile.Template > src/Makefile
else
	sed "s|/path_to_conda_mlpack_env/|$conda_env_path|g" src/Makefile.Template > src/Makefile
fi

mkdir bin
cd src
make
