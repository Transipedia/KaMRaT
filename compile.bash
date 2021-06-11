#!/bin/bash

conda_env_path=$1

if [ ! $conda_env_path ]; then
	sed 's|-I/path_to_conda_mlpack_env/include||g' src/Makefile.Template > src/Makefile
else
	sed "s|/path_to_conda_mlpack_env/|$conda_env_path|g" src/Makefile.Template > src/Makefile
fi

mkdir bin
make -C src

make -C related-tools/prepare_kmer_table/dekupl-joinCounts # dekupl-joinCounts in related-tools/
mv related-tools/prepare_kmer_table/dekupl-joinCounts/joinCounts bin/
mv related-tools/prepare_kmer_table/splitCV.bash bin/
mv related-tools/prepare_kmer_table/revCompFastq.pl bin/
chmod +x bin/*
