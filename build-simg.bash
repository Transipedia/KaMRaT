#!/bin/bash

mkdir bin/
singularity build --fakeroot bin/kamrat.sif kamrat.simg.def
singularity run bin/kamrat.sif
