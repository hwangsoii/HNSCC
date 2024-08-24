#!/bin/bash
#$ -S /bin/bash

TEMP=/data/project/TRIUMPH/temp/sigprofiler_240513
if [ ! -d $TEMP ]
then
    mkdir -p $TEMP
fi

qsub -e $TEMP -o $TEMP -pe smp 10 -N runsig /data/project/TRIUMPH/script/TRIUMPH/run_sigprofiler.sh