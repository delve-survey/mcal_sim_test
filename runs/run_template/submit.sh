#!/bin/bash

for f in `ls job_*.sh`
do
    echo $f
    sbatch $f
done