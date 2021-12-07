
#source your conda environment first
#Then source this file as "source bash_profile.sh"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/oracle_instant_client/instantclient_21_1


#TEMP SOLUTION. DELETE LATER.
export MEDS_DIR=/home/dhayaa/Desktop/DECADE/TEMP_MEDS/

#Need this to point python to desmeds master file
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/desmeds_master:$PYTHONPATH
