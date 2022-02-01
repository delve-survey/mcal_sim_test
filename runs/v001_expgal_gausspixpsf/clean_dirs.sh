
#A convenience script that cleans up all directories
#that we load data into. Also removes all job.sh
#and log files from the current directory

rm -r /home/dhayaa/scratch-midway2/TMP_DIR/*
rm -r /home/dhayaa/scratch-midway2/MEDS_DIR/*
rm -r /home/dhayaa/scratch-midway2/PREP_OUTPUTS/*


rm //home/dhayaa/Desktop/DECADE/mcal_sim_test/runs/v001_expgal_gausspixpsf/job_*
rm //home/dhayaa/Desktop/DECADE/mcal_sim_test/runs/v001_expgal_gausspixpsf/*.log