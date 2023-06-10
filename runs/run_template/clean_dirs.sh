
#A convenience script that cleans up all directories
#that we load data into. Also removes all job.sh
#and log files from the current directory

rm -r $TMPDIR/*
rm -r $MEDS_DIR*
rm -r $PREP_DIR/*


rm //home/dhayaa/Desktop/DECADE/mcal_sim_test/runs/v55_TestrunBeforeFull/job_*
rm //home/dhayaa/Desktop/DECADE/mcal_sim_test/runs/v55_TestrunBeforeFull/*.log