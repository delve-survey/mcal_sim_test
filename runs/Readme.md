# Procedure to run

This is a (maybe incomplete) guide on how to use the files in this folder to quickly start image sim runs.

1. Change the `config.yaml.temp` file (the galsim input) according to the complexity you need in sims
2. Change the header lines in `job.sh.temp` to choose the Midway partition, account, mail-user and etc settings
3. Do `vim render_jobs.py` (or via other editors), and change the specific tiles you want to simulate.
4. Then do `python render_jobs.py`. This does four things: 

    (i) runs the `prep` version of the pipeline for every tile of interest, 
    
    (ii) creates a galsim config file with both the `g1 = 0.02` and `g1 = -0.02` version
    
    (iii) creates a `clean_dirs.sh` file that can be used to cleanup directories as needed.
    
    (iv) creates a `job_*.sh` file for each tile. One for g+ and g- version. Each bash file starts a job that continues the full pipeline to the end.
    
5. Now you should have a list of jobs in your directory. Run `source submit.sh` and this will submit all jobs.
6. Note that at the very end of each job, the code will move the metacal output to a different directory (make sure to set yours), and delete all other outputs from the pipeline (eg. galsim galaxies, meds files, etc.)
