import os, io, subprocess as sp
import jinja2
import yaml
import pandas as pd
import argparse
import datetime as dt
import numpy as np
import glob

my_parser = argparse.ArgumentParser()

my_parser.add_argument('--Initialize', action='store_true', default = False)
my_parser.add_argument('--Maintain',   action='store_true', default = False)

my_parser.add_argument('--MaxConcurrentJobs', action='store', type = int, default = 10)
my_parser.add_argument('--MaxCutoffTime',     action='store', type = int, default = 3600*48) #Maxcutoff time in seconds

args = vars(my_parser.parse_args())

#Print args for debugging state
print('-------INPUT PARAMS----------')
for p in args.keys():
    print('%s : %s'%(p.upper(), args[p]))
print('-----------------------------')
print('-----------------------------')

TILENAME_SEED = 42
NUM_OF_TILES  = 40
tiles = pd.read_csv(os.environ['RUN_DIR'] + '/data/Tilelist_DR3_1_1.csv')
tilenames = list(np.random.default_rng(TILENAME_SEED).choice(tiles['TILENAME'].values, size = NUM_OF_TILES, replace = False))

tilenames = ['DES0821+0626', 'DES0849+0252']

# tilenames = ['DES0849+0252']

print(tilenames)
if __name__ == '__main__':
    
    #Automatically get folder name (assuming certain folder structure)
    name = os.path.basename(os.path.dirname(__file__))

    #Create output directory for metacal
    MCAL_DIR = os.environ['MCAL_DIR']
    os.makedirs(MCAL_DIR +'/'+ name, exist_ok=True)

    #Create the config file
    with open('config.yaml.temp', 'r') as fp:
        tmp = jinja2.Template(fp.read())


    for g_sign, g_multiplier in zip(['plus', 'minus'], [1, -1]): 
        with open('config_%s.yaml'%g_sign, 'w') as fp:
            fp.write(tmp.render(g1 = 0.02*g_multiplier))

        
    ########################################################################3
    
    def is_plus_finished(tilename):
        
        args  = {'name' : name, 'tile' : tilename}
        plus  = os.path.join(os.environ['MCAL_DIR'], "/%(name)s/metacal_%(tile)s_gplus.fits" % args)
        
        return os.path.isfile(plus)
    
    def is_minus_finished(tilename):
         
        args  = {'name' : name, 'tile' : tilename}
        minus = os.path.join(os.environ['MCAL_DIR'], "/%(name)s/metacal_%(tile)s_gminus.fits" % args)
        
        return os.path.isfile(minus)
    
    def is_finished(tilename):

        return is_plus_finished(tilename) and is_minus_finished(tilename)
    
    def is_plus_job(tilename):
        
        return os.path.isfile('job_%s_plus.sh' % tilename)
    
    def is_minus_job(tilename):
         
        return os.path.isfile('job_%s_minus.sh' % tilename)

        
    def create_job(tilename, gal_seed):
        
        #Now create all job.sh files for running sims
        with open('job.sh.temp', 'r') as fp:
            tmp = jinja2.Template(fp.read())
        
        if not is_plus_finished(tilename):
            with open('job_%s_plus.sh' % tilename, 'w') as fp:
                fp.write(tmp.render(tilename=tilename, model_name = name, plus_or_minus = "plus", seed_galsim=gal_seed, seed_mcal=42))
            os.system('chmod u+x job_%s_plus.sh' % tilename)
        
        if not is_minus_finished(tilename):
            with open('job_%s_minus.sh' % tilename, 'w') as fp:
                fp.write(tmp.render(tilename=tilename, model_name = name, plus_or_minus = "minus", seed_galsim=gal_seed, seed_mcal=42))

            os.system('chmod u+x job_%s_minus.sh' % tilename)
            
    def prep_job(tilename):
        
        os.system('python $RUN_DIR/run_sims.py prep --tilename="%s" --bands="riz" --output-desdata="$PREP_DIR/%s/outputs_%s_gplus" --config-file="config_plus.yaml"'%(tilename, name, tilename))
        os.system('cp -r $PREP_DIR/%s/outputs_%s_gplus $PREP_DIR/%s/outputs_%s_gminus'%(name, tilename, name, tilename))
        
        
    def current_job_count():
        
        x = sp.check_output("squeue --format='%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R' --sort=+i -u dhayaa", shell = True)
        j = pd.read_csv(io.StringIO(x.decode("utf-8")), delim_whitespace=True)
        
        count = 0
        tiles = []
        
        #metacal_DES0849+0252_seed0_gminus
        for i in range(len(j)):
            
            condition1 = 'metacal' in j['NAME'].values[i]
            condition2 = j['STATE'].values[i] == 'RUNNING'
            
            if condition1:
                count += 1
                tiles.append(j['NAME'].values[i][8:8+12])
                
        for t in tiles:
            print("CURRENTLY RUNNING: %s" % t)
                
        return count
            
    def tilename_from_jobfilename(job):
        
        #job_DES0849+0252_minus.sh
            
        return job[4:4+12]
    
                
    ################################################################
    
    #In initial step we just create every job that we need
    if args['Initialize']:
        
        for i, tilename in enumerate(tilenames):
        
            create_job(tilename, i)
    
    #In next step we keep track of all jobs and add to queue when we can
    elif args['Maintain']:
        
        job_list = sorted(glob.glob('job_*.sh'))
        start = dt.datetime.now()

        while (dt.datetime.now() - start).seconds < 3600*48: #48 hours max
            
            if len(job_list) == 0:
                print("ALL JOBS HAVE BEEN STARTED. NOTHING TO MAINTAIN")
                break
                
            if current_job_count() >= args['MaxConcurrentJobs']:
            
                time.sleep(120)

            else:
                
                j = job_list[0]
                t = tilename_from_jobfilename(j)
                
                prep_job(t)
                
                os.system('sbatch %s' % j)
                os.system('rm %s' % j)
                
                print("SUBMITTED JOB %s" % j)
                
                #This just removes the top item from the list
                #which we do as we just ran that job
                job_list = job_list[1:]
                
        else:
            
            print("GONE OVERTIME. SHUTTING DOWN SCRIPT")            
                
    