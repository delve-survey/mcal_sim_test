import os
import jinja2
import yaml
import pandas as pd
import argparse
import numpy as np

my_parser = argparse.ArgumentParser()

my_parser.add_argument('--NoPrep', action='store_true', default = False)
my_parser.add_argument('--OverwriteJobs', action='store_true', default = False)

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

tilenames = ['DES0849+0252']

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



    #Create clean_dirs command for easy use
    with open('clean_dirs.sh.temp', 'r') as fp:
        tmp = jinja2.Template(fp.read())

    with open('clean_dirs.sh', 'w') as fp:
        fp.write(tmp.render(current_directory = os.getcwd()))



    #Now create all job.sh files for running sims
    with open('job.sh.temp', 'r') as fp:
        tmp = jinja2.Template(fp.read())

    for i, tilename in enumerate(tilenames):
        gal_seed = i
        if os.path.isfile('job_%s_minus.sh' % tilename) & os.path.isfile('job_%s_plus.sh' % tilename) & (args['OverwriteJobs'] == False):
            print("----------------------------")
            print("%s ALREADY PROCESSED"% tilename)
            print("SKIPPING PROCESSING")
            print("----------------------------")
            continue

        if args['NoPrep'] == False:
            os.system('python $RUN_DIR/run_sims.py prep --tilename="%s" --bands="riz" --output-desdata="$PREP_DIR/%s/outputs_%s_seed%d_gplus" --config-file="config_plus.yaml"'%(tilename, name, tilename, gal_seed))
            os.system('cp -r $PREP_DIR/%s/outputs_%s_seed%d_gplus $PREP_DIR/%s/outputs_%s_seed%d_gminus'%(name, tilename, gal_seed, name, tilename, gal_seed))

 
        with open('job_%s_plus.sh' % tilename, 'w') as fp:
            fp.write(tmp.render(tilename=tilename, model_name = name, plus_or_minus = "plus", seed_galsim=gal_seed, seed_mcal=42))
        with open('job_%s_minus.sh' % tilename, 'w') as fp:
            fp.write(tmp.render(tilename=tilename, model_name = name, plus_or_minus = "minus", seed_galsim=gal_seed, seed_mcal=42))
        os.system('chmod u+x job_%s_plus.sh' % tilename)
        os.system('chmod u+x job_%s_minus.sh' % tilename)
