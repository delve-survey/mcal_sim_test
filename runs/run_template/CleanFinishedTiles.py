'''
This file searches the output directory
and finds all tiles that were run succesfully.
Then it renames both the jobscript and log files
as finished files so they are not run again.
'''

import os
import sys

from render_jobs import tilenames

print('--------------------------')
print('USING TILES:')
print('--------------------------')
for t in tilenames:
    print(t)
print('--------------------------')

N = len(tilenames)

name     = os.path.basename(os.path.dirname(__file__))
MCAL_DIR = os.environ['MCAL_DIR']
PATH     = MCAL_DIR + '/' + name
print('GETTING MCAL FILES FROM:')
print(PATH)

seed_list = [i for i in range(len(tilenames))]

#Loop over tiles (and also a combination of all tiles together)
#and check which ones exist
for tilename, seed in zip(tilenames, seed_list):
    
    for typ in ['plus', 'minus']:
        fname = PATH + '/metacal_%s_seed%d_g%s.fits'%(tilename, seed, typ)
        if os.path.isfile(fname): 
            
            jdir = os.path.dirname(__file__) + '/job_%s_%s.sh'%(tilename, typ) #job file directory
            if os.path.isfile(jdir): os.system(r'mv %s %s'%(jdir, jdir[:-3] + '.finished')) #Only run if file wasnt already processed
            
            ldir = os.path.dirname(__file__) + '/metacal_%s_seed%d_g%s.log'%(tilename, seed, typ) #log file directory
            if os.path.isfile(ldir): os.system(r'mv %s %s'%(ldir, ldir[:-4] + '.finished')) #Only run if file wasnt already processed
            
        else:
            print("NOT COMPUTED: %s (%s version). KEEPING SH FILES."%(tilename, typ))
