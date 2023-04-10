'''
Routines for moving files around, deleting remaining files
and saving relevant quantities.
'''

import shutil
import yaml
import os

def initialize_files(tilename, bands, output_desdata, config):
    
    #test limit:
    
#     for name, value in os.environ.items():
#         print("{0}: {1}".format(name, value))
        
#     for i in range(10000):
#         os.system("cp -v /scratch/midway2/dhayaa/PREP_OUTPUTS/v55_TestrunBeforeFull/outputs_DES0849+0252_seed0_gplus/simple_des_y3_sims/y3v02/band_info_files/DES0849+0252_r_info.yaml %s/tmp%d.txt" % (os.environ['TMPDIR'], i))
#         if i % 1000 == 0: os.system("du -sh %s" % (os.environ['TMPDIR'] ))
            
#     for i in range(200):
#         os.system("cp -v /scratch/midway2/dhayaa/TMP_DIR/meds-DES0821+0626-r/sources-r/DEC_Taiga/multiepoch/delve/r5918/DES0821+0626/p02/seg/DES0821+0626_r5918p02_r_segmap.fits %s/tmp%d.txt" % (os.environ['TMPDIR'], i))
#         os.system("du -sh %s" % (os.environ['TMPDIR'] ))
        
    move_prepdir_files(tilename, output_desdata)
    move_tmpdir_files(tilename, output_desdata)

    for b in bands:
        rewrite_band_infoing(tilename, b, output_desdata)
    
#Helper functions to run the above cleanup/re-organization code
def move_tmpdir_files(tile, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile}
    
    old_path = os.environ['SHEAR_TMPDIR'] + "/*%(tile)s*" % args
    new_path = os.environ['TMPDIR']
    
    print(old_path, new_path)
    os.system("cp -r %s %s" % (old_path, new_path))
    
    print("---------------------------------")
    os.system("du -sh %s" % (new_path))
    print("---------------------------------")
    
    os.system("tree %s" % (new_path))

#     shutil.copy(old_path, new_path)
        
    return True

#Helper functions to run the above cleanup/re-organization code
def move_prepdir_files(tile, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'basename': os.path.basename(output_desdata),
            'tile' : tile}
    
    old_path = os.environ['PREP_DIR'] + "/%(name)s/%(basename)s" % args
    new_path = os.environ['TMPDIR'] + "/%(name)s/%(basename)s" % args

    os.makedirs(os.environ['TMPDIR'] + "/%(name)s" % args)
    
    print(old_path, new_path)
    os.system("mv %s %s" % (old_path, new_path))
    
    
    print("---------------------------------")
    os.system("du -sh %s" % (new_path))
    print("---------------------------------")
    
    os.system("tree %s" % (new_path))
        
    return True


def rewrite_band_infoing(tile, band, output_desdata):
    
    def rewrite(X, a, b):
        
        for k in X.keys():
            
            if isinstance(X[k], str):
                X[k] = X[k].replace(a, b)
         
            if isinstance(X[k], dict):
                X[k] = rewrite(X, a, b)
                
            if k == 'src_info':
                for i in range(len(X[k])):
                    X[k][i] = rewrite(X[k][i], a, b)
                
        return X
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile,
            'basename': os.path.basename(output_desdata),
            'band' : band}
    
    #This is always going to be path in our runs so just sort of hardcode this assumption
    band_path = os.environ['TMPDIR'] + '/%(name)s/%(basename)s/simple_des_y3_sims/y3v02/band_info_files/%(tile)s_%(band)s_info.yaml'%args
    
    with open(band_path, 'r') as fp:
        band_info = yaml.load(fp, Loader=yaml.Loader)
    
    old_path = os.environ['SHEAR_TMPDIR']
    new_path = os.environ['TMPDIR']
    
    band_info_new = rewrite(band_info, old_path, new_path)
        
    with open(band_path, 'w') as fp:
        yaml.dump(band_info, fp)
         
    return True