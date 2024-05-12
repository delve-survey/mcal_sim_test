from astropy.io import fits
import numpy as np, healpy as hp
import os, sys, joblib
from tqdm import tqdm
import yaml, argparse
import fitsio

from SimButler import tilenames as tilenames_orig

#Some other packages for doing SOM stuff
sys.path.append('/home/dhayaa/Desktop/DECADE/CosmicShearPhotoZ/SOMPZ')
from SOM import Classifier


class SuppressPrint:
    def __enter__(self):
        self._original_stdout = sys.stdout  # Save a reference to the original standard output
        sys.stdout = open(os.devnull, 'w')  # Redirect standard output to a null device

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()  # Close the null device
        sys.stdout = self._original_stdout  # Restore the original standard output


def parallel_concatenator(col_caller, N, *args):
    
    def _func_(i): return i, col_caller(*tuple([a[i] for a in args]) )

    final_out = [0]*N
    with joblib.parallel_backend("loky"):
        jobs    = [joblib.delayed(_func_)(i) for i in range(N)]
        outputs = joblib.Parallel(n_jobs = 40, verbose=10,)(jobs)
        for o in outputs: final_out[o[0]] = o[1]
    
    return np.concatenate(final_out, axis = 0)

    
    
summary = np.mean

my_parser = argparse.ArgumentParser()

my_parser.add_argument('--TileCount',  action='store', type = int, default = 10)
my_parser.add_argument('--NewClassification',  action='store_true', default = False)

args = vars(my_parser.parse_args())

#Print args for debugging state
print('-------INPUT PARAMS----------')
for p in args.keys():
    print('%s : %s'%(p.upper(), args[p]))
print('-----------------------------')
print('-----------------------------')

tilenames = tilenames_orig[0:args['TileCount']]

#Flag maps
Foreground_map = hp.read_map('/project/chihway/dhayaa/DECADE/Foreground_Masks/GOLD_Ext0.2_Star5.fits', dtype = int)
Badcolor_map   = hp.read_map('/project2/kadrlica/chinyi/DELVE_DR3_1_bad_colour_mask.fits', dtype = int)
EBV            = hp.read_map('/project/chihway/dhayaa/DECADE/Imsim_Inputs/ebv_sfd98_fullres_nside_4096_ring_equatorial.fits')


# print('--------------------------')
# print('USING TILES:')
# print('--------------------------')
# for t in tilenames:
#     print(t)
# print('--------------------------')

N = len(tilenames)

name     = os.path.basename(os.path.dirname(__file__))
MCAL_DIR = os.environ['MCAL_DIR']
PATH     = MCAL_DIR + '/' + name
print('GETTING MCAL FILES FROM:')
print(PATH)

tiles_list = [[x] for x in tilenames]
tiles_list.append(tilenames)

seed_list = [[i] for i in range(len(tilenames))]
seed_list.append([i for i in range(len(tilenames))])

tiles_list = [tiles_list[-1]]
seed_list  = [seed_list[-1]]

with open(os.path.dirname(__file__) + '/config_plus.yaml', 'r') as fp:
    config = yaml.load(fp, Loader=yaml.Loader)
        
def _get_subset_array(res, i, string, b):
    '''
    Convenience function for getting reading out specific arrays
    and only for specific objects that are part of a jackknife "patch"
    '''
    
    if b is None:
        ind = np.where(res['batch_num_%s'%string[3:]] != i)[0] #Select only objects in given batch/subset
    else:
        ind = np.where( (res['batch_num_%s'%string[3:]] != i) & 
                        (res['bin_%s'%string[3:]] == b) )[0] #Select only objects in given batch/subset and tomobin
    
    return res[string][ind]


def _get_cat_path(tile, band, seed, mode = 'plus'):
    
    args = {'name' : name,
            'seed' : seed,
            'mode' : mode,
            'tile' : tile,
            'band' : band}
    
    
    return os.environ['MCAL_DIR'] + r"/%(name)s/SrcExtractor_%(tile)s_g%(mode)s_%(band)s-cat.fits" % args
    
#Loop over tiles (and also a combination of all tiles together)
#and get m and c for objects in each tile
for tiles, seed in zip(tiles_list, seed_list):
    
    _t = []
    _s = []
    
    for t, i in zip(tiles, seed):
        
        try:
            #Open (only ext = 0 since no data is loaded then), save to variable, and delete. Latter two steps to make sure memory is cleared
            x = fitsio.read(PATH + '/metacal_%s_gplus.fits'%(t), ext = 0); del x
            x = fitsio.read(PATH + '/metacal_%s_gminus.fits'%(t), ext = 0); del x
            _t.append(t)
            _s.append(i)
        except:
            print("SKIPPING", t)
        
        

    tiles = _t
    seed  = _s

    if len(tiles) == 0: continue

    N = len(tiles) 
    print("N_TILES:", N)    
    #if N == 1: continue #Just quick "hack" so we skip all single tile calculations
    
    print("--------------------------------------------------------")
    print("TILES:", tiles)
    print("--------------------------------------------------------")
    
    tinds  = range(len(tiles))
    gplus  = 0 #[fitsio.read(PATH + '/metacal_%s_gplus.fits'%(t)) for t in tiles]
    gminus = 0 #[fitsio.read(PATH + '/metacal_%s_gminus.fits'%(t)) for t, i in zip(tiles, seed)]
    
    if 'truedet' not in config['gal_kws']['truth_type']:
        print(_get_cat_path(tiles[0], 'r', seed[0], mode = 'plus'))        
        
        mask_SEflag_plus = ( (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'FLAGS'), 
                                                    N, tiles, tinds) <= 3) & 
                             (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'i', i, mode = 'plus'), columns = 'FLAGS'), 
                                                    N, tiles, tinds) <= 3) & 
                             (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'z', i, mode = 'plus'), columns = 'FLAGS'), 
                                                    N, tiles, tinds) <= 3)
                           )
        
        mask_snr_plus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'FLUX_AUTO'), 
                                              N, tiles, tinds) /
                         parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'FLUXERR_AUTO'), 
                                              N, tiles, tinds) > 5
                        )
        
        
        mask_size_plus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'FLUX_RADIUS'), 
                                              N, tiles, tinds) * 0.236 > 0.5
                         )
                
        
        mask_SEflag_minus = ( (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'FLAGS'), 
                                                    N, tiles, tinds) <= 3) & 
                             (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'i', i, mode = 'minus'), columns = 'FLAGS'), 
                                                    N, tiles, tinds) <= 3) & 
                             (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'z', i, mode = 'minus'), columns = 'FLAGS'), 
                                                    N, tiles, tinds) <= 3)
                           )
        
        mask_snr_minus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'FLUX_AUTO'), 
                                              N, tiles, tinds) /
                         parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'FLUXERR_AUTO'), 
                                              N, tiles, tinds) > 5
                        )
        
        
        mask_size_minus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'FLUX_RADIUS'), 
                                              N, tiles, tinds) * 0.236 > 0.5
                         )
        
        
        number_plus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'NUMBER') + i*int(1e6), 
                                              N, tiles, tinds)
                         )
        
        
        number_minus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'NUMBER') + i*int(1e6), 
                                              N, tiles, tinds)
                         )

        
        print("Size before cuts", len(number_plus))
        
        ra_plus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'ALPHA_J2000'), 
                                              N, tiles, tinds)
                  )
        
        ra_minus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'ALPHA_J2000'), 
                                              N, tiles, tinds)
                  )
        
        dec_plus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'plus'), columns = 'DELTA_J2000'), 
                                              N, tiles, tinds)
                  )
        
        dec_minus = (parallel_concatenator(lambda t,i: fitsio.read(_get_cat_path(t, 'r', i, mode = 'minus'), columns = 'DELTA_J2000'), 
                                              N, tiles, tinds)
                  )
        
        area_mask_plus  = ((Foreground_map == 0) & (Badcolor_map == 0))[hp.ang2pix(4096, ra_plus,  dec_plus,  lonlat = True)]
        area_mask_minus = ((Foreground_map == 0) & (Badcolor_map == 0))[hp.ang2pix(4096, ra_minus, dec_minus, lonlat = True)]
        
        
    ngrid = 2
    
    spacing = 10_000/ngrid
    Results = {}
    
    for g_name, g, pos in zip(['p', 'm'], [gplus, gminus], [[ra_plus, dec_plus], [ra_minus, dec_minus]]):
        
        n = 'gplus' if g_name == 'p'  else 'gminus'
        
        Flags = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_flags'), N, tiles, tinds)
        ID    = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'id') + i*int(1e6), N, tiles, tinds)

        if 'truedet' not in config['gal_kws']['truth_type']:
                
            if g_name == 'p':
                inds    = np.intersect1d(ID, number_plus, return_indices = True)[1:]
                SE_mask = (mask_SEflag_plus & mask_snr_plus & mask_size_plus & area_mask_plus)[inds[1]][np.argsort(inds[0])]
                assert np.allclose(ID, number_plus[inds[1]][np.argsort(inds[0])] ), "ID is not matched properly"
            elif g_name == 'm':
                inds    = np.intersect1d(ID, number_minus, return_indices = True)[1:]
                SE_mask = (mask_SEflag_minus & mask_snr_minus & mask_size_minus & area_mask_minus)[inds[1]][np.argsort(inds[0])]
                assert np.allclose(ID, number_minus[inds[1]][np.argsort(inds[0])] ), "ID is not matched properly"
        
        A = EBV[hp.ang2pix(4096, pos[0], pos[1], lonlat = True)][inds[1]][np.argsort(inds[0])]
        
        
        for s in tqdm(['noshear', '1p', '1m', '2p', '2m']):
              
            
            SNR    = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_s2n_%s'%s), N, tiles, tinds)
            Tratio = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_T_ratio_%s'%s), N, tiles, tinds)
            g1, g2 = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_g_%s'%s), N, tiles, tinds).T
            T      = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_T_%s'%s), N, tiles, tinds)
            mag_r  = 30 - 2.5*np.log10(parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_flux_%s'%s)[:, 0], N, tiles, tinds))
            
            mag_r -= A * 2.140

            Mask = (Flags == 0) & (SNR > 10) & (SNR < 1000) & (T < 10) & (Tratio > 0.5) & np.invert((T > 2) & (SNR < 30)) & np.invert((np.log10(T) < (22.25 - 0)/3.5) & (g1**2 + g2**2 > 0.8**2))
            
            Mask = Mask & SE_mask
            
            
            #######################################
            #Now setup the SOM classifications
            #######################################
            
            cell_path = PATH + '/wide_cells_%s_%s_%s.npy'%(t, n, s)
            tomobins  = np.load('/project/chihway/dhayaa/DECADE/SOMPZ/Runs/20240408/TomoBinAssign.npy')
                
            
            if os.path.isfile(cell_path) & (args['NewClassification'] == False):
                
                wide_cells = np.load(cell_path).astype(int)
                wide_bins  = tomobins[wide_cells.astype(int)]
                
                assert wide_bins.size == np.sum(Mask), f"Classified {wide_bins.size} out of {np.sum(Mask)} needed objects"
                
            
            else:
                
                flux  = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_flux_%s'%s), N, tiles, tinds)
                flux_err = parallel_concatenator(lambda t,i: fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'mcal_flux_err_%s'%s), N, tiles, tinds)

                flux[:, 0] *= 10**(0.4 * A * 2.140)
                flux[:, 1] *= 10**(0.4 * A * 1.569)
                flux[:, 2] *= 10**(0.4 * A * 1.196)

                flux_err[:, 0] *= 10**(0.4 * A * 2.140)
                flux_err[:, 1] *= 10**(0.4 * A * 1.569)
                flux_err[:, 2] *= 10**(0.4 * A * 1.196)
                
                flux     = flux[Mask]
                flux_err = flux_err[Mask]

                weight_path = '/project/chihway/dhayaa/DECADE/SOMPZ/Runs/20240408/WIDE_SOM_weights.npy'
                som = Classifier(som_weights = np.load(weight_path)) #Build classifier first, as its faster
                Nproc = np.max([os.cpu_count(), flux.shape[0]//100_000 + 1])
                inds  = np.array_split(np.arange(flux.shape[0]), Nproc)

                #Temp func to run joblib
                def _func_(i):

                    with SuppressPrint():
                        cell_id_i = som.classify(flux[inds[i], :], flux_err[inds[i], :])

                    return i, cell_id_i

                with joblib.parallel_backend("loky"):
                    jobs    = [joblib.delayed(_func_)(i) for i in range(Nproc)]
                    outputs = joblib.Parallel(n_jobs = 30, verbose=10,)(jobs)

                    cell_id = np.zeros(flux.shape[0])
                    for o in outputs: cell_id[inds[o[0]]] = o[1][0]

#                 np.save(cell_path, cell_id)
                wide_bins = tomobins[cell_id.astype(int)]
            
            
            Results['e1_%s_%s'%(g_name, s)] = g1[Mask]
            Results['e2_%s_%s'%(g_name, s)] = g2[Mask]
                       
            Results['bin_%s_%s'%(g_name, s)] = wide_bins #Already masked since we only pass in masked fluxes
            
            
            def batch_num_maker(t, i):
                
                x = fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'x')
                y = fitsio.read(PATH + '/metacal_%s_%s.fits'%(t, n), columns = 'y')
                
                index = x//spacing + y//spacing * ngrid + i * ngrid**2
                
                return index
                
            Results['batch_num_%s_%s'%(g_name, s)] = parallel_concatenator(batch_num_maker, N, tiles, tinds).astype(int)[Mask]


    Npatch = ngrid**2*len(tilenames)
    N = np.arange(Npatch)
#     N = np.random.choice(N, 200, replace = False)
    e1_plus  = np.zeros(N.size)
    e1_minus = np.zeros(N.size)
    e2_plus  = np.zeros(N.size)
    e2_minus = np.zeros(N.size)

    R11_plus  = np.zeros(N.size)
    R11_minus = np.zeros(N.size)
    R22_plus  = np.zeros(N.size)
    R22_minus = np.zeros(N.size)

    for b_i, b in enumerate([None, 0, 1, 2, 3]):
        
        def _func_(i):
            
            e1_plus  = summary(_get_subset_array(Results, N[i], 'e1_p_noshear', b))
            e1_minus = summary(_get_subset_array(Results, N[i], 'e1_m_noshear', b))
            e2_plus  = summary(_get_subset_array(Results, N[i], 'e2_p_noshear', b))
            e2_minus = summary(_get_subset_array(Results, N[i], 'e2_m_noshear', b))

            R11_plus  = ( summary( _get_subset_array(Results, N[i], 'e1_p_1p', b) ) - 
                          summary( _get_subset_array(Results, N[i], 'e1_p_1m', b) ) )/0.02
            R11_minus = ( summary( _get_subset_array(Results, N[i], 'e1_m_1p', b) ) - 
                          summary( _get_subset_array(Results, N[i], 'e1_m_1m', b) ) )/0.02
            R22_plus  = ( summary( _get_subset_array(Results, N[i], 'e2_p_2p', b) ) - 
                          summary( _get_subset_array(Results, N[i], 'e2_p_2m', b) ) )/0.02
            R22_minus = ( summary( _get_subset_array(Results, N[i], 'e2_m_2p', b) ) - 
                          summary( _get_subset_array(Results, N[i], 'e2_m_2m', b) ) )/0.02
            
            return i, (e1_plus, e1_minus, e2_plus, e2_minus), (R11_plus, R11_minus, R22_plus, R22_minus)
        
        
        with joblib.parallel_backend("loky"):
            jobs    = [joblib.delayed(_func_)(i) for i in range(Npatch)]
            outputs = joblib.Parallel(n_jobs = 40, verbose=10,)(jobs)
            for o in outputs: 
                e1_plus[o[0]]  = o[1][0]
                e1_minus[o[0]] = o[1][1]
                e2_plus[o[0]]  = o[1][2]
                e2_minus[o[0]] = o[1][3]
                
                R11_plus[o[0]]  = o[2][0]
                R11_minus[o[0]] = o[2][1]
                R22_plus[o[0]]  = o[2][2]
                R22_minus[o[0]] = o[2][3]

        R11 = (R11_plus + R11_minus)/2
        R22 = (R22_plus + R22_minus)/2

        g1 = (e1_plus - e1_minus)/(2*R11)
        g2 = (e2_plus + e2_minus)/(2*R22)

        m = g1/0.02 - 1
        c = g2

        if b == None: b = -99
        print('\n---------------------------------\n')
        print("----------- COMBINED (bin %d)-------" % b)
        print('size      N = %d (Masked %0.2f)'%(np.sum(Results['bin_p_noshear'] == b), np.average(Mask)))
        print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
        print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))

        print('\n---------------------------------\n')

        with open(os.path.dirname(__file__) + '/README.md', 'w') as f:

            f.write("TEST RESULTS FOR %s\n" %name)
            f.write("---------------------------------\n")
            f.write("size      N = %d (Masked %0.2f)\n"%(np.sum(Results['bin_p_noshear'] == b), np.average(Mask)))
            f.write("recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]\n"%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
            f.write("recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]\n"%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))

        R11 = R11_plus
        R22 = R22_plus

        g1 = e1_plus/R11
        g2 = e2_plus/R22

        m = g1/0.02 - 1
        c = g2

        print("----------- PLUS ---------------")
        print('size      N = %d'%Results['e1_p_noshear'].shape[0])
        print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
        print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))

        R11 = R11_minus
        R22 = R22_minus

        g1 = e1_minus/R11
        g2 = e2_minus/R22

        m = g1/(-0.02) - 1
        c = g2

        print("----------- MINUS ---------------")
        print('size      N = %d'%Results['e1_p_noshear'].shape[0])
        print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)*np.sqrt(Npatch)))
        print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)*np.sqrt(Npatch)))


        print("------------------------------------------------------------")
