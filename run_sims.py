#!/usr/bin/env python

#copied directly from https://github.com/beckermr/misc/blob/68b42ea226913b97236ce5145378c6cf6fe05b8b/bin/run-simple-des-y3-sim
import logging
import sys

import click
import yaml

from band_infoing import make_band_info
#from matts_misc.simple_des_y3_sims.simulating import End2EndSimulation
#from matts_misc.simple_des_y3_sims.true_detecting import make_true_detections
#from matts_misc.simple_des_y3_sims.medsing import make_meds_files
#from matts_misc.simple_des_y3_sims.run_metacal import run_metacal

for lib in ['matts_misc.simple_des_y3_sims']:
    lgr = logging.getLogger(lib)
    hdr = logging.StreamHandler(sys.stdout)
    hdr.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    lgr.setLevel(logging.DEBUG)
    lgr.addHandler(hdr)


@click.group()
def cli():
    """Run simple DES Y3 end-to-emd simulations."""
    pass


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to prep for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--n-files', type=int, default=None,
              help='number of SE images to keep - useful for testing')
def prep(tilename, bands, output_desdata, n_files):
    """Prepare a tile for simulating."""
    make_band_info(
        tilename=tilename,
        bands=[b for b in bands],
        output_meds_dir=output_desdata,
        n_files=n_files)
'''

@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--seed', type=int, required=True,
              help='the base RNG seed')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def galsim(tilename, bands, output_desdata, seed, config_file):
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    sim = End2EndSimulation(
        seed=seed,
        output_meds_dir=output_desdata,
        tilename=tilename,
        bands=[b for b in bands],
        gal_kws=config['gal_kws'],
        psf_kws=config['psf_kws'])
    sim.run()

#STOPPING here since we want to get the 2 first steps above done at first

@cli.command('true-detection')
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to run detection for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
def true_detection(tilename, bands, output_desdata, config_file):
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    make_true_detections(
        tilename=tilename,
        bands=[b for b in bands],
        output_meds_dir=output_desdata,
        box_size=config['true_detection']['box_size'])


@cli.command('swarp')
def swarp():
    click.echo('run swarp for a tile')


@cli.command('source-extractor')
def source_extractor():
    click.echo('run source-extractor for a tile')


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to make MEDS files for as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--config-file', type=str, required=True,
              help='the YAML config file')
@click.option('--meds-config-file', type=str, required=True,
              help='the YAML config file for MEDS making')
def meds(tilename, bands, output_desdata, config_file, meds_config_file):
    with open(config_file, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.Loader)
    with open(meds_config_file, 'r') as fp:
        meds_config = yaml.load(fp, Loader=yaml.Loader)
    make_meds_files(
        tilename=tilename,
        bands=[b for b in bands],
        output_meds_dir=output_desdata,
        psf_kws=config['psf_kws'],
        meds_config=meds_config)


@cli.command()
@click.option('--tilename', type=str, required=True,
              help='the coadd tile to simulate')
@click.option('--bands', type=str, required=True,
              help=('a list of bands to simulate as '
                    'a concatnated string (e.g., "riz")'))
@click.option('--output-desdata', type=str, required=True,
              help='the output DESDATA directory')
@click.option('--seed', type=int, required=True,
              help='the base RNG seed')
def metacal(tilename, bands, output_desdata, seed):
    run_metacal(
        tilename=tilename,
        output_meds_dir=output_desdata,
        bands=[b for b in bands],
        seed=seed)
'''

if __name__ == '__main__':
    cli()
