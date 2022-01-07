# mcal_sim_test
Testing metacal measurements on simulated images

## Basic instructions:
1. create a conda environment using ```conda env create --file environment.yml```

2. activate the environment by doing ```conda activate shear``` (remember to export several needed directories, in particular TMPDIR for now)

3. that's probably not a complete set of all packages you need; in particular you should download desmeds from https://github.com/esheldon/desmeds (because the conda version doesn't have all we need); you might or might not need to install desmeds manually following the instructions in that link; you'll also need easyaccess and other simpler conda-installable packages

4. to run the basic prep stage, do for instance: ```python run_sims.py prep --tilename="DES0544-2249" --bands="riz" --output-desdata="outputs-DES0544-2249"```

5. after running prep, do ```python run_sims.py galsim --tilename="DES0544-2249" --bands="riz" --output-desdata="outputs-DES0544-2249" --seed="42" --config-file="./runs/v000_no_detection/config.yaml"```

6. then, ```python run_sims.py true-detection --tilename="DES0544-2249" --bands="riz" --output-desdata="outputs-DES0544-2249"  --config-file="./runs/v000_no_detection/config.yaml"```

7. then, ```python run_sims.py meds --tilename="DES0544-2249" --bands="riz" --output-desdata="outputs-DES0544-2249"  --config-file="./runs/v000_no_detection/config.yaml" --meds-config-file="./runs/v000_no_detection/meds.yaml"```

8. home stretch: ```python run_sims.py metacal --tilename="DES0544-2249" --bands="riz" --output-desdata="outputs-DES0544-2249"  --seed="1"```


#### Some notes and caveats (from Lucas): 

0. For testing, it's way better to set `--bands="r"` above so you do only one band

1. Starting when I picked up the code again in January, the kicpaa job node was expecting ```python3``` as opposed to ```python``` in the command line

2. I've had weird problems when some files were already existing from previous runs. Make sure you delete everything in the directory set by `output-desdata`, which above is simply `outputs-DES0544-2249`, as well as everything in the directory you're exporting as `TMP_DIR` and also everything in `MEDS_DIR`

3. The `meds` step generates a lot of data (several gigabytes at least)


 
