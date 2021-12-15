# mcal_sim_test
Testing metacal measurements on simulated images

## Basic instructions:
1. create a conda environment using ```conda env create --file environment.yml```

2. activate the environment by doing ```conda activate shear``` (remember to export several needed directories, in particular TMPDIR for now)

3. that's probably not a complete set of all packages you need; in particular you should download desmeds from https://github.com/esheldon/desmeds (because the conda version doesn't have all we need); you might or might not need to install desmeds manually following the instructions in that link; you'll also need easyaccess and other simpler conda-installable packages

4. to run the basic prep stage, do for instance: ```python run_sims.py prep --tilename="DES0544-2249" --bands="i" --output-desdata="outputs-DES0544-2249"```

5. after running prep, do ```python run_sims.py galsim --tilename="DES0544-2249" --bands="i" --output-desdata="outputs-DES0544-2249" --seed="42" --config-file="./runs/v000_no_detection/config.yaml"```

6. then, ```python run_sims.py true-detection --tilename="DES0544-2249" --bands="i" --output-desdata="outputs-DES0544-2249"  --config-file="./runs/v000_no_detection/config.yaml"```
