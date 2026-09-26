# Full pipeline: create the MMC, PCM and MC run folders, then run them all.

python mmc_runs_setup.py
python pcm_runs_setup.py
python preprocessing.py

sh run_pcm.sh
sh run_cases.sh
sh run_mc.sh