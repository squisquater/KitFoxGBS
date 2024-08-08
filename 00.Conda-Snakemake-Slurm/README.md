## Install Conda

I would recommend using [mamba](https://mamba.readthedocs.io/en/latest/) instead. It has all the same functionality but is **WAY** faster. I use micromamba which is a specific “flavor” of mamba but would recommend using just mamba. So basically I use conda/mamba /micromamba interchangeably below. Use the command associated with what you decide to install. They’ll all work.

**Conda Cheat Sheet:** [https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)

**Note:** Depending on how you want to run snakemake, you may have to make some mods to your .bashrc or .bash_profile to ensure you can activate conda from both the head node and the interactive nodes. Anytime you make changes to those files you need to log out of farm and log back in for them to take effect.

## Create & Activate the Conda Environment + Install Necessary Packages

To create a conda environment that will contain all the software you need to process your GBS data (including snakemake!) you can use this config file [GBS_config.yml](00.Conda-Snakemake-Slurm/GBS_config.yml)

Once you have this file you can create your conda environment
```
mamba env create -f GBS_config.yml
```
