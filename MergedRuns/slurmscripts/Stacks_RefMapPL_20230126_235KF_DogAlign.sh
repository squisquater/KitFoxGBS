#!/bin/bash -l
#SBATCH -t 2:00:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=10G
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Stacks_RefMapPL_20230126_235KF_DogAlign.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Stacks_RefMapPL_20230126_235KF_DogAlign.err

module load stacks

ref_map.pl -T 1 --samples /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/bamfiles/CanFamAlign_SE --popmap /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/popmaps/KF235popmap.txt --unpaired -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE
