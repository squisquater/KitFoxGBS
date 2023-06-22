#!/bin/bash -l
#SBATCH -t 10:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=50M
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_missingness.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_missingness.err

module load plink

plink --file /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/populations.plink --allow-extra-chr --dog --recode --missing --out /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Plink_20230126_235KF_DogAlign_maf0.02moh0.6
