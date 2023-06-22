#!/bin/bash -l
#SBATCH -t 10:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=50M
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmscripts/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.2.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmscripts/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.2.err

module load plink

plink --file /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7 --allow-extra-chr --dog --geno 0.2 --recode --out /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.2
