#!/bin/bash -l
#SBATCH -t 10:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=50M
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1_SJKFkeep.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1_SJKFkeep.err

module load plink

plink --file /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1 --allow-extra-chr --dog --keep /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Keep_SJKF157.txt --recode --make-bed --out /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1_SJKFkeep
