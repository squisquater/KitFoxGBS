#!/bin/bash -l
#SBATCH -t 10:00:00
#SBATCH -p high
#SBATCH --mem=5G
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Admixture_20230204_SJKF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Admixture_20230204_SJKF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1.err

for K in {1,2,3,4,5,6,7,8,9,10}; do ~/admixture_linux-1.3.0/admixture -B2000 --cv=10 /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Plink_20230126_235KF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1_SJKFkeep.bed $K | tee /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6/Admixture/Admixture_20230204_SJKF_DogAlign_maf0.02moh0.6_mind0.7_geno0.1_K${K}.out

done
