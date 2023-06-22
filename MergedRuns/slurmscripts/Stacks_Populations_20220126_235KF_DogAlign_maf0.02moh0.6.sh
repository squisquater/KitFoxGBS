#!/bin/bash -l
#SBATCH -t 2:00:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=10G
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Stacks_Populations_20230126_235KF_DogAlign_maf0.02moh0.6.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmlogs/Stacks_Populations_20230126_235KF_DogAlign_maf0.02moh0.6.err

module load stacks

populations -P /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE -O /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/Stacks_Populations_20230126_235KF_maf0.02moh0.6 -M /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/popmaps/KF235popmap2.txt \
--fstats \
--min-maf 0.02 --max-obs-het 0.60 \
--structure --vcf --plink --radpainter --genepop --phylip-var --treemix \
