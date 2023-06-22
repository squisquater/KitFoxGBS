#!/bin/bash -l
#SBATCH -t 5:00:00
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH --mem=10G
#SBATCH -o /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmscripts/Stacks_Populations_20230128_169KF_DogAlign_SNPwhitelist15066_Subspp.out
#SBATCH -e /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/slurmscripts/Stacks_Populations_20230128_169KF_DogAlign_SNPwhitelist15066_Subspp.err

module load stacks

populations -P /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF235_CanFamAlign_SE/ -O /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF169_CanFamAlign_SE/Stacks_Populations_20230128_169KF_DogAlign_SNPwhitelist15066_Subspp -M /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/popmaps/PopMapSubspp_KF169_maf0.02moh0.6_mind0.7_geno0.2.txt \
-W /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/stacks_outputs/MergedRuns_KF169_CanFamAlign_SE/Plink_20230128_235KF_DogAlign_maf0.02moh0.6_mind0.5_geno0.2_SNPwhitelist.txt \
--fstats -k --smooth --sigma 150000 --bootstrap \
--structure --vcf --plink --radpainter --genepop --phylip-var --treemix
