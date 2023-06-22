#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem 8GB
#SBATCH -o "copybams.out"

cp /group/ctbrowngrp2/sophiepq/KitFoxGBS/GBS1_20220914/bamfiles/CanFamAlign_SE/*R1.sort.bam.bai  /group/ctbrowngrp2/sophiepq/KitFoxGBS/MergedRuns/bamfiles/CanFamAlign_SE/

