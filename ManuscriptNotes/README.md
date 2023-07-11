## Read Processing and Locus Assembly

Raw paired-end reads were demultiplexed using the 'process_radtags' function in Stacks v. 1.30 [cite] using the options -q (discard low quality reads) and -c (remove reads with uncalled bases); this program was used because it could handle combinatorial barcodes. A schematic of the bioinformatic methods is provided in Fig 2 for clarity. To check for and merge overlapping reads, the demultiplexed reads were run through PEAR v. 0.9.8 [54] using default settings and a 16 bp minimum overlap.

Stacks demultiplexing and filtering retained on average 95.6% (93.1–97.1%) of reads per sample. <- this stat includes some non kit fox samples.
If I actually want this data I need to extract it from the log files.
GBS1: 68.8%
GBS2: 94.1%
GBS3: 94.7%
GBS4: 97.1%
GBS5: 
GBS6: 96.1%
GBS7: 97.0%

PEAR merged on average 56.2% (44.6–68.4%) of demultiplexed reads.
