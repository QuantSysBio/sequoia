from glob import glob
from os.path import join

SAMPLES, = glob_wildcards(join(reads_fq, '{sample}_1.fq.gz'))
PATTERN_R1 = '{sample}_1.fq.gz'
PATTERN_R2 = '{sample}_2.fq.gz'

PATTERN_R1_trimmed = '{sample}_1.fq.gz'
PATTERN_R2_trimmed = '{sample}_2.fq.gz'
