rule process_erc:
    input:
        ercc_fasta = config["reference"]['external_sequences'],
        reference_transcriptome = config["reference"]['transcriptome_fasta'],
        reference_genome = config["reference"]['genome_fasta'],
        reference_gtf = config["reference"]['reference_gtf'],
        reference_gff3 = config["reference"]['reference_gff3']
    output:
        ercc_gtf_out = 'data/references/ERCC.gff3',
        reference_gtf_out = temp('data/references/gencode.v45.primary_assembly.annotation_copy.gtf'),
        reference_gff3_out = temp('data/references/gencode.v45.primary_assembly.annotation_copy.gff3'),
        reference_transcriptom_out = temp('data/references/gencode.v45.transcripts_copy.fa'),
        reference_genome_out = temp('data/references/GRCh38.primary_assembly.genome_copy.fa')
    conda:
        'preprocess_ercc.yaml'
    script:
        'preprocess_ercc.py'

