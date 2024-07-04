'''
written by Sina
'''

from Bio import SeqIO
import shutil



print(snakemake.input['ercc_fasta'])
fasta_sequences = SeqIO.parse(open(snakemake.input['ercc_fasta']),'fasta')

##convert fasta to gtf file 
with open(snakemake.output['ercc_gtf_out'], 'w') as gff3:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        tags = name.split('|')
        gff3.write('\n{}\t'.format(tags[0])) ##seqid
        gff3.write('spike_in\tgene\t1\t')
        gff3.write(str(len(sequence))) ##end position
        gff3.write('\t.\t+\t.\t')
        gff3.write('gene_id={};gene_name={};transcript_id={}'.format(tags[1], tags[0], tags[0]))

##append fasta to referenc transcriptom fasta

shutil.copyfile(snakemake.input['reference_transcriptome'], snakemake.output['reference_transcriptom_out'])
with open(snakemake.output['reference_transcriptom_out'], 'a') as tr_ref:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        print(name)
        tr_ref.write('\n' + name + '\n')
        tr_ref.write(sequence)

##append fasta to reference genome fasta
shutil.copyfile(snakemake.input['reference_genome'], snakemake.output['reference_genome_out'])
with open(snakemake.output['reference_genome_out'], 'a') as tr_ref:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        print(name)
        tr_ref.write('\n' + name + '\n')
        tr_ref.write(sequence)

##append ercc gtf to reference gff3
shutil.copyfile(snakemake.input['reference_gff3'], snakemake.output['reference_gff3_out'])
with open(snakemake.output['reference_gff3_out'], 'a') as tr_ref:
    with open(snakemake.output['ercc_gtf_out'], 'r') as tr_ercc:
        tr_ref.write('\n' + tr_ercc.read())

##append ercc gtf to reference gtf
shutil.copyfile(snakemake.input['reference_gtf'], snakemake.output['reference_gtf_out'])
with open(snakemake.output['reference_gtf_out'], 'a') as tr_ref:
    with open(snakemake.output['ercc_gtf_out'], 'r') as tr_ercc:
        tr_ref.write('\n' + tr_ercc.read())