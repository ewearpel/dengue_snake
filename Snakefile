configfile: "config.yaml"

import os

accession_file = config["accessionfile"]
prefix = config["output_prefix"]
work = config.get("working_directory", os.getcwd())
samples = [line.strip() for line in open(accession_file)]
threads = config["threads"]

rule all:
    input:
        "results/tree/virus_aligned.treefile"

rule download_data:
    output:
        "{work}/data/genomes/{sample}.fasta"
    shell:
        ' wget -q "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={wildcards.sample}&db=nuccore&report=fasta&retmode=text" -O {output}'

rule concatenate_genomes:
    input:
        expand("{work}/data/genomes/{sample}.fasta",sample=samples)
    output:
        "{work}/intermediate/genomes/{prefix}_genomes.fasta"
    shell:
        'cat {input} > {output}'

rule align_genomes:
    input:
        rules.concatenate_genomes.output
    output:
        "{work}/intermediate/aligned/{prefix}_aligned.fasta"
    conda:
        "yaml/mafft_env.yaml"
    shell:
        'mafft --auto {input} > {output}'

rule lower_to_upper_nucleotides:
    input:
        rules.align_genomes.output
    output:
        "{work}/intermediate/aligned/{prefix}_aligned_clean.fasta"
    shell:
        'cat {input} | tr [:lower:] [:upper:] > {output}'


rule convert_fasta_to_phy:
    input:
        rules.lower_to_upper_nucleotides.output
    output:
        "{work}/intermediate/aligned/{prefix}_aligned.phy"
    conda:
        "yaml/emboss_env.yaml"
    shell:
        'seqret -sequence {input} -outseq {output} -osformat2 phylip'

rule clean_phylip_file:
    input:
        rules.convert_fasta_to_phy.output
    output:
        "{work}/intermediate/aligned/{prefix}_aligned_clean.phy"
    shell:
        "sed -E 's/\.[0-9]+/  /g; s/\./ /g' {input} > {output}"

rule maximum_likelihood_tree:
    input:
        rules.clean_phylip_file.output
    output:
        "{work}/results/tree/{prefix}_aligned.treefile"
    conda:
        "yaml/iqtree_env.yaml"
    params:
        threads=threads
    shell:
        '''
        [ -d {wildcards.work}/results/tree ] && rm -rf {wildcards.work}/results/tree
        mkdir -p {wildcards.work}/results/tree
        iqtree2 -s {input} --prefix "{work}/results/tree/{prefix}_aligned" -nt {params.threads}
        '''


