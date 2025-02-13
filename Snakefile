configfile: "config.yaml"

accession_file = config["accessionfile"]
prefix = config["output_prefix"]
samples = [line.strip() for line in open(accession_file)]
threads = config["threads"]

rule all:
    input:
        "results/tree/virus_aligned.treefile"

rule download_data:
    output:
        "data/genomes/{sample}.fasta"
    shell:
        ' wget -q "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={wildcards.sample}&db=nuccore&report=fasta&retmode=text" -O {output}'

rule concatenate_genomes:
    input:
        expand("data/genomes/{sample}.fasta",sample=samples)
    output:
        "intermediate/genomes/{prefix}_genomes.fasta"
    shell:
        'cat {input} > {output}'

rule align_genomes:
    input:
        rules.concatenate_genomes.output
    output:
        "intermediate/aligned/{prefix}_aligned.fasta"
    conda:
        "yaml/mafft_env.yaml"
    shell:
        'mafft --auto {input} > {output}'

rule lower_to_upper_nucleotides:
    input:
        rules.align_genomes.output
    output:
        "intermediate/aligned/{prefix}_aligned_clean.fasta"
    shell:
        'cat {input} | tr [:lower:] [:upper:] > {output}'


rule convert_fasta_to_phy:
    input:
        rules.lower_to_upper_nucleotides.output
    output:
        "intermediate/aligned/{prefix}_aligned.phy"
    conda:
        "yaml/emboss_env.yaml"
    shell:
        'seqret -sequence {input} -outseq {output} -osformat2 phylip'

rule clean_phylip_file:
    input:
        rules.convert_fasta_to_phy.output
    output:
        "intermediate/aligned/{prefix}_aligned_clean.phy"
    shell:
        "sed -E 's/\.[0-9]+/  /g; s/\./ /g' {input} > {output}"

rule maximum_likelihood_tree:
    input:
        rules.clean_phylip_file.output
    output:
        "results/tree/{prefix}_aligned.treefile"
    conda:
        "yaml/iqtree_env.yaml"
    params:
        threads=threads
    shell:
        '''
        [ -d results/tree ] && rm -rf results/tree
        mkdir -p results/tree
        iqtree2 -s {input} --prefix "results/tree/{wildcards.prefix}_aligned" -nt {params.threads}
        '''


