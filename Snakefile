samples = [
    "NC_026433.1",
    "PP151906",
    "PP151907",
    "PP151908",
    "PP151909",
    "PP151910",
    "PP151911",
    "PP151912",
    "PP151913",
    "PP151914",
    "PP151915",
    "PP151916",
    "PP151917",
    "PP151918",
    "PP151919",
    "PP151920",
    "PP151921",
    "PP151922",
    "PP151923",
    "PP151924",
    "PP151925"
    ]

rule all:
    input:
        "intermediate/tree/virus_aligned.treefile"


rule download_data:
    output:
        "data/genomes/{sample}.fasta"
    shell:
        ' wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={wildcards.sample}&db=nuccore&report=fasta&retmode=text" -O {output}'

rule concatenate_genomes:
    input:
        expand("data/genomes/{sample}.fasta",sample=samples)
    output:
        "intermediate/genomes/virus_genomes.fasta"
    shell:
        'cat {input} > {output}'

rule align_genomes:
    input:
        rules.concatenate_genomes.output
    output:
        "intermediate/aligned/virus_aligned.fasta"
    conda:
        "yaml/mafft_env.yaml"
    shell:
        'mafft --auto {input} > {output}'

rule lower_to_upper_nucleotides:
    input:
        rules.align_genomes.output
    output:
        "intermediate/aligned/virus_aligned_clean.fasta"
    shell:
        'cat {input} | tr [:lower:] [:upper:] > {output}'


rule convert_fasta_to_phy:
    input:
        rules.lower_to_upper_nucleotides.output
    output:
        "intermediate/aligned/virus_aligned.phy"
    conda:
        "yaml/emboss_env.yaml"
    shell:
        'seqret -sequence {input} -outseq {output} -osformat2 phylip'

rule clean_phylip_file:
    input:
        rules.convert_fasta_to_phy.output
    output:
        "intermediate/aligned/virus_aligned_clean.phy"
    shell:
        "sed -E 's/\.[0-9]+/  /g; s/\./ /g' {input} > {output}"

rule maximum_likelihood_tree:
    input:
        rules.clean_phylip_file.output
    output:
        "intermediate/tree/virus_aligned.treefile"
    conda:
        "yaml/iqtree_env.yaml"
    shell:
        '''
        mkdir -p intermediate/tree
        iqtree2 -s {input} --prefix "intermediate/tree/virus_aligned"
        '''


