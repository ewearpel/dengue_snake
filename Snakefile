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
        "intermediate/aligned/dengue_aligned.fasta"

rule download_data:
    output:
        "intermediate/genomes/{sample}.fasta"
    shell:
        ' wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={wildcards.sample}&db=nuccore&report=fasta&retmode=text" -O {output}'

rule concatenate_genomes:
    input:
        expand("intermediate/genomes/{sample}.fasta",sample=samples)
    output:
        "intermediate/genomes/dengue_genomes.fasta"
    shell:
        'cat {input} > {output}'

rule align_genomes:
    input:
        rules.concatenate_genomes.output
    output:
        "intermediate/aligned/dengue_aligned.fasta"
    conda:
        "mafft.yaml"
    shell:
        'mafft --auto {input} > {output}'

