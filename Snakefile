SPECIES = "mCanLor2"

rule all:
    input:
        expand("aligned_bam/{species}/{species}_merged_sorted.bam.bai", species = SPECIES),
        expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt", species = SPECIES),
        expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt", species = SPECIES),  
        expand("gene_position/{species}/combined/{species}_IGL_pos.sorted.bed", species = SPECIES)

rule dataPrepAutomate:
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        expand("aligned_bam/{species}/{species}_merged_sorted.bam.bai", species = SPECIES)
    params:
        species = expand("{species}", species = SPECIES)
    shell:
        """
        sbatch --partition=qcb {input.script} -s {params.species}
        """

rule runIGdetective:
    input:
        script = "code/igDetective.sh",
        pri_genome = expand("/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.pri.fasta", species = SPECIES),
        alt_genome = expand("/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.alt.fasta", species = SPECIES)
    output:
        pri_out = expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt", species = SPECIES),
        alt_our = expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt", species = SPECIES)
    params:
        pri_outdir = expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/", species = SPECIES),
        alt_ourdir = expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/", species = SPECIES)
    shell:
        """
        sbatch --partition=qcb {input.script} {input.pri_genome} {params.pri_outdir}
        sbatch --partition=qcb {input.script} {input.alt_genome} {params.alt_ourdir}
        """
       
rule geneLociAutomate:
    input:
        script = "code/geneLociAutomated.sh",
    output:
        expand("gene_position/{species}/alt/{species}_IGL_pos.sorted.bed", species = SPECIES),  
        expand("gene_position/{species}/pri/{species}_IGL_pos.sorted.bed", species = SPECIES)
    params:
        species = expand("{species}", species = SPECIES)
    shell:
        """
        {input.script} -s {params.species} -g pri
        {input.script} -s {params.species} -g alt
        """

rule mergeLoci:
    input:
        pri = expand("gene_position/{species}/pri/{species}_IGL_pos.sorted.bed", species = SPECIES),
        alt = expand("gene_position/{species}/alt/{species}_IGL_pos.sorted.bed", species = SPECIES),
        script = "code/mergePriAlt.sh"
    output:
        expand("gene_position/{species}/combined/{species}_IGL_pos.sorted.bed", species = SPECIES)
    params:
        species = expand("{species}", species = SPECIES)
    shell:
        """
        {input.script} -s {params.species}
        """

#####TODOs

#rule mismatch:

#rule qv0: