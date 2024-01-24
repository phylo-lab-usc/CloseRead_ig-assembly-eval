SPECIES = "mCanLor2"

rule all:
    input:
        expand("aligned_bam/{species}/{species}_merged_sorted.bam.bai", species = SPECIES),
        expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt", species = SPECIES),
        expand("/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt", species = SPECIES),  
        expand("gene_position/{species}/combined/{species}_IGL_pos.sorted.bed", species = SPECIES),
        expand("errorStats/{species}/IGH.txt", species = SPECIES),
        expand("errorStats/{species}/IGK.txt", species = SPECIES),
        expand("errorStats/{species}/IGL.txt", species = SPECIES),
        expand("errorStats/{species}/nonIG.txt", species = SPECIES)

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

rule convertPrimaryBam:
    input:
        bam = expand("aligned_bam/{species}/{species}_merged_sorted.bam", species = SPECIES)
    output:
        outputbam = expand("aligned_bam/{species}/{species}_merged_sorted_primary.bam", species = SPECIES)
    shell:
        """
        source /etc/profile.d/modules.sh
        module load gcc/11.3.0
        module load samtools/1.17
        samtools view -b -F 0x800 -F 0x100 -@ 60 {input.bam} > {output.outputbam}
        """

rule cigarProcessing:
    input:
        script = "code/cigar_processing.py",
        bam = expand("aligned_bam/{species}/{species}_merged_sorted_primary.bam", species = SPECIES),
        IGHbed = expand("gene_position/{species}/combined/{species}_IGH_pos.sorted.bed", species = SPECIES),
        IGKbed = expand("gene_position/{species}/combined/{species}_IGK_pos.sorted.bed", species = SPECIES),
        IGLbed = expand("gene_position/{species}/combined/{species}_IGL_pos.sorted.bed", species = SPECIES)
    output:
        expand("errorStats/{species}/IGH.txt", species = SPECIES),
        expand("errorStats/{species}/IGK.txt", species = SPECIES),
        expand("errorStats/{species}/IGL.txt", species = SPECIES),
        expand("errorStats/{species}/nonIG.txt", species = SPECIES)
    params:
        species = expand("{species}", species = SPECIES)
    shell:
        """
        python {input.script} {input.bam} {input.IGHbed} {input.IGKbed} {input.IGLbed} {params.species}
        """

#####TODOs

#rule stats analysis:
#rule qv0: