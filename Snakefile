SPECIES = ["mCanLor1", "mCanLor2"]

rule all:
    input:
        expand("aligned_bam/{species}/{species}_merged_sorted.bam", species = SPECIES),
        expand("aligned_bam/{species}/{species}_merged_sorted.bam.bai", species = SPECIES),
        expand("errorStats/{species}/IGH.txt", species = SPECIES),
        expand("errorStats/{species}/IGK.txt", species = SPECIES),
        expand("errorStats/{species}/IGL.txt", species = SPECIES),
        expand("errorStats/{species}/nonIG.txt", species = SPECIES)

rule dataPrepAutomate:
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        "aligned_bam/{species}/{species}_merged_sorted.bam.bai",
        "aligned_bam/{species}/{species}_merged_sorted.bam"
    params:
        species = "{species}"
    shell:
        """
        sbatch --partition=qcb {input.script} -s {params.species}
        """

rule runIGdetective:
    input:
        script = "code/igDetective.sh",
        pri_genome = "/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.pri.fasta",
        alt_genome = "/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.alt.fasta"
    output:
        pri_out = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt",
        alt_our = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt"
    params:
        pri_outdir = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/",
        alt_ourdir = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/"
    shell:
        """
        sbatch --partition=qcb {input.script} {input.pri_genome} {params.pri_outdir}
        sbatch --partition=qcb {input.script} {input.alt_genome} {params.alt_ourdir}
        """
       
rule geneLociAutomate:
    input:
        "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt",
        "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt",
        script = "code/geneLociAutomated.sh"
    output:
        "gene_position/{species}/alt/{species}_IGL_pos.sorted.bed",  
        "gene_position/{species}/pri/{species}_IGL_pos.sorted.bed"
    params:
        species = "{species}"
    shell:
        """
        {input.script} -s {params.species} -g pri
        {input.script} -s {params.species} -g alt
        """

rule mergeLoci:
    input:
        pri = "gene_position/{species}/pri/{species}_IGL_pos.sorted.bed",
        alt = "gene_position/{species}/alt/{species}_IGL_pos.sorted.bed",
        script = "code/mergePriAlt.sh"
    output:
        "gene_position/{species}/combined/{species}_IGL_pos.sorted.bed",
        "gene_position/{species}/combined/{species}_IGH_pos.sorted.bed",
        "gene_position/{species}/combined/{species}_IGK_pos.sorted.bed"
    params:
        species = "{species}"
    shell:
        """
        {input.script} -s {params.species}
        """

rule convertPrimaryBam:
    input:
        bam = "aligned_bam/{species}/{species}_merged_sorted.bam"
    output:
        outputbam = "aligned_bam/{species}/{species}_merged_sorted_primary.bam"
    shell:
        """
        source /etc/profile.d/modules.sh
        module load gcc/11.3.0
        module load samtools/1.17
        samtools view -b -F 0x800 -F 0x100 -@ 30 {input.bam} > {output.outputbam}
        """

rule cigarProcessing:
    input:
        script = "code/cigar_processing.py",
        bam = "aligned_bam/{species}/{species}_merged_sorted_primary.bam", 
        IGHbed = "gene_position/{species}/combined/{species}_IGH_pos.sorted.bed",
        IGKbed = "gene_position/{species}/combined/{species}_IGK_pos.sorted.bed", 
        IGLbed = "gene_position/{species}/combined/{species}_IGL_pos.sorted.bed"
    output:
        "errorStats/{species}/IGH.txt",
        "errorStats/{species}/IGK.txt",
        "errorStats/{species}/IGL.txt",
        "errorStats/{species}/nonIG.txt"
    params:
        species = "{species}"
    shell:
        """
        python {input.script} {input.bam} {input.IGHbed} {input.IGKbed} {input.IGLbed} {params.species}
        """

#####TODOs

#rule stats analysis:
#rule qv0: