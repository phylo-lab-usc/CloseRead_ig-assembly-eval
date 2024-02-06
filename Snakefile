SPECIES = ["mCanLor1"]
CORRECTED  = ["True"]


files = dict()
species = SPECIES[0]
if(CORRECTED[0] == "True"):
    files["bam"] = "aligned_bam/{species}.corrected/{species}_merged_sorted.bam"
    files["bai"] = "aligned_bam/{species}.corrected/{species}_merged_sorted.bam.bai"
    files['igDetec_pri'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/combined_genes_IGL.txt"
    files['igDetec_alt'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/combined_genes_IGK.txt"
    files['geneOut_pri'] = "gene_position/{species}/lja/{species}_IGL_pos.sorted.bed"
    files['geneOut_alt'] = "gene_position/{species}/lja/{species}_IGH_pos.sorted.bed"
    files['geneOut_alt2'] = "gene_position/{species}/lja/{species}_IGK_pos.sorted.bed"
    files['outputbam'] = "aligned_bam/{species}.corrected/{species}_merged_sorted_primary.bam"
    files['IGHbed'] = "gene_position/{species}/lja/{species}_IGH_pos.sorted.bed"
    files['IGKbed'] = "gene_position/{species}/lja/{species}_IGK_pos.sorted.bed" 
    files['IGLbed'] = "gene_position/{species}/lja/{species}_IGL_pos.sorted.bed"
    files['final_IGH'] = "errorStats/{species}.lja/IGH.txt"
    files['final_IGK'] = "errorStats/{species}.lja/IGK.txt"
    files['final_IGL'] = "errorStats/{species}.lja/IGL.txt"
    files['final_nonIG'] = "errorStats/{species}.lja/nonIG.txt"
else:
    files["bam"] = "aligned_bam/{species}/{species}_merged_sorted.bam"
    files["bai"] = "aligned_bam/{species}/{species}_merged_sorted.bam.bai"
    files['igDetec_pri'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt"
    files['igDetec_alt'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt"
    files['geneOut_pri'] = "gene_position/{species}/pri/{species}_IGL_pos.sorted.bed"
    files['geneOut_alt'] = "gene_position/{species}/alt/{species}_IGL_pos.sorted.bed"
    files['geneOut_alt2'] = "gene_position/{species}/alt/{species}_IGK_pos.sorted.bed"
    files['outputbam'] = "aligned_bam/{species}/{species}_merged_sorted_primary.bam"
    files['IGHbed'] = "gene_position/{species}/combined/{species}_IGH_pos.sorted.bed"
    files['IGKbed'] = "gene_position/{species}/combined/{species}_IGK_pos.sorted.bed" 
    files['IGLbed'] = "gene_position/{species}/combined/{species}_IGL_pos.sorted.bed"
    files['final_IGH'] = "errorStats/{species}/IGH.txt"
    files['final_IGK'] = "errorStats/{species}/IGK.txt"
    files['final_IGL'] = "errorStats/{species}/IGL.txt"
    files['final_nonIG'] = "errorStats/{species}/nonIG.txt"


rule all:
    input:
        expand(files['final_IGH'], species = SPECIES, corrected = CORRECTED),
        expand(files['final_IGK'], species = SPECIES, corrected = CORRECTED),
        expand(files['final_IGL'], species = SPECIES, corrected = CORRECTED),
        expand(files['final_nonIG'], species = SPECIES, corrected = CORRECTED)


rule dataPrepAutomate:
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        files["bam"],
        files["bai"]
    params:
        species = "{species}",
        corrected = CORRECTED[0]
    shell:
        """
        sbatch --partition=qcb {input.script} -s {params.species} -c {params.corrected}
        """

rule runIGdetective_noCorr:
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

rule runIGdetective_Corr:
    input:
        script = "code/igDetective.sh",
        corr_genome = "/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.ljacorr.merged.fasta"
    output:
        corr_out1 = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/combined_genes_IGL.txt",
        corr_out2 = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/combined_genes_IGK.txt"
    params:
        corr_outdir = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/"
    shell:
        """
        sbatch --partition=qcb {input.script} {input.corr_genome} {params.corr_outdir}
        """

rule geneLociAutomate:
    input:
        files['igDetec_pri'],
        files['igDetec_alt'],
        script = "code/geneLociAutomated.sh"
    output:
        files['geneOut_pri'],  
        files['geneOut_alt'],
        files['geneOut_alt2']
    params:
        species = "{species}",
        corrected = CORRECTED[0]
    shell:
        """
        if [ {params.corrected} == "True" ]; then
            {input.script} -s {params.species} -g lja
        else
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
        bam = files["bam"]
    output:
        outputbam = files['outputbam']
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
        bam = files['outputbam'], 
        IGHbed = files['IGHbed'],
        IGKbed = files['IGKbed'], 
        IGLbed = files['IGLbed']
    output:
        files['final_IGH'],
        files['final_IGK'],
        files['final_IGL'],
        files['final_nonIG']
    params:
        species = "{species}"
    shell:
        """
        python {input.script} {input.bam} {input.IGHbed} {input.IGKbed} {input.IGLbed} {params.species}
        """

#####TODOs

#rule stats analysis:
#rule qv0: