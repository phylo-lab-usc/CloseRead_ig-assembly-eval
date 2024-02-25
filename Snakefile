SPECIES = ["mCynVol1"]
CORRECTED  = ["False"]
SOURCE = ["hifi_fastq"]

files = dict()
if(CORRECTED[0] == "True"):
    files["bam"] = "aligned_bam/{species}.corrected/{species}_merged_sorted.bam"
    files["bai"] = "aligned_bam/{species}.corrected/{species}_merged_sorted.bam.bai"
    files['igDetec_pri'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/combined_genes_IGL.txt"
    files['igDetec_alt'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.lja.igdetective/combined_genes_IGK.txt"
    files['geneOut_pri'] = "gene_position/{species}/lja/{species}_IGL_pos.sorted.bed"
    files['geneOut_alt'] = "gene_position/{species}/lja/{species}_IGH_pos.sorted.bed"
    files['geneOut_alt2'] = "gene_position/{species}/lja/{species}_IGK_pos.sorted.bed"
    files['outputbam'] = "aligned_bam/{species}.corrected/{species}_merged_sorted_primary.bam"
    files['outputbai'] = "aligned_bam/{species}.corrected/{species}_merged_sorted_primary.bam.bai"
    files['IGHbed'] = "gene_position/{species}/lja/{species}_IGH_pos.sorted.bed"
    files['IGKbed'] = "gene_position/{species}/lja/{species}_IGK_pos.sorted.bed" 
    files['IGLbed'] = "gene_position/{species}/lja/{species}_IGL_pos.sorted.bed"
    files['final_IGH'] = "errorStats/{species}.lja/IGH.txt"
    files['final_IGK'] = "errorStats/{species}.lja/IGK.txt"
    files['final_IGL'] = "errorStats/{species}.lja/IGL.txt"
    files['final_nonIG'] = "errorStats/{species}.lja/nonIG.txt"
    files['genome'] = "/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.ljacorr.merged.fasta"
    files['pileup_IGH'] = "errorStats/{species}.lja/IGH_alt_pileup.txt"
    files['pileup_IGK'] = "errorStats/{species}.lja/IGK_alt_pileup.txt"
    files['pileup_IGL'] = "errorStats/{species}.lja/IGL_alt_pileup.txt"
else:
    files["bam"] = "aligned_bam/{species}/{species}_merged_sorted.bam"
    files["bai"] = "aligned_bam/{species}/{species}_merged_sorted.bam.bai"
    files['igDetec_pri'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt"
    files['igDetec_alt'] = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt"
    files['geneOut_pri'] = "gene_position/{species}/pri/{species}_IGL_pos.sorted.bed"
    files['geneOut_alt'] = "gene_position/{species}/alt/{species}_IGL_pos.sorted.bed"
    files['geneOut_alt2'] = "gene_position/{species}/alt/{species}_IGK_pos.sorted.bed"
    files['outputbam'] = "aligned_bam/{species}/{species}_merged_sorted_primary.bam"
    files['outputbai'] = "aligned_bam/{species}/{species}_merged_sorted_primary.bam.bai"
    files['IGHbed'] = "gene_position/{species}/combined/{species}_IGH_pos.sorted.bed"
    files['IGKbed'] = "gene_position/{species}/combined/{species}_IGK_pos.sorted.bed" 
    files['IGLbed'] = "gene_position/{species}/combined/{species}_IGL_pos.sorted.bed"
    files['final_IGH'] = "errorStats/{species}/IGH.txt"
    files['final_IGK'] = "errorStats/{species}/IGK.txt"
    files['final_IGL'] = "errorStats/{species}/IGL.txt"
    files['final_nonIG'] = "errorStats/{species}/nonIG.txt"
    files['genome'] = "/home1/zhuyixin/sc1/AssmQuality/assemblies/{species}.merged.fasta"
    files['pileup_IGH'] = "errorStats/{species}/IGH_alt_pileup.txt"
    files['pileup_IGK'] = "errorStats/{species}/IGK_alt_pileup.txt"
    files['pileup_IGL'] = "errorStats/{species}/IGL_alt_pileup.txt"

rule all:
    input:
        expand(files['final_IGH'], species = SPECIES, corrected = CORRECTED),
        expand(files['final_IGK'], species = SPECIES, corrected = CORRECTED),
        expand(files['final_IGL'], species = SPECIES, corrected = CORRECTED),
        expand(files['final_nonIG'], species = SPECIES, corrected = CORRECTED),
        expand(files['pileup_IGH'], species = SPECIES, corrected = CORRECTED),
        expand(files['pileup_IGK'], species = SPECIES, corrected = CORRECTED),
        expand(files['pileup_IGL'], species = SPECIES, corrected = CORRECTED)


rule dataPrepAutomate:
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        files["bam"],
        files["bai"]
    params:
        species = "{species}",
        corrected = CORRECTED[0],
        source = SOURCE[0]
    shell:
        """
        sbatch --partition=qcb {input.script} -s {params.species} -c {params.corrected} -w {params.source}
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
        fi
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
        outputbam = files['outputbam'],
        outputbai = files['outputbai']
    shell:
        """
        source /etc/profile.d/modules.sh
        module load gcc/11.3.0
        module load samtools/1.17
        samtools view -b -F 0x800 -F 0x100 -@ 30 {input.bam} > {output.outputbam}
        samtools index -@ 30 {output.outputbam}
        """

rule cigarProcessing:
    input:
        script = "code/cigar_processing_region.py",
        bam = files['outputbam'],
        bai = files['outputbai'], 
        pri = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt",
        alt = "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt",
        bed1 = files['geneOut_pri'],  
        bed2 = files['geneOut_alt'],
        bed3 = files['geneOut_alt2']
    output:
        files['final_IGH'],
        files['final_IGK'],
        files['final_IGL'],
        files['final_nonIG']
    params:
        species = "{species}",
        corrected = CORRECTED[0],
        genes=["IGH", "IGK", "IGL"],
        bam = files['outputbam'],
        bai = files['outputbai']
    shell:
        """
        declare -A gene_loci
        if [ {params.corrected} == "True" ]; then
            mkdir -p errorStats/{params.species}.lja
            for gene in {params.genes}; do
                most_freq_chr_pri=$(awk '{{print $1}}' gene_position/{params.species}/lja/{params.species}_${{gene}}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 1 | awk '{{print $2}}')
                most_freq_chr_sec=$(awk '{{print $1}}' gene_position/${params.species}/lja/${params.species}_${{gene}}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 2 | tail -n 1 | awk '{{print $2}}')
                #most_freq_chr_third=$(awk '{{print $1}}' gene_position/${params.species}/lja/${params.species}_${{gene}}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 3 | tail -n 1 | awk '{{print $2}}')
                priloc=$(awk -v chr="$most_freq_chr_pri" 'BEGIN{{min=""; max=0}} $1 == chr {{if (min=="" || $2<min) min=$2; if ($3>max) max=$3}} END {{print chr":"min"-"max}}' gene_position/{params.species}/lja/{params.species}_${{gene}}_pos.sorted.bed)
                altloc=$(awk -v chr="$most_freq_chr_sec" 'BEGIN{{min=""; max=0}} $1 == chr {{if (min=="" || $2<min) min=$2; if ($3>max) max=$3}} END {{print chr":"min"-"max}}' gene_position/{params.species}/lja/{params.species}_${{gene}}_pos.sorted.bed)
                gene_loci[${{gene}}_pri]=$priloc
                gene_loci[${{gene}}_alt]=$altloc
            done
            python {input.script} {params.bam} ${{gene_loci[IGH_pri]}} ${{gene_loci[IGH_alt]}} ${{gene_loci[IGK_pri]}} ${{gene_loci[IGK_alt]}} ${{gene_loci[IGL_pri]}} ${{gene_loci[IGL_alt]}} {params.species}.lja
        else
            mkdir -p errorStats/{params.species}
            for gene in {params.genes}; do
                most_freq_chr_pri=$(awk '{{print $1}}' gene_position/{params.species}/pri/{params.species}_${{gene}}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 1 | awk '{{print $2}}')
                most_freq_chr_alt=$(awk '{{print $1}}' gene_position/{params.species}/alt/{params.species}_${{gene}}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 1 | awk '{{print $2}}')
                priloc=$(awk -v chr="$most_freq_chr_pri" 'BEGIN{{min=""; max=0}} $1 == chr {{if (min=="" || $2<min) min=$2; if ($3>max) max=$3}} END {{print chr":"min"-"max}}' gene_position/{params.species}/pri/{params.species}_${{gene}}_pos.sorted.bed)
                altloc=$(awk -v chr="$most_freq_chr_alt" 'BEGIN{{min=""; max=0}} $1 == chr {{if (min=="" || $2<min) min=$2; if ($3>max) max=$3}} END {{print chr":"min"-"max}}' gene_position/{params.species}/alt/{params.species}_${{gene}}_pos.sorted.bed)
                gene_loci[${{gene}}_pri]=$priloc
                gene_loci[${{gene}}_alt]=$altloc
            done
            echo ${{gene_loci[IGH_pri]}}
            python {input.script} {params.bam} ${{gene_loci[IGH_pri]}} ${{gene_loci[IGH_alt]}} ${{gene_loci[IGK_pri]}} ${{gene_loci[IGK_alt]}} ${{gene_loci[IGL_pri]}} ${{gene_loci[IGL_alt]}} {params.species}
        fi
        """

rule coverageAnalysis:
    input:
        "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.pri.igdetective/combined_genes_IGL.txt",
        "/home1/zhuyixin/sc1/AssmQuality/igGene/{species}.alt.igdetective/combined_genes_IGL.txt",
        bam = files['outputbam'],
        bai = files['outputbai'],
        bed1 = files['geneOut_pri'],  
        bed2 = files['geneOut_alt'],
        bed3 = files['geneOut_alt2']
    output:
        files['pileup_IGH'],
        files['pileup_IGK'],
        files['pileup_IGL']
    params:
        species = "{species}",
        assemblies = files['genome'],
        script = "code/coverage_snake.sh",
        genes=["IGH", "IGK", "IGL"],
        bam = files['outputbam'],
        bai = files['outputbai'],
    shell:
        """
        sbatch --partition=qcb {params.script} -s {params.species} -a {params.assemblies} -b /home1/zhuyixin/sc1/AssmQuality/{params.bam}
        """
#####TODOs

#rule stats analysis:
#rule qv0: