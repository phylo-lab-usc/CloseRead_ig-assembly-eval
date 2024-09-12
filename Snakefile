SPECIES = ["fScoMax1"]
fastqdir = ["hifi_fastq"]
HAPLOID = ["False"]
knownLoci = config.get("knownLoci", False)
known_lociPos = config.get("loci_file", "NA")  # Default value if not provided

HOME = "/home1/zhuyixin/zhuyixin_proj/AssmQuality"
files = dict()
files["merged.bam"] = "{HOME}/aligned_bam/{species}/{species}_merged_sorted.bam"
files["merged.csi"] = "{HOME}/aligned_bam/{species}/{species}_merged_sorted.bam.csi"
files['merged.genome'] = "{HOME}/assemblies/{species}.merged.fasta"
files['pri.genome'] = "{HOME}/assemblies/{species}.pri.fasta"
files['alt.genome'] = "{HOME}/assemblies/{species}.alt.fasta"
files['priRead.bam'] = "{HOME}/aligned_bam/{species}/{species}_merged_sorted_primary.bam"
files['priRead.csi'] = "{HOME}/aligned_bam/{species}/{species}_merged_sorted_primary.bam.csi"
files['final.genePos_IG'] = "{HOME}/gene_position/{species}.final.Ig_loci.txt"
files['read_IGH.out'] = "{HOME}/errorStats/{species}/IGH.txt"
files['read_IGK.out'] = "{HOME}/errorStats/{species}/IGK.txt"
files['read_IGL.out'] = "{HOME}/errorStats/{species}/IGL.txt"
files['read_nonIG.out'] = "{HOME}/errorStats/{species}/nonIG.txt"
files['pileup_IGH.out'] = "{HOME}/errorStats/{species}/IGH_pri_pileup.txt"
files['pileup_IGK.out'] = "{HOME}/errorStats/{species}/IGK_pri_pileup.txt"
files['pileup_IGL.out'] = "{HOME}/errorStats/{species}/IGL_pri_pileup.txt"
files['cigarend'] = "{HOME}/errorStats/{species}/cigar.end"
files['pileupend'] = "{HOME}/errorStats/{species}/pileup.end"
files['igDetect.pri'] = "{HOME}/igGene/{species}.pri.txt"
files['igDetect.alt'] = "{HOME}/igGene/{species}.alt.txt"

rule all:
    input:
        expand(files['cigarend'], species = SPECIES, haploid = HAPLOID, HOME = HOME),
        expand(files['pileupend'], species = SPECIES, haploid = HAPLOID, HOME = HOME)

rule dataPrepAutomate:
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        files["merged.bam"],
        files["merged.csi"]
    params:
        species = "{species}",
        source = fastqdir[0],
        haploid = HAPLOID[0]
    shell:
        """
        sbatch --partition=gpu {input.script} -s {params.species} -w {params.source} -h {params.haploid} -d {HOME}
        """

rule convertPrimaryBam:
    input:
        bam = files["merged.bam"],
    output:
        outputbam = files['priRead.bam'],
        outputcsi = files['priRead.csi']
    shell:
        """
        source /etc/profile.d/modules.sh
        module load gcc/11.3.0
        module load samtools/1.17
        samtools view -b -F 0x800 -F 0x100 -@ 30 {input.bam} > {output.outputbam}
        samtools index -c -@ 30 {output.outputbam}
        """

rule lociLocation:
    input:
        pri_genome = files['pri.genome'],
        alt_genome = files['alt.genome'] if (HAPLOID[0]=='False') else [],
        lociScript = "code/igDetective.sh"
    output:
        files['igDetect.pri'],
        files['igDetect.alt']
    params:
        pri_outdir = "{HOME}/igGene/{species}.pri.igdetective/",
        alt_outdir = "{HOME}/igGene/{species}.alt.igdetective/",
        species = "{species}"
    shell:
        """
        if [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/{params.species}.pri.txt" ]; then
            sbatch --partition=gpu {input.lociScript} {input.pri_genome} {params.pri_outdir} {params.species} pri {HOME}
        fi
        if [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/{params.species}.alt.txt" ]; then
            sbatch --partition=gpu {input.lociScript} {input.alt_genome} {params.alt_outdir} {params.species} alt {HOME}
        fi
        """

#process igDetective result into a loci file, note this only reflect the most frequency chromosome. Need further modification to reflect distal loci.
rule finalIGLoci:
    input:
        script = "code/finalGene.py",
        pri = files['igDetect.pri'],
        alt = files['igDetect.alt']
    output:
        finalout = files['final.genePos_IG']
    params:
        species = "{species}"
    shell:
        """
        rm -rf {output.finalout}
        touch {output.finalout}
        python {input.script} {params.species} {HOME}
        """

rule cigarProcessing:
    input:
        script = "code/cigar_processing_region.py",
        bam = files['priRead.bam'],
        csi = files['priRead.csi'],
        finalout = files['final.genePos_IG'] if not knownLoci else known_lociPos
    output:
        files['cigarend']
    params:
        IGH_out = files['read_IGH.out'],
        IGK_out = files['read_IGK.out'],
        IGL_out = files['read_IGL.out'],
        nonIG_out = files['read_nonIG.out'],
        species = "{species}"
    shell:
        """
        mkdir -p {HOME}/errorStats/{params.species}
        rm -rf {params.IGH_out}
        rm -rf {params.IGK_out}
        rm -rf {params.IGL_out}
        rm -rf {params.nonIG_out}
        python {input.script} {input.bam} {input.finalout} {params.species}
        touch {HOME}/errorStats/{params.species}/cigar.end
        """

rule coverageAnalysis:
    input:
        finalout = files['final.genePos_IG'] if not knownLoci else known_lociPos,
        bam = files['priRead.bam'],
        csi = files['priRead.csi'],
        script = "code/coverage_snake.sh"
    output:
        files['pileupend']
    params:
        species = "{species}",
        assemblies = files['merged.genome'],
        IGH_out = files['pileup_IGH.out'],
        IGK_out = files['pileup_IGK.out'],
        IGL_out = files['pileup_IGL.out'],
    shell:
        """
        mkdir -p {HOME}/errorStats/{params.species}
        rm -rf {params.IGH_out}
        rm -rf {params.IGK_out}
        rm -rf {params.IGL_out}
        rm -rf {HOME}/errorStats/{params.species}/*_pileup.txt
        sbatch --partition=gpu {input.script} -s {params.species} -a {params.assemblies} -b {input.bam} -f {input.finalout} -d {HOME}
        """
