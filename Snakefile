configfile: "config.yaml"

SPECIES = config['speciesList']
fastqdir = config['fastqdir']
HAPLOID = config['haploid']
HOME = config['home']

knownLoci = config.get("knownLoci", False)
known_lociDir= config.get("loci_dir", "NA")  

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
files['logMinimap'] = "{HOME}/log/dataPrepAutomated.{species}.log"
files['logindex'] = "{HOME}/log/index.{species}.log"
files['logigDetective'] = "{HOME}/log/igDetective.{species}.log"
files['logcigarProcessing'] = "{HOME}/log/cigarProcessing.{species}.log"
files['logcoverageAnalysis'] = "{HOME}/log/coverageAnalysis.{species}.log"

if not knownLoci:
    igAnnotation = files['final.genePos_IG']
else:
    igAnnotation = lambda wildcards: f"{known_lociDir}/{wildcards.species}.customIG.txt"


rule all:
    input:
        expand(files['cigarend'], species = SPECIES, haploid = HAPLOID, HOME = HOME),
        expand(files['pileupend'], species = SPECIES, haploid = HAPLOID, HOME = HOME)


rule dataPrepAutomate:
    resources:
       mem="60G",
    threads: 32
    log: files["logMinimap"]
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        files["merged.bam"],
        files["merged.csi"]
    params:
        species = "{species}",
        source = fastqdir[0],
        haploid = HAPLOID[0],
        conda = config["condaPath"] 
    shell:
        """
        {input.script} -s {params.species} -w {params.source} -h {params.haploid} -d {HOME} -c {params.conda}
        """

rule convertPrimaryBam:
    resources:
       mem="60G",
    threads: 30
    log: files["logindex"]
    input:
        bam = files["merged.bam"],
    output:
        outputbam = files['priRead.bam'],
        outputcsi = files['priRead.csi']
    shell:
        """
        samtools view -b -F 0x800 -F 0x100 -@ 30 {input.bam} > {output.outputbam}
        samtools index -c -@ 30 {output.outputbam}
        """

rule lociLocation:
    resources:
       mem="30G",
    threads: 10
    log: files['logigDetective']
    input:
        pri_genome = files['pri.genome'],
        alt_genome = files['alt.genome'] if (HAPLOID[0]=='False') else [],
        lociScript = "code/igDetective.sh"
    output:
        files['igDetect.pri'],
        out = files['igDetect.alt'] if (HAPLOID[0]=='False') else []
    params:
        pri_outdir = "{HOME}/igGene/{species}.pri.igdetective/",
        alt_outdir = "{HOME}/igGene/{species}.alt.igdetective/",
        species = "{species}",
        igdetective_home = config["igdetective_home"],
        conda = config["condaPath"],
        haploid = HAPLOID[0],
        condaEnv = config["condaEnvPath"]
    shell:
        """
        if [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/{params.species}.pri.txt" ]; then
            {input.lociScript} {input.pri_genome} {params.pri_outdir} {params.species} pri {HOME} {params.igdetective_home} {params.conda} {params.condaEnv}
        fi
        if [ "{params.haploid}" == "False" ] && [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/{params.species}.alt.txt" ]; then
            {input.lociScript} {input.alt_genome} {params.alt_outdir} {params.species} alt {HOME} {params.igdetective_home} {params.conda} {params.condaEnv}
        fi
        """

#process igDetective result into a loci file, note this only reflect the most frequency chromosome. Need further modification to reflect distal loci.
rule finalIGLoci:
    resources:
       mem="30G",
    threads: 10
    input:
        script = "code/finalGene.py",
        pri = files['igDetect.pri'],
        alt = files['igDetect.alt'] if (HAPLOID[0]=='False') else []
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
    resources:
       mem="30G",
    threads: 10
    log: files['logcigarProcessing']
    input:
        script = "code/cigar_processing_region.py",
        bam = files['priRead.bam'],
        csi = files['priRead.csi'],
        finalout = igAnnotation
    output:
        files['cigarend']
    params:
        IGH_out = files['read_IGH.out'],
        IGK_out = files['read_IGK.out'],
        IGL_out = files['read_IGL.out'],
        nonIG_out = files['read_nonIG.out'],
        species = "{species}",
        conda = config["condaPath"],
        condaEnv = config["condaEnvPath"]
    shell:
        """
        conda init
        source {params.conda}
        conda activate ig-assembly-eval
        which python
        mkdir -p {HOME}/errorStats/{params.species}
        rm -rf {params.IGH_out}
        rm -rf {params.IGK_out}
        rm -rf {params.IGL_out}
        rm -rf {params.nonIG_out}
        {params.condaEnv}/ig-assembly-eval/bin/python {input.script} {input.bam} {input.finalout} {params.species}
        touch {HOME}/errorStats/{params.species}/cigar.end
        """

rule coverageAnalysis:
    resources:
       mem="30G",
    threads: 10
    log: files['logcoverageAnalysis']
    input:
        finalout = igAnnotation,
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
        conda = config["condaPath"]
    shell:
        """
        mkdir -p {HOME}/errorStats/{params.species}
        rm -rf {params.IGH_out}
        rm -rf {params.IGK_out}
        rm -rf {params.IGL_out}
        rm -rf {HOME}/errorStats/{params.species}/*_pileup.txt
        {input.script} -s {params.species} -a {params.assemblies} -b {input.bam} -f {input.finalout} -d {HOME} -c {params.conda}
        """
