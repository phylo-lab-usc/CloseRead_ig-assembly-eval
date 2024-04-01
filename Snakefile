SPECIES = ["mUrsArc2", "mPerMan1", "mThoBot1"]
SOURCE = ["hifi_fastq"]
HAPLOID = ["True"]

files = dict()
files["merged.bam"] = "aligned_bam/{species}/{species}_merged_sorted.bam"
files["merged.csi"] = "aligned_bam/{species}/{species}_merged_sorted.bam.csi"
files['merged.genome'] = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/{species}.merged.fasta"
files['pri.genome'] = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/{species}.pri.fasta"
files['alt.genome'] = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/{species}.alt.fasta"
files['priRead.bam'] = "aligned_bam/{species}/{species}_merged_sorted_primary.bam"
files['priRead.csi'] = "aligned_bam/{species}/{species}_merged_sorted_primary.bam.csi"
files['ref.genePos_IG'] = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/gene_position/{species}/ref_loci_details.txt"
files['Ig.genePos_IG'] = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/gene_position/{species}/Ig_loci_details.txt"
files['final.genePos_IG'] = "gene_position/{species}/final.Ig_loci.txt"
files['read_IGH.out'] = "errorStats/{species}/IGH.txt"
files['read_IGK.out'] = "errorStats/{species}/IGK.txt"
files['read_IGL.out'] = "errorStats/{species}/IGL.txt"
files['read_nonIG.out'] = "errorStats/{species}/nonIG.txt"
files['pileup_IGH.out'] = "errorStats/{species}/IGH_alt_pileup.txt"
files['pileup_IGK.out'] = "errorStats/{species}/IGK_alt_pileup.txt"
files['pileup_IGL.out'] = "errorStats/{species}/IGL_alt_pileup.txt"

rule all:
    input:
        expand(files['read_IGH.out'], species = SPECIES, haploid = HAPLOID),
        expand(files['read_IGK.out'], species = SPECIES, haploid = HAPLOID),
        expand(files['read_IGL.out'], species = SPECIES, haploid = HAPLOID),
        expand(files['read_nonIG.out'], species = SPECIES, haploid = HAPLOID),
        expand(files['pileup_IGH.out'], species = SPECIES, haploid = HAPLOID),
        expand(files['pileup_IGK.out'], species = SPECIES, haploid = HAPLOID),
        expand(files['pileup_IGL.out'], species = SPECIES, haploid = HAPLOID)


rule dataPrepAutomate:
    input:
        script = "code/dataPrepAutomated.sh"
    output:
        files["merged.bam"],
        files["merged.csi"]
    params:
        species = "{species}",
        source = SOURCE[0],
        haploid = HAPLOID[0]
    shell:
        """
        sbatch --partition=gpu {input.script} -s {params.species} -w {params.source} -h {params.haploid}
        """

rule convertPrimaryBam:
    input:
        bam = files["merged.bam"]
    output:
        outputbam = files['priRead.bam'],
        outputcsi = files['priRead.csi']
    shell:
        """
        source /etc/profile.d/modules.sh
        module load gcc/11.3.0
        module load samtools/1.17
        samtools view -b -F 0x800 -F 0x100 -@ 40 {input.bam} > {output.outputbam}
        samtools index -c -@ 40 {output.outputbam}
        """

rule lociLocation:
    input:
        pri_genome = files['pri.genome'],
        alt_genome = files['alt.genome'] if (HAPLOID[0]=='False') else [],
        lociScript = "code/refGene.sh"
    output:
        refout = files['ref.genePos_IG'],
        Igout = files['Ig.genePos_IG'],
        finalout = files['final.genePos_IG']
    params:
        pri_outdir = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/{species}.pri.igdetective/",
        alt_outdir = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/{species}.alt.igdetective/",
        haploid = HAPLOID[0],
        species = "{species}"
    shell:
        """
        read primary_result alternate_result <<< $({input.lociScript} -s {params.species} -h {params.haploid} -p {params.pri_outdir} -a {params.alt_outdir})
        if [ $alternate_result == "false" ] && [ {params.haploid} == "False" ]
        then
            sbatch --partition=qcb code/igDetective.sh {input.alt_genome} {params.alt_outdir} {params.species} alt
        elif [ $primary_result == "false" ]
        then
            sbatch --partition=qcb code/igDetective.sh {input.pri_genome} {params.pri_outdir} {params.species} pri
        else
            echo "exist, skip IGdetective"
        fi
        code/finalGene.sh -s {params.species} -h {params.haploid} -p $primary_result -a $alternate_result
        """


rule cigarProcessing:
    input:
        script = "code/cigar_processing_region.py",
        bam = files['priRead.bam'],
        csi = files['priRead.csi'], 
        finalout = files['final.genePos_IG']
    output:
        IGH_out = files['read_IGH.out'],
        IGK_out = files['read_IGK.out'],
        IGL_out = files['read_IGL.out'],
        nonIG_out = files['read_nonIG.out']
    params:
        species = "{species}",
    shell:
        """
        mkdir -p errorStats/{params.species}
        rm -rf {output.IGH_out}
        rm -rf {output.IGK_out}
        rm -rf {output.IGL_out}
        rm -rf {output.nonIG_out}
        python {input.script} {input.bam} {input.finalout} {params.species}
        """

rule coverageAnalysis:
    input:
        finalout = files['final.genePos_IG'],
        bam = files['priRead.bam'],
        csi = files['priRead.csi'],
    output:
        IGH_out = files['pileup_IGH.out'],
        IGK_out = files['pileup_IGK.out'],
        IGL_out = files['pileup_IGL.out']
    params:
        species = "{species}",
        assemblies = files['merged.genome'],
        script = "code/coverage_snake.sh"
    shell:
        """
        rm -rf {output.IGH_out}
        rm -rf {output.IGK_out}
        rm -rf {output.IGL_out}
        sbatch --partition=qcbr {params.script} -s {params.species} -a {params.assemblies} -b {input.bam} -f {input.finalout}
        """
