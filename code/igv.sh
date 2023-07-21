function igv {
    mkdir /home1/zhuyixin/sc1/ImmAssm/igv
    cat /home1/zhuyixin/sc1/ImmAssm/gene_position/mCanLor1/gene_IGH_pos_sorted.bed | while read chrom start end sv score strand
    do
    echo $sv
    mkdir /home1/zhuyixin/sc1/ImmAssm/igv/snapshot
    bam="/home1/zhuyixin/sc1/ImmAssm/aligned_bam/mCanLor1/mCanLor1_merged.bam"
    echo "new"
    echo "preference SAM.SHOW_ALL_BASES 0"
    echo "genome /home1/zhuyixin/sc1/ImmAssm/assemblies/mCanLor1.pri.cur.20210315.fasta"
    echo "snapshotDirectory /home1/zhuyixin/sc1/ImmAssm/igv/snapshot"
    echo "load ${bam}"
    echo "load /home1/zhuyixin/sc1/ImmAssm/gene_position/mCanLor1/gene_IGH_pos_sorted.bed"
    echo "goto ${chrom}:${start}-${end}"
    echo "snapshot mCanLor1_merged_${sv}.png"
    done    > /home1/zhuyixin/sc1/ImmAssm/igv/igv.txt
}

igv
