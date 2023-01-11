cd ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/
# 3 day time max
for each in ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/*_1_qc.fastq;
do export FNAME="${each##*/}";export BASE="${FNAME%_1_qc.fastq*}"; qsub -N "Sample_${BASE}" ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads_mapping_stats_qc.sh ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/${FNAME}; done
#--
# map subsetted reads to linear references
for each in ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/vgmap_comparisons/*_1_qc_subset.fastq;
for each in ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/vgmap_comparisons/SRR7471499_1_qc_subset.fastq;
do export FNAME="${each##*/}";export BASE="${FNAME%_R1.fastq*}"; qsub -N "Sample_${BASE}" ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads_mapping_stats_qc_subset.sh ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/${FNAME}; done


# move everything into appropriate directory - bwa outputs to stdout current dir (where fastq files were)
DATAPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/"
cp $DATAPATH*.fastq ./


##################################
#!/bin/bash
#$ -l mem_free=20G
#$ -l h_vmem=21G
#$ -l h_fsize=20G
#$ -m e
#$ -M dduchen3@jhu.edu
#SBATCH --time=72:0:0
#SBATCH --partition=shared

# Map fastq files to all 44 HBV references using BWA mem, calculate stats, and generate QC report.
#TESTING:
#BASE=SRR7471499_1_qc.fastq
#
OUTPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/qc_reads/linear_ref/"
REFSPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/reference/all44_refs/individual_refs/"
DATAPATH=$OUTPATH
#
cd $OUTPATH
#
BASE=$1
#
echo file num_reads mapped_reads secondary_aln mapped_percent > ${BASE%*_1_qc.fastq}"_map_summary.txt"
#
ml bwa
ml samtools
ml java
#
for each in $REFSPATH*.fasta;
do export REF="${each##*/}"; echo $REFSPATH${REF};
export specificRef=${REF##*Sequence_}; export specificRef=${specificRef%*.fasta};
cd $OUTPATH;
#export fullBASE="$DATAPATH${BASE}"; echo $fullBASE; bwa mem -t 16 -M $REFSPATH${REF} $fullBASE | samtools view -Sb - | samtools sort - > ${BASE%*_1_qc.fastq}"_raw_qc_"${specificRef}.bam;
export fullBASE="${BASE}"; echo $fullBASE;cd $OUTPATH; bwa mem -t 16 -M $REFSPATH${REF} $fullBASE ${fullBASE%1_qc.fastq}"2_qc.fastq" | samtools view -Sb - | samtools sort - > ${BASE%*_1_qc.fastq}"_raw_qc_"${specificRef}.bam;
# stats
samtools flagstat ${BASE%*_1_qc.fastq}"_raw_qc_${specificRef}.bam" > ${BASE%*_1_qc.fastq}"_raw_qc_"${specificRef}.stats;
done
# - Just combine into table after everything:
#echo file num_reads mapped_reads secondary_aln mapped_percent > C4_raw_qc_map_summary.txt
#for each in *.stats; do sh ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/map_summary_stats.sh $each >> C4_raw_qc_map_summary.txt;
#done
# updated script: ~/work/dduchen3/projects/HBV/reference/linear_mapping_stats.sh
#--
# ---Now: run once upon completion of batch script
#-- qualimap reports
#ml java
#for each in *.bam;do export BAM="${each##*/}"; echo $BAM;~/code/dduchen_tools/qualimap_v2.2.1/qualimap bamqc -bam $BAM -gd HUMAN -ip -outfile ${BAM%.bam}_${specificRef}"_qualimap_report";done
# -- multiqc report
#multiqc .
