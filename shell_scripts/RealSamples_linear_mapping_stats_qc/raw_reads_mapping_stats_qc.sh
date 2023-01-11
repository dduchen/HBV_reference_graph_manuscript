cd ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/raw_reads/linear_ref/
# 3 day time max
for each in ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/raw_reads/linear_ref/*_1.fastq;
do export FNAME="${each##*/}";export BASE="${FNAME%_1.fastq*}"; qsub -N "Sample_raw_${BASE}" ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/raw_reads_mapping_stats_qc.sh ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/raw_reads/linear_ref/${FNAME}; done
#--
# move everything into appropriate directory - bwa outputs to stdout current dir (where fastq files were)
DATAPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/"
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
OUTPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/raw_reads/linear_ref/"
REFSPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/reference/all44_refs/individual_refs/"
DATAPATH=$OUTPATH
#
cd $OUTPATH
#
BASE=$1
#
echo file num_reads mapped_reads secondary_aln mapped_percent > ${BASE%*_1.fastq}"_raw_map_summary.txt"
#
ml bwa
ml samtools
ml java
#
for each in $REFSPATH*.fasta;
do export REF="${each##*/}"; echo $REFSPATH${REF};
export specificRef=${REF##*Sequence_}; export specificRef=${specificRef%*.fasta};
cd $OUTPATH;
export fullBASE="${BASE}"; echo $fullBASE;cd $OUTPATH; bwa mem -t 16 -M $REFSPATH${REF} $fullBASE ${fullBASE%1.fastq}"2.fastq" | samtools view -Sb - | samtools sort - > ${BASE%*_1.fastq}"_raw_"${specificRef}.bam;
# stats
samtools flagstat ${BASE%*_1.fastq}"_raw_${specificRef}.bam" > ${BASE%*_1.fastq}"_raw_"${specificRef}.stats;
done
# - combine into table:
for each in *.stats; do sh ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/map_summary_stats.sh $each >> ${BASE%*_1.fastq}"_raw_map_summary.txt";
done
#--
# ---Now: run once upon completion of batch script
#-- qualimap reports
#ml java
#for each in *.bam;do export BAM="${each##*/}"; echo $BAM;~/code/dduchen_tools/qualimap_v2.2.1/qualimap bamqc -bam $BAM -gd HUMAN -ip -outfile ${BAM%.bam}_${specificRef}"_qualimap_report" -outformat PDF;done
# -- multiqc report
