cd ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/graph44_ref

# move everything into appropriate directory - bwa outputs to stdout current dir (where fastq files were)
DATAPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/"
cp $DATAPATH*_pear.assembled.150bp.fastq ./

# 3 day time max
for each in /home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/*1_qc.fastq;
do export FNAME="${each##*/}";export BASE="${FNAME%_1_qc.fastq}"; qsub -N "${BASE}" ~/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/alignment/C4/graph_mapping_and_stats.sh /home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/${FNAME}; done
#--
# .sh script -->

### testing new gbz haplotype-informed threads
#vg giraffe -Z $graphdir'HBV_44refs.smooth_norm.giraffe.gbz' -H $graphdir'HBV_44refs.smooth_norm.augmented.gbwt' -m $graphdir'HBV_44refs.smooth_norm.min' -d $graphdir'HBV_44refs.smooth_norm.dist' -f ${DATAPATH}SRR7471501_1_qc.fastq -f ${DATAPATH}SRR7471501_2_qc.fastq -o gam -p --paired-distance-limit 3 > SRR7471501_smooth_giraffe_new.gam

##################################
#!/bin/bash
#$ -l mem_free=75G
#$ -l h_vmem=75G
#$ -l h_fsize=75G
#$ -m e
#$ -M dduchen3@jhu.edu
#SBATCH --time=168:0:0
#SBATCH --partition=unlimited

# Map fastq files to all 44 HBV references using BWA mem, calculate stats, and generate QC report.
#TESTING:
#BASE=SRR7471499_1_qc.fastq
#
OUTPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/graph44_ref/"
graphdir="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/reference/pggb_ref/HBV_44Refs_pggb/"
DATAPATH="/home-1/dduchen3@jhu.edu/work/dduchen3/projects/HBV/Longitudinal_CHB_SRP152039/NonSeroconverters/C4/C4_qc/"
#
cd $OUTPATH
#
#BASE=$1
#
# giraffe -- completed
#vg giraffe -Z $graphdir'HBV_44refs.smooth_norm.giraffe.gbz' -H $graphdir'HBV_44refs.smooth_norm.augmented.gbwt' -m $graphdir'HBV_44refs.smooth_norm.min' -d $graphdir'HBV_44refs.smooth_norm.dist' -f $DATAPATH${FNAME} -f $DATAPATH${FNAME%_1_qc.fastq}"_2_qc.fastq" -o gam -p > ${FNAME%_1_qc.fastq}"_smooth_giraffe.gam"
vg giraffe -Z $graphdir'HBV_44refs.seqwish_norm.giraffe.gbz' -H $graphdir'HBV_44refs.seqwish_norm.augmented.gbwt' -m $graphdir'HBV_44refs.seqwish_norm.min' -d $graphdir'HBV_44refs.seqwish_norm.dist' -f $DATAPATH${FNAME} -f $DATAPATH${FNAME%_1_qc.fastq}"_2_qc.fastq" -o gam -p > ${FNAME%_1_qc.fastq}"_seqwish_giraffe.gam"
#vg stats -a $OUTPATH${FNAME%_1_qc.fastq}"_smooth_giraffe.gam" > $OUTPATH${FNAME%_1_qc.fastq}"_smooth_giraffe.stats"
vg stats -a $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_giraffe.gam"  > $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_giraffe.stats"
#
# more exact/accurate alignments -->
# vg map
#vg map -x $graphdir'HBV_44refs.smooth_norm.xg' -g $graphdir'HBV_44refs.smooth_norm.gcsa' -f $DATAPATH${FNAME} -f $DATAPATH${FNAME%_1_qc.fastq}"_2_qc.fastq" --log-time > $OUTPATH${FNAME%_1_qc.fastq}"_smooth_map.gam"
#vg map -x $graphdir'HBV_44refs.seqwish_norm.xg' -g $graphdir'HBV_44refs.seqwish_norm.gcsa' -f $DATAPATH${FNAME} -f $DATAPATH${FNAME%_1_qc.fastq}"_2_qc.fastq" --log-time > $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_map.gam"
#vg stats -a $OUTPATH${FNAME%_1_qc.fastq}"_smooth_map.gam" > $OUTPATH${FNAME%_1_qc.fastq}"_smooth_map.stats"
#vg stats -a $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_map.gam" > $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_map.stats"
# vg mpmap --> much.... much slower
#vg mpmap -x $graphdir'HBV_44refs.smooth_norm.xg' -g $graphdir'HBV_44refs.smooth_norm.gcsa' -d $graphdir'HBV_44refs.smooth_norm.dist' -n DNA -f $DATAPATH${FNAME} -f $DATAPATH${FNAME%_1_qc.fastq}"_2_qc.fastq" -F GAM > $OUTPATH${FNAME%_1_qc.fastq}"_smooth_mpmap.gam"
#vg mpmap -x $graphdir'HBV_44refs.seqwish_norm.xg' -g $graphdir'HBV_44refs.seqwish_norm.gcsa' -d $graphdir'HBV_44refs.seqwish_norm.dist' -n DNA -f $DATAPATH${FNAME} -f $DATAPATH${FNAME%_1_qc.fastq}"_2_qc.fastq" -F GAM > $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_mpmap.gam"
#vg stats -a $OUTPATH${FNAME%_1_qc.fastq}"_smooth_mpmap.gam" > $OUTPATH${FNAME%_1_qc.fastq}"_smooth_mpmap.stats"
#vg stats -a $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_mpmap.gam" > $OUTPATH${FNAME%_1_qc.fastq}"_seqwish_mpmap.stats"
#
#
done
