#
# JHPCE saved: ~/Data/HBV/reference/genotyping/BestReference44_PGGBedyeetSeqwish.R
# Updated version: ~/Data/HBV/reference/genotyping/BestReference44_PGGBwfmashSeqwish.R
#########################################################################################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
options(stringsAsFactors=F)
best_ref <- function(input_tsv,marcc){
# MARCC or JHPCE: :
if(marcc==TRUE){
  load("~/work/dduchen3/projects/HBV/genotyping/pangenome_ref/HBV_44Refs_pggb_wfmash_norm_Path_Specific_Nodes_wWeights.R")
  } else {
  load("~/Data/HBV/reference/genotyping/HBV_44Refs_pggb_wfmash_norm_Path_Specific_Nodes_wWeights.R")
  }
# load appropriate alignment files
dat<-fread(input=input_tsv) # more sequence mapped to graph, use mpmap
#
# base-level coverage - get median coverage/node
datsum<-ddply(dat, c("node.id"), summarise,
               bplen    = length(node.id),
               median = median(coverage),
               mean = mean(coverage),
               sd   = sd(coverage))
#-- add the calculated 'weight' for each node
node_weights<-node_weights[order(node_weights$node),]
datsum<-datsum[order(datsum$node.id),]
if(all(node_weights$node==datsum$node.id)){
#	print("Reference graph nodes match alignment nodes")
	datsum$node.weight<-node_weights$weight
	predict<-hbv_44ref_complete
	mapped<-as.character(datsum[datsum$median>0,]$node.id)
	map_summary<-data.frame()
	for(i in 1:length(predict)){
		ref<-names(predict[i])
		mapped_nodes<-mapped[which(mapped %in% as.character(predict[[i]]$nodes))]
		#print(paste0(ref,": ",sum(datsum[datsum$node.id %in% mapped_nodes,]$bplen),"bp mapped"))
		#print(paste0(ref,": ",sum(datsum[datsum$node.id %in% mapped_nodes,]$node.weight)," weighted score"))
		a<-ref
		b<-sum(datsum[datsum$node.id %in% mapped_nodes,]$bplen)
		c<-sum(datsum[datsum$node.id %in% mapped_nodes,]$node.weight)
		map_summary<- rbind(map_summary,c(a,b,c))
		}
	colnames(map_summary)<-c("Reference_Path","Mapped_BP","Weighted_Score")
	map_summary[,2]<-as.numeric(map_summary[,2]);map_summary[,3]<-as.numeric(map_summary[,3])
# save output, then forward best ref to stdout
  fwrite(map_summary,file=paste0(gsub("_coverage.tsv","",input_tsv),"_44ref_AlnSummary.txt"),sep="\t",col.names=T,row.names=F,quote=F)
  write.table(map_summary[map_summary$Weighted_Score==max(map_summary$Weighted_Score),]$Reference_Path,file="best_reference.txt",sep="\t",col.names=F,row.names=F,quote=F)
  return(cat(map_summary[map_summary$Weighted_Score==max(map_summary$Weighted_Score),]$Reference_Path))} else {
	return(cat("Reference graph nodes do not match alignment nodes, did you use graph/index from:
  HBV_44refs.seqwish_norm"))}
}
best_ref(input_tsv=args[1],marcc=args[2])
