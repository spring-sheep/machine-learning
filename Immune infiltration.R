rm(list = ls())
library(estimate)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")  
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='LIHC'
scores=estimate(dat,pro)
head(scores)
library(tibble)
scores = rownames_to_column(as.data.frame(scores),"id")
scores$id = str_replace_all(scores$id,"\\.","-")
#save(scores,file = 'Rdata/scores.Rdata')

#免疫浸润
#浸润计算
geneset = rio::import("zip/m28.xlsx",skip = 2)
geneset = split(geneset$Metagene,geneset$`Cell type`)
lapply(geneset[1:3], head)
library(GSVA)
re <- gsva(as.matrix(dat), geneset, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)
#save(re,file = 'Rdata/re.Rdata')

#共识聚类分析
library('ConsensusClusterPlus')
#ConsensusClusterPlus
results = ConsensusClusterPlus(as.matrix(re),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='pictures/ConsensusCluster',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot=NULL)
#write.csv(re,file = 'fig.csv')

icl <- calcICL(results,title = 'pictures/ConsensusCluster',plot = 'pdf')

## PAC = Proportion of ambiguous clustering 模糊聚类比例
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK

PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw(base_rect_size = 1.5)+
  geom_point(size=4,shape=21,color='darkred',fill='orange')+
  ggtitle('Proportion of ambiguous clustering')+
  xlab('Cluster number K')+ylab(NULL)+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=13))
ggsave(filename = 'pictures/ConsensusCluster/PAC.pdf',width = 3.8,height = 4)