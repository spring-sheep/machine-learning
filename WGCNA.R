rm(list = ls())
library(WGCNA)
library(ggplot2)
collectGarbage()
load("D:/aTCGA机器学习mRNA/Rdata/tumor样本与分组信息（用于wgcna）.Rdata")
load("D:/aTCGA机器学习mRNA/Rdata/免疫分组.Rdata")

match(names(Cluster),rownames(Group))
Group$Cluster <- Cluster
Group$Cluster <-ifelse(Group$Cluster == 1,0,1)
table(Group$Cluster)
Group <- Group[,c(9,1:8)]

datExpr <- exp_tumor

#2.2.1选择软阈值功率: 网络拓扑分析
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 5,
                        RsquaredCut = 0.8,
                        networkType = "unsigned")

FitIndices <- sft[["fitIndices"]]
write.csv(FitIndices,file = 'submit data/fig2c.csv')

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#结果如图所示。我们选择功率14，它是无标度拓扑匹配指数达到0.90的最低功率。
sft[["powerEstimate"]]#:查看最佳功率

#2.2.2共表达相似性和邻接性
#现在我们使用软阈值功率6计算邻接:
softPower = 14 ;
adjacency = adjacency(datExpr, power = softPower,type = "unsigned")
collectGarbage()

#2.2.3拓扑重叠矩阵
#为了最小化噪声和虚假关联的影响，我们将邻接转换成拓扑重叠矩阵，并计算相应的不相似度
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = "unsigned");
#save(TOM,file = 'Rdata/WGCNA-TOM.Rdata')
load("D:/aTCGA机器学习mRNA/Rdata/WGCNA-TOM.Rdata")

rm(adjacency)
collectGarbage()
dissTOM = 1-TOM
rm(TOM)
collectGarbage()
#2.2.4使用 TOM 进行聚类
#我们现在使用分层聚类来产生基因的分层聚类树(树状图)。注意，我们使用的函数 hclust 提供了比标准 hclust 函数更快的分层集群例程
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
#最后一个命令绘制的聚类树状图如图所示。在聚类树(树状图)中，每一片叶子，也就是一条短的垂直线，对应着一个基因。
#树状图组的分支密集相互连接，高度共表达的基因。模块识别相当于单个分支的识别(“从树状图中剪下分支”)。
#有几种方法用于树枝切割，我们的标准方法是从包装 DynamicTreeCut 动态树切割。下一段代码说明了它的用法
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, #该值越高（或者如果为TRUE），将产生更多更小的簇
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            #method = "tree"   结果不好，使用默认数字
);
table(dynamicMods)
#该函数返回标记为1-22的22个模块，从最大到最小。标签0是为未分配的基因保留的。上面的命令列出了模块的大小。
#我们现在在基因树状图下绘制模块分配
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#2.2.5合并表达式配置文件非常相似的模块
#动态树切割可以识别表达式配置文件非常相似的模块。合并这些模块可能是谨慎的，因为它们的基因是高度共表达的。
#为了量化整个模块的共表达相似性，我们计算了它们的特征基因，并将它们按照相关性进行聚类
# Calculate eigengenes

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes


# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#我们选择0.25的高度切割，对应于0.75的相关性，以进行合并
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#为了了解合并对我们的模块颜色造成了什么影响，我们再次绘制基因树状图，原始和合并的模块颜色位于下面
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

#在随后的分析中，我们将在 mergedColors 中使用合并的模块颜色。我们保存相关的变量，以便在本教程的后续部分中使用
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "Rdata/wgcna模块筛选输入数据.RData")
#load(file = "Rdata/wgcna模块筛选输入数据.RData")


#挑选模块
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

clinical <- Group[,c(-4,-5)]
clinical$Cluster <- ifelse(clinical$Cluster==1,2,1)
moduleTraitCor = cor(MEs, clinical, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#由于我们有相当多的模块和特性，一个合适的图形表示将有助于阅读表格。我们用相关值对每个关联进行颜色编码
sizeGrWindow(15,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))