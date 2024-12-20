rm(list = ls()); gc()
setwd("/home/chenxu_OD/perturb")
library(Seurat)
library(Matrix)
library(stringr)
data <- Read10X(data.di = '03_scCRISPR_report_single/THP_1_HSV_4/sample_filtered_feature_bc_matrix')
data1 <- Read10X(data.dir = '03_scCRISPR_report_single/THP-1_HSV_1/sample_filtered_feature_bc_matrix')
data2 <- Read10X(data.dir = '03_scCRISPR_report_single/THP_1_HSV_2/sample_filtered_feature_bc_matrix')
data3 <- Read10X(data.dir = '03_scCRISPR_report_single/THP_1_HSV_3/sample_filtered_feature_bc_matrix')
str_sub(colnames(data2$`Gene Expression`), -1, -1) <- "2"
str_sub(colnames(data3$`Gene Expression`), -1, -1) <- "3"
str_sub(colnames(data$`Gene Expression`), -1, -1) <- "4"
str_sub(colnames(data2$`CRISPR Guide Capture`), -1, -1) <- "2"
str_sub(colnames(data3$`CRISPR Guide Capture`), -1, -1) <- "3"
str_sub(colnames(data$`CRISPR Guide Capture`), -1, -1) <- "4"
data_all<-cbind(data1$`Gene Expression`,data2$`Gene Expression`,data3$`Gene Expression`,data$`Gene Expression`)

data_all_cri<-cbind(data1$`CRISPR Guide Capture`,data2$`CRISPR Guide Capture`,data3$`CRISPR Guide Capture`,data$`CRISPR Guide Capture`)

data_a<-list(data_all,data_all_cri)
object <- CreateSeuratObject(counts = data_all, project = "perturb3",min.cells=3,min.features=200)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
object <- subset(object, subset = nFeature_RNA > 200 & percent.mt < 20)
object<-NormalizeData(object,normalization.method = "LogNormalize",scale.factor = 10000)

object<-ScaleData(object,feature=rownames(object))

object<-FindVariableFeatures(object,selection.method="mean.var.plot")

object<-RunPCA(object, features = VariableFeatures(object = object))

object<-FindNeighbors(object,dim=1:15)

object<-FindClusters(object,resolution=0.5)

DimPlot(object,reduction="pca")
object <- RunUMAP(object, dims = 1:10)
DimPlot(object, reduction = "umap")
object <- RunTSNE(object, dims = 1:10)
FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "Infection_Level",pt.size = 0.2)
FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Infection_Level",pt.size = 0.2)
######################
hsv1<-read.table("hsv1.txt",sep = '\t')
hsv2<-read.table("hsv2.txt",sep = '\t')
hsv3<-read.table("hsv3.txt",sep = '\t')
hsv4<-read.table("hsv4.txt",sep = '\t')
hsv_a<-rbind(hsv1,hsv2,hsv3,hsv4)
hsv_a<-hsv_a[colnames(object),]
# 假设 'seurat_obj' 是你的Seurat对象，'hsv1' 是包含细胞感染数据的DataFrame
# 确保hsv1的行名或列名与Seurat对象中的细胞名相匹配
hsv_a$Infected <- ifelse(hsv_a$value > 0, "Infected", "Not Infected")
hsv_a$Infection_Level <- ifelse(hsv_a$value > 1158, "High",
                                ifelse(hsv_a$value > 0, "Low", "Uninfected"))

# 将感染状态添加到 Seurat 对象的元数据中
object <- AddMetaData(object, metadata = hsv_a$Infection_Level, col.name = "Infection_Level")
object <- AddMetaData(object, metadata = hsv_a$Infected, col.name = "Infected")

object[["virus_load"]] <- log1p(hsv_a$value)
library(ggplot2)
cols <- colorRampPalette(c("lightgrey", "yellow", "red"))(100)
# 使用95%或99%分位数作为颜色映射的上限
FeaturePlot(object, features = 'virus_load',reduction = "tsne") +
  scale_color_gradientn(colors = c('lightgrey', 'yellow', 'red'),
                        limits = c(0, 10),    # 设置颜色映射的数据范围
                        breaks = seq(0, 10, 2),  # 设置颜色条上的标记点
                        labels = seq(0, 10, 2))  # 设置颜色条标签
FeaturePlot(object, features = 'virus_load',reduction = "umap") +
  scale_color_gradientn(colors = c('lightgrey', 'yellow', 'red'),
                        limits = c(0, 10),    # 设置颜色映射的数据范围
                        breaks = seq(0, 10, 2),  # 设置颜色条上的标记点
                        labels = seq(0, 10, 2))  # 设置颜色条标签
############################
mat<-data_all_cri
mat<-as.matrix(mat)
mat[mat < 20] = 0
mat<-mat[,colnames(object)]
rownames(mat) <- c(substr(rownames(mat)[1:76441], 1, nchar(rownames(mat)) - 2),rep('Non-Targeting',1000))

# 将行名相同的行相加
reduced_mat <- aggregate(mat, by=list(rownames(mat)), FUN=sum)

# 将结果中的Group.1列（包含行名）设为新矩阵的行名
rownames(reduced_mat) <- reduced_mat$Group.1
reduced_mat<- reduced_mat[, -1]  # 移除第一列，因为它现在已经成为行名

# 打印最终矩阵
jun<-colSums(reduced_mat[c(1, 3, 5, 7), ])
nov<-colSums(reduced_mat[c(2, 4, 6, 8), ])
reduced_mat<-reduced_mat[-c(1:8),]
reduced_mat<-rbind(reduced_mat,jun,nov)
###############################
library(ggplot2)
library(patchwork)

# 每个基因mRNA数量
sumGrna<-apply(reduced_mat,1,sum)
df <- data.frame(gRNA = sumGrna)

# 绘制带有紫色、绿色、黄色渐变的直方图
p_combined <- ggplot(df, aes(x = gRNA)) +
  # 绘制直方图
  geom_histogram(
    aes(fill = after_stat(count)),  # 将填充颜色映射到计数值
    binwidth = 1000,  # 设置柱子宽度
    color = "black",
    alpha = 0.6  # 设置透明度，方便看到箱线图
  ) +
  scale_fill_gradientn(
    colors = c("#FFD700", "green", "#6A0DAD")  # 设置颜色渐变
  ) +
  # 绘制箱线图
  geom_boxplot(
    aes(y = 1950),  # 将箱线图的位置放置在 y = 0
    width = 150,  # 设置箱线图的宽度
    fill = "#FFD700",  # 设置箱线图的填充颜色
    color = "darkblue",  # 设置箱线图的边框颜色
    notch = TRUE  # 使用缺口
  ) +
  # 添加标题和标签
  labs(
    title = "Distribution of each gRNA Numbers",
    x = "gRNA mRNA",
    y = "gRNA Number",
    fill = "gRNA mRNA"
  ) +
  # 设置x轴显示范围和刻度
  scale_x_continuous(
    limits = c(0, 70000),  # 设置x轴显示范围
    breaks = seq(0, 70000, by = 10000)  # 设置x轴刻度每10000一个刻度
  ) +
  # 使用minimal主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # 设置标题的大小和加粗
    axis.title.x = element_text(size = 16),  # 设置X轴标签的大小
    axis.title.y = element_text(size = 16),  # 设置Y轴标签的大小
    axis.text.x = element_text(size = 14),   # 设置X轴刻度文本的大小
    axis.text.y = element_text(size = 14),   # 设置Y轴刻度文本的大小
    legend.title = element_text(size = 16),  # 设置图例标题的大小
    legend.text = element_text(size = 14)    # 设置图例标签的大小
  )

# 显示图形
###########

sumGrna<-apply(reduced_mat,2,sum)
df <- data.frame(sumGrna = sumGrna)

# 绘制带有紫色、绿色、黄色渐变的直方图


# 创建直方图和箱线图结合的图形
p_combined2 <- ggplot(df, aes(x = sumGrna)) +
  # 绘制直方图
  geom_histogram(
    aes(fill = after_stat(count)),  # 将填充颜色映射到计数值
    binwidth = 1000,  # 设置柱子宽度
    color = "black",
    alpha = 0.6  # 设置透明度，方便看到箱线图
  ) +
  scale_fill_gradientn(
    colors = c("#FFD700", "green", "#6A0DAD")  # 设置颜色渐变
  ) +
  # 绘制箱线图（不显示异常值）
  geom_boxplot(
    aes(y = 11250),  # 将箱线图的位置放置在 y = 1850
    width = 800,  # 设置箱线图的宽度
    fill = "#FFD700",  # 设置箱线图的填充颜色为金色
    color = "darkblue",  # 设置箱线图的边框颜色为深金色
    notch = TRUE,  # 使用缺口
  ) +
  # 添加标题和标签
  labs(
    title = "Distribution of gRNA Numbers in each cell",
    x = "All gRNA Expression in each cell",
    y = "Cell Numbers",
    fill = "Cell Number"
  ) +
  # 设置x轴显示范围和刻度
  scale_x_continuous(
    limits = c(0, 70000),  # 设置x轴的显示范围
    breaks = seq(0, 70000, by = 10000)  # 设置x轴刻度每10000一个刻度
  ) +
  # 使用minimal主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # 设置标题的大小和加粗
    axis.title.x = element_text(size = 16),  # 设置X轴标签的大小
    axis.title.y = element_text(size = 16),  # 设置Y轴标签的大小
    axis.text.x = element_text(size = 14),   # 设置X轴刻度文本的大小
    axis.text.y = element_text(size = 14),   # 设置Y轴刻度文本的大小
    legend.title = element_text(size = 16),  # 设置图例标题的大小
    legend.text = element_text(size = 14)    # 设置图例标签的大小
  )

# 显示图形


########

sumGrna <- colSums(reduced_mat != 0)
df <- data.frame(sumGrna = sumGrna)

# 绘制带有紫色、绿色、黄色渐变的直方图
p_combined3 <- ggplot(df, aes(x = sumGrna)) +
  # 绘制直方图
  geom_histogram(
    aes(fill = after_stat(count)),  # 将填充颜色映射到计数值
    binwidth = 2,  # 设置柱子宽度
    color = "black",
    alpha = 0.6  # 设置透明度，方便看到箱线图
  ) +
  scale_fill_gradientn(
    colors = c("#FFD700", "green", "#6A0DAD")  # 设置颜色渐变
  ) +
  # 绘制箱线图（不显示异常值）
  geom_boxplot(
    aes(y = 24000),  
    width = 1000,  # 设置箱线图的宽度
    fill = "#FFD700",  # 设置箱线图的填充颜色为金色
    color = "darkblue",  
    notch = TRUE,  # 使用缺口
  ) +
  # 添加标题和标签
  labs(
    title = "Distribution of gRNA type coverage",
    x = "Cell Number",
    y = "gRNA Numbers",
    fill = "gRNA Number"
  ) +
  # 设置x轴显示范围和刻度
  scale_x_continuous(
    limits = c(0, 20),  # 设置x轴显示范围
    breaks = seq(0, 20, by = 2)  # 设置x轴刻度每30一个刻度
  ) +
  # 使用minimal主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # 设置标题的大小和加粗
    axis.title.x = element_text(size = 16),  # 设置X轴标签的大小
    axis.title.y = element_text(size = 16),  # 设置Y轴标签的大小
    axis.text.x = element_text(size = 14),   # 设置X轴刻度文本的大小
    axis.text.y = element_text(size = 14),   # 设置Y轴刻度文本的大小
    legend.title = element_text(size = 16),  # 设置图例标题的大小
    legend.text = element_text(size = 14)    # 设置图例标签的大小
  )
###########
library(patchwork)
sumGrna <- rowSums(reduced_mat != 0)
df <- data.frame(sumGrna = sumGrna)

# 绘制带有紫色、绿色、黄色渐变的直方图
p_combined4 <- ggplot(df, aes(x = sumGrna)) +
  # 绘制直方图
  geom_histogram(
    aes(fill = after_stat(count)),  # 将填充颜色映射到计数值
    binwidth = 2,  # 设置柱子宽度
    color = "black",
    alpha = 0.6  # 设置透明度，方便看到箱线图
  ) +
  scale_fill_gradientn(
    colors = c("#FFD700", "green", "#6A0DAD")  # 设置颜色渐变
  ) +
  # 绘制箱线图（不显示异常值）
  geom_boxplot(
    aes(y = 3900),  
    width = 200,  # 设置箱线图的宽度
    fill = "#FFD700",  # 设置箱线图的填充颜色为金色
    color = "darkblue",  
    notch = TRUE,  # 使用缺口
  ) +
  # 添加标题和标签
  labs(
    title = "Distribution of gRNA type coverage",
    x = "Cell Number",
    y = "gRNA Numbers",
    fill = "gRNA Number"
  ) +
  # 设置x轴显示范围和刻度
  scale_x_continuous(
    limits = c(0, 50),  # 设置x轴显示范围
    breaks = seq(0, 50, by = 5)  # 设置x轴刻度每30一个刻度
  ) +
  # 使用minimal主题
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # 设置标题的大小和加粗
    axis.title.x = element_text(size = 16),  # 设置X轴标签的大小
    axis.title.y = element_text(size = 16),  # 设置Y轴标签的大小
    axis.text.x = element_text(size = 14),   # 设置X轴刻度文本的大小
    axis.text.y = element_text(size = 14),   # 设置Y轴刻度文本的大小
    legend.title = element_text(size = 16),  # 设置图例标题的大小
    legend.text = element_text(size = 14)    # 设置图例标签的大小
  )

# 显示图形

(p_combined+p_combined2)/(p_combined3+p_combined4)
#############################
#############################
#############################
gc<-which(reduced_mat[10790,]!=0)
reduced_matt<-reduced_mat[-10790,gc]
zero_columns <- which(colSums(reduced_matt != 0) == 0)
cell_1<-hsv_a$cell
gc<-subset(hsv_a,cell_1 %in% colnames(reduced_matt)[zero_columns],value)
######对照矩阵
#####感染水平
p2<-gc$value#对照细胞感染情况
cell_non_g<- ifelse(colnames(reduced_mat) %in% rownames(gc), "NT","non")
cell_non_g<-t(cell_non_g)
cell_non_g<-as.character(cell_non_g)
idx = 0
find_g<-function(x){
  idx<<-idx+1
  k<-which(x!=0)#导入gRNA的细胞位置
  if (rownames(reduced_mat)[idx] %in% rownames(object[["RNA"]]) && length(k)>4){
    p1<-hsv_a$value[k]#导入gRNA的细胞感染情况
    hsv_p<-wilcox.test(p1,p2)
    cell_non_g[k]<-"gRNA"
    object <- AddMetaData(object, metadata = cell_non_g, col.name = "grna")
    gRNA_p<-FindMarkers(object,logfc.threshold = 0.25,only.pos=T,features  = rownames(reduced_mat)[idx],ident.1 ="NT" ,ident.2 ="gRNA",group.by = "grna")
    if(hsv_p$p.value<0.05 && nrow(gRNA_p)!=0) {
      gRNA_p$HSV_median<-median(p1)
      gRNA_p$wilcox_p<-hsv_p$p.value
      gRNA_p$cellNum<-length(p1)
      gRNA_p$infPct<-length(p1 != 0)/length(p1)
      print(idx)
      return(gRNA_p)
    }
  }
}
result<-apply(reduced_mat,1,find_g)
result <- do.call(rbind, result)

infP<-c()
for (i in c(1:105)){ 
  k<-which(reduced_matt[i,]!=0)
  p1<-hsv_a$value[k]
  infP<-c(infP,length(which(p1 != 0))/length(p1))
}
result$infPct<-infP

############################
############################
#############################
gc<-which(reduced_mat[10790,]!=0)
reduced_matt<-reduced_mat[-10790,gc]
zero_columns <- which(colSums(reduced_matt != 0) == 0)
cell_1<-hsv_a$cell
gc<-subset(hsv_a,cell_1 %in% colnames(reduced_matt)[zero_columns],value)

######对照矩阵
#####感染水平
p2<-gc$value#对照细胞感染情况
table(gc$value==0)


##########################
# 加载必要的 R 包
library(limma)

# 读取 FPKM 数据矩阵（假设行是基因，列是样本）
geo <- read.csv("GSE249664_THP-1_HSV-1_24H_fpkm.csv", row.names = 1)
rownames(geo)<-geo[,1]
geo<-geo[,-1]
geo<-as.matrix(geo)
library(limma)

# 假设 geo 是包含 FPKM 的基因表达数据，行为基因，列为样本
# 对数据进行预处理和标准化
geo_filtered <- geo[rowMeans(geo) > 1, ]  # 筛选出平均表达值大于 1 的基因，减少低表达的噪声基因
geo_normalized <- normalizeBetweenArrays(geo_filtered, method = "quantile")  # 进行样本间标准化

# 创建分组信息
group_geo <- factor(c("HC", "HC", "HSV", "HSV"))

# 创建设计矩阵
design <- model.matrix(~ 0 + group_geo)
colnames(design) <- levels(group_geo)

# 对数据进行差异表达分析
fit <- lmFit(geo_normalized, design)

# 设置对比：HSV 组 vs. HC 组
contrast_matrix <- makeContrasts(HSV_vs_HC = HSV - HC, levels = design)

# 应用对比矩阵
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_genes <- topTable(fit2, adjust = "BH", number = 100, sort.by = "logFC")

# 筛选出上调基因（logFC > 0）
hsvg_g <- rownames(all_genes[all_genes$logFC > 0, ])
################################################
library(AUCell)
expr_matrix <- as.matrix(object@assays$RNA@layers$data)
rownames(expr_matrix)<-rownames(object@assays$RNA)
colnames(expr_matrix)<-colnames(object@assays$RNA)
#gene_sets <- list(hsvg = hsvg,cgas=cgas,rece=rece,apop=apop,jak=jak,nlrp3=nlrp3,rigi=rigi,tlr=tlr,nfkb=nfkb)
gene_sets <- list(hsvg = hsvg_g)
# Step 1: 构建排名矩阵
rankings <- AUCell_buildRankings(expr_matrix, nCores = 1, plotStats = TRUE)

# Step 2: 计算基因集的 AUC
auc_scores <- AUCell_calcAUC(gene_sets, rankings,aucMaxRank = nrow(expr_matrix) * 0.1)

object@meta.data$hsvg<- as.numeric(getAUC(auc_scores)["hsvg", ])
ggplot(object@meta.data, aes(x =Infection_Level , y =hsvg , fill = Infection_Level)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "ISG Gene set  Score  ", x = "Group", y = "Score")
########################################################


hsvg_g<-intersect(hsvg_g,rownames(expr_matrix))
gene_sets <- list(hsvg = hsvg_g)
auc_scores <- AUCell_calcAUC(gene_sets, rankings,aucMaxRank = nrow(expr_matrix)*0.05)

object@meta.data$hsvg<- as.numeric(getAUC(auc_scores)["hsvg", ])
ggplot(object@meta.data, aes(x =Infection_Level , y =hsvg , fill = Infection_Level)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "hsvg AUC Score", x = "Group", y = "AUC Score")
