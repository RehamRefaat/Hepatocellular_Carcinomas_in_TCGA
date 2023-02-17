library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-LIHC")
query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")
lihc_res = getResults(query_TCGA) # make results as table
head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table
############################
head(lihc_res$sample_type)
summary(lihc_res$sample_type)
###################################
query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary solid Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)
dim(tcga_data)
colnames(colData(tcga_data))
table(tcga_data@colData$vital_status)
table(tcga_data@colData$tumor_stage)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)
dim(assay(tcga_data))
head(assay(tcga_data)[,1:10])
head(rowData(tcga_data))
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

tcga_data = readRDS(file = "tcga_data.RDS")
########################################################3
limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}
###########################################################################
limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)
########################################################################
# Save The File 
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)
#############################
gender_limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="gender",
  reference_group="female"
)

#####################################################
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("topright", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}
res_pca = plot_PCA(limma_res$voomObj, "definition")
#############################################################################
# Transpose and make it into a matrix object
d_mat = as.matrix(t(limma_res$voomObj$E))

# As before, we want this to be a factor
d_resp = as.factor(limma_res$voomObj$targets$definition)

# Divide data into training and testing set

# Set (random-number-generator) seed so that results are consistent between runs
set.seed(42)
train_ids = createDataPartition(d_resp, p=0.75, list=FALSE)

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

y_train = d_resp[train_ids]
y_test  = d_resp[-train_ids]
x = x_train
y = y_train
res = cv.glmnet(
  x,
  y ,
  alpha = 0.5,
  family="binomial"
)

y_pred = predict(res, newx=x_test, type="class", s="lambda.min")

confusion_matrix = table(y_pred, y_test)

# Evaluation statistics
print(confusion_matrix)

print(paste0("Sensitivity: ",sensitivity(confusion_matrix)))
print(paste0("Specificity: ",specificity(confusion_matrix)))
print(paste0("Precision: ",precision(confusion_matrix)))

res_coef = coef(res, s="lambda.min") # the "coef" function returns a sparse matrix
dim(res_coef$`Primary solid Tumor`)
head(res_coef)
res_coef = res_coef[res_coef[,1] != 0,]
# note how performing this operation changed the type of the variable
head(res_coef)

res_coef = res_coef[-1]

relevant_genes = names(res_coef) # get names of the (non-zero) variables.
length(relevant_genes) # number of selected genes
head(relevant_genes)

head(limma_res$voomObj$genes)

relevant_gene_names = limma_res$voomObj$genes[relevant_genes,"external_gene_name"]

head(relevant_gene_names) # few select genes (with readable names now)

print(intersect(limma_res$topGenes$ensembl_gene_id, relevant_genes))

########################################################################################
# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)

# Show green color for genes that also show up in DE analysis
colorLimmaGenes = ifelse(
  # Given a vector of boolean values
  (relevant_genes %in% limma_res$topGenes$ensembl_gene_id),
  "green", # if true, return green for that value
  "white" # if false, return white for that value
)

# As you've seen a good looking heatmap involves a lot of parameters
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labRow=relevant_gene_names, # use gene names instead of ensembl annotation
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)
##########################################################################



# Using the same method as in Day-2, get the dendrogram from the heatmap
# and cut it to get the 2 classes of genes

# Extract the hierarchical cluster from heatmap to class "hclust"
hc = as.hclust(gene_heatmap$rowDendrogram)

# Cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters = cutree(hc, k=2)
table(clusters)



# selecting just a few columns so that its easier to visualize the table
gprofiler_cols = c("significant","p.value","overlap.size","term.id","term.name")

# make sure the URL uses https
set_base_url("https://biit.cs.ut.ee/gprofiler")

# Group 1, up in tumor
gprofiler(names(clusters[clusters %in% 1]))[, gprofiler_cols]
install.packages("gprofiler2")
library("gprofiler2")


########################################################################################


library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
tcga_data = readRDS(file = "tcga_data.RDS")
limma_res = readRDS(file = "limma_res.RDS")

clinical = tcga_data@colData

dim(clinical)

clin_df = clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "ajcc_pathologic_stage")]


# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)

Surv(clin_df$overall_survival, clin_df$deceased)
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)

print(fit)

fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)

print(fit)

ggsurvplot(fit, data=clin_df)
ggsurvplot(fit, data=clin_df, pval=T)
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, risk.table.col="strata")
# remove any of the letters "a", "b" or "c", but only if they are at the end
# of the name, eg "stage iiia" would become simply "stage iii"
clin_df$tumor_stage = gsub("[abc]$", "", clin_df$ajcc_pathologic_stage)

# we remove those with stage "not reported", since they are unknown
clin_df[which(clin_df$tumor_stage == "not reported"), "ajcc_pathologic_stage"] = NA

# finally, we also remove those with tumor stage 4, since they are too few
clin_df[which(clin_df$tumor_stage == "stage iv"), "ajcc_pathologic_stage"] = NA

table(clin_df$tumor_stage)
############################################
######################################





fit = survfit(Surv(overall_survival, deceased) ~ tumor_stage, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T)


# let's extract the table of differential expression we got earlier
expr_df = limma_res$topGenes

# print the first row, to see the gene name, the logFC value and the p-value
print(expr_df[1, ])

# get the ensembl gene id of the first row
gene_id = expr_df[1, "ensembl_gene_id"]

# also get the common gene name of the first row
gene_name = expr_df[1, "external_gene_name"]




# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)
expr_diseased = d_mat[rownames(clin_df), gene_id]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")



# get the expression values for the selected gene
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]

# find the median value of the gene and print it
median_value = median(clin_df$gene_value)
print(median_value)




# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it
# among the "up-regulated", otherwise among the "down-regulated"
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)


# and finally, we produce a Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, title=paste(gene_name))
