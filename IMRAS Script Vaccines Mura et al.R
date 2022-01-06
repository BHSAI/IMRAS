---
  title: 'Immunoprofiling identifies functional B and T cell subsets induced by an attenuated whole parasite malaria vaccine as correlates of sterile protection'
author: "Marie Mura, Pinyi Lu, Tanmaya Atre, Jessica B. Bolton, Elizabeth H. Duncan, Sidharta Chaudhury, Elke S. Bergmann-Leitner"
date: "1 décembre 2021"
---


# Define the function pkgTest to check if a package has been installed.
## If it hasn't, then install and load it.
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

##Install binary versions of packages
pkgTest("plyr")
pkgTest("caret")
pkgTest("corrplot")
pkgTest("RColorBrewer")
pkgTest("tinytex")
pkgTest("readxl")
pkgTest("tidyverse")
pkgTest("xlsx")
pkgTest("ggplot2")
pkgTest("ggfortify")
pkgTest("FactoMineR")
pkgTest("factoextra")

# Choose work directory and download files

setwd("~/EN COURS/IMRAS_Rstudio/Data")

load("~/EN COURS/IMRAS_Rstudio/data/Tfh_output.rda")
load("~/EN COURS/IMRAS_Rstudio/data/subject_info.rda")
load("~/EN COURS/IMRAS_Rstudio/data/Serology_MSD.rda")
load("~/EN COURS/IMRAS_Rstudio/data/Serology_ELISA.rda")
load("~/EN COURS/IMRAS_Rstudio/data/CYTOK_output.rda")
load("~/EN COURS/IMRAS_Rstudio/data/CD8_output.rda")
load("~/EN COURS/IMRAS_Rstudio/data/CD4_output.rda")
load("~/EN COURS/IMRAS_Rstudio/data/Bflow_output.rda")
load("~/EN COURS/IMRAS_Rstudio/data/BELI_output.rda")


# remove mock subjects
BELI_output <- BELI_output[c(1:5,7,9:18),]
CYTOK_output<-CYTOK_output[c(1:5,7,9:14),]
subject_info <- subject_info[c(1:5,7,9:18),]

# combined all dataset
dataset <- list(Tfh_output, Serology, MSD, Bflow_output, BELI_output,CYTOK_output, CD4_output, CD8_output)
combinedData <- purrr::reduce(dataset, dplyr::full_join, by = "Subject_ID")
combinedData <- merge (combinedData, subject_info, by = "Subject_ID")

# univariate analysis: compare T0 vs T1
df.list <- list(BELI_output,CYTOK_output, Serology, MSD, Tfh_output, CD4_output, CD8_output, Bflow_output)

comput.p <- function(X,Y,a){
  if(diff(range(X)) == 0 | diff(range(Y)) == 0 ){
    return(wilcox.test(X, Y, 
                       paired=TRUE, 
                       conf.level=0.95)$p.value)
  }
  else if(shapiro.test(X)$p.value > a & shapiro.test(Y)$p.value > a){
    return(t.test(X, Y, 
                  paired=TRUE, 
                  conf.level=0.95)$p.value)
  }else{
    return(wilcox.test(X, Y, 
                       paired=TRUE, 
                       conf.level=0.95)$p.value)
  }
}

results <- data.frame(matrix(ncol = 2, nrow = 0))
names(results) <- c("name","pvalue")
time_point <- 2
for(z in 1:length(df.list)){
  l <- which(grepl("T0",names(df.list[[z]])))
  for (i in l) {
    for (j in 1:(time_point- 1)){
      p <- comput.p(df.list[[z]][,i],df.list[[z]][,(i+j)],0.05)
      results[nrow(results) + 1,] = list(names(df.list[[z]])[i+j],p)
    }
  }
}
results$BH = p.adjust(results$pvalue, method = "BH")
results
VIR_complete1 <- results[which(results$pvalue < 0.05 & results$BH < 0.05),]
VIR_complete2 <- results[which(results$pvalue < 0.05 & results$BH < 0.2),]

# build matrice of  correlation
M <-cor(combinedData[,VIR_complete1$name], use="pairwise.complete.obs")

# Add significance to the matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(combinedData[,VIR_complete1$name])

rownames(M)<-c("B ELIspot.AMA",      "B ELIspot.CSP",    "B ELIspot.MSP",           
               "B ELIspot.SPZ",       "IFA.SPZ",             "IL8.AMA1",            
               "ELISA.AMA1" ,         "ELISA.CSP_FL" ,
               "MSD.CSP_NANP"  ,       "cTfh CCR6+CXCR3+.CSP",     
               "cTfh CXCR3+TNFa+.M",        "T CD4+CXCR3+TNFa+IFNg+.CSP" ,       
               "T CD4+CXCR3+TNFa+.SPZ",     "T CD4+CXCR3-TNFa+IFNg+.CSP",        
               "T CD4+CXCR3-TNFa+IFNg+.SPZ", "T CD4+CXCR3-TNFa+.SPZ",             
               "T CD4+TNFa+IFNg+.CSP",     "lymphocytes.CSP",
               "lymphocytes.M",          "lymphocytes.SPZ",
               "BC.SPZ" , "Naive_BC.SPZ",     "MBC.SPZ" ,"IgD-MBC.SPZ",  
               "IgD+MBC.SPZ","transitional_BC.SPZ")
colnames(M)<-c("B ELIspot.AMA",      "B ELIspot.CSP",    "B ELIspot.MSP",           
               "B ELIspot.SPZ",       "IFA.SPZ",             "IL8.AMA1",            
               "ELISA.AMA1" ,         "ELISA.CSP_FL" ,
               "MSD.CSP_NANP"  ,       "cTfh CCR6+CXCR3+.CSP",     
               "cTfh CXCR3+TNFa+.M",        "T CD4+CXCR3+TNFa+IFNg+.CSP" ,       
               "T CD4+CXCR3+TNFa+.SPZ",     "T CD4+CXCR3-TNFa+IFNg+.CSP",        
               "T CD4+CXCR3-TNFa+IFNg+.SPZ", "T CD4+CXCR3-TNFa+.SPZ",             
               "T CD4+TNFa+IFNg+.CSP",     "lymphocytes.CSP",
               "lymphocytes.M",          "lymphocytes.SPZ",
               "BC.SPZ" , "Naive_BC.SPZ",     "MBC.SPZ" ,"IgD-MBC.SPZ",  
               "IgD+MBC.SPZ","transitional_BC.SPZ")

corrplot::corrplot(M, type="upper", order="hclust", tl.cex=0.45,cl.cex=0.8,tl.srt=45,
                   tl.col="Black",col=RColorBrewer::brewer.pal(n=8, name="RdYlBu"),p.mat = p.mat, sig.level = 0.05, insig = "blank")

dev.print(jpeg,filename="CorMat vaccine-incuced factors.jpg",quality=100,units="px", width=2000,res=300)

# prepare dataset of significant factors for humoral and cellular response

ImmuneSign<-combinedData[c(1:9,13:15),c(VIR_complete1$name,"Subject_ID")]

matrixT1<-combinedData[c(1:9,13:15),c(VIR_complete2$name,"Subject_ID")]
match(VIR_complete2$name, names(combinedData))
matrixT0<-combinedData[c(1:9,13:15),c(1,362,  364,  366,  368,  370,  372,  382,  404,  418,  420,  422,  462,  464,  476,  486,  488,  170,  174,  178,   10,   18,   38,   66,  122,  126,  146,  162,
                        518,  554,  586,  598,  610,  618,  630,  650,  694,  698,  706,  714,  726,  742,  750,  762,  779,  783,  787,  803,  819,  855,  867,  879,  935,  943, 1043,
                        1047,  214,  238,  262,  286,  310,  334)]

matrix <- merge (matrixT0, matrixT1, by = "Subject_ID")

#Order data by column names
matrix[, 2:ncol(matrix)] <- matrix[, 2:ncol(matrix)][, order(names(matrix[, 2:ncol(matrix)]))]
colnames(matrix)[2:ncol(matrix)] <- names(matrix[, 2:ncol(matrix)][, order(names(matrix[, 2:ncol(matrix)]))])

combinedData <- matrix

# univariate analysis: compare cohort 1 and 2
comput.p <- function(X,Y,a){
  combinedData <- c(X,Y)
  groups <- rep(
    c(1,2),
    times = c(length(X),length(Y))
  )
  if(diff(range(X, na.rm = TRUE)) == 0 | diff(range(Y, na.rm = TRUE)) == 0 ){
    return(kruskal.test(combinedData, groups)$p.value)
  }
  else if(shapiro.test(X)$p.value > a & shapiro.test(Y)$p.value > a){
    return(t.test(X, Y, 
                  paired=FALSE, 
                  conf.level=0.95)$p.value)
  }else{
    return(kruskal.test(combinedData, groups)$p.value)
  }
}


Subject_ID_One <- subject_info$Subject_ID[which(subject_info$Cohort == "One")]
Subject_ID_Two <- subject_info$Subject_ID[which(subject_info$Cohort == "Two")]

ImmuneMeasure_ID <- colnames(combinedData)[-1]
results <- data.frame(matrix(ncol = 2, nrow = length(ImmuneMeasure_ID )))
colnames(results) <- c("Parameter", "pvalue")
for (i in 1:length(ImmuneMeasure_ID)) {
  ONE <- combinedData[which(combinedData$Subject_ID %in% Subject_ID_One),ImmuneMeasure_ID[i]]
  TWO <- combinedData[which(combinedData$Subject_ID %in% Subject_ID_Two),ImmuneMeasure_ID[i]]
  results[i,1] <- ImmuneMeasure_ID[i]
  results[i,2] <- comput.p(as.numeric(ONE), as.numeric(TWO), 0.05)}
results$fdr <- p.adjust(results$pvalue, method = "BH")
results
VIR_complete <- results[which(results$pvalue < 0.05 & results$fdr < 0.05),]
VIR_complete

write.xlsx(VIR_complete, file="cohort_effect.xlsx")

Cohortdiff<-combinedData[,c(VIR_complete$Parameter,"Subject_ID")]

M2 <-cor(Cohortdiff[,1:16], use="pairwise.complete.obs")

# Add significance to the matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat2 <- cor.mtest(Cohortdiff[,1:16])

corrplot::corrplot(M2, type="upper", order="hclust", tl.cex=0.45,cl.cex=0.8,tl.srt=45,
                   tl.col="Black",col=RColorBrewer::brewer.pal(n=8, name="RdYlBu"),p.mat = p.mat, sig.level = 0.05, insig = "blank")



# PCA plot

X <- merge (combinedData, subject_info, by = "Subject_ID")

PCA(X[,c(2:49,51:123)], scale.unit = TRUE, ncp = 5, graph = TRUE)

res.pca <- PCA(X[,c(2:49,51:123)], graph = FALSE)
print(res.pca)

fviz_pca_ind (res.pca,
              col.ind = X$Cohort, palette = c("darkorange1","brown4"),
              addEllipses = TRUE, label = "var",
              col.var = "black", repel = TRUE,
              legend.title = "Cohort")

dev.print(jpeg,filename="PCA plot cohort 1vs2.jpg",quality=100,units="px", width=1200,res=300)


# compare P/NP from immunization induced factors
matrixCD4 <- matrix[,c(1,54:83)]
matrixCD8 <- matrix[,c(1,84:101)]
matrixHR <- matrix[,c(1:51)]

dataset <- list(matrixCD4,subject_info)
combinedData <- purrr::reduce(dataset, dplyr::full_join, by = "Subject_ID")


# univariate analysis: compute p to compare P/NP
comput.p <- function(X,Y,a){
  combinedData <- c(X,Y)
  groups <- rep(
    c(1,2),
    times = c(length(X),length(Y))
  )
  if(diff(range(X, na.rm = TRUE)) == 0 | diff(range(Y, na.rm = TRUE)) == 0 ){
    return(kruskal.test(combinedData, groups)$p.value)
  }
  else if(shapiro.test(X)$p.value > a & shapiro.test(Y)$p.value > a){
    return(t.test(X, Y, 
                  paired=FALSE, 
                  conf.level=0.95)$p.value)
  }else{
    return(kruskal.test(combinedData, groups)$p.value)
  }
}

Subject_ID_P <- subject_info$Subject_ID[which(subject_info$Status == "Protected")]
Subject_ID_NP <- subject_info$Subject_ID[which(subject_info$Status == "Not Protected")]

ImmuneMeasure_ID <- colnames(combinedData)[-1]
results <- data.frame(matrix(ncol = 2, nrow = length(ImmuneMeasure_ID )))
colnames(results) <- c("Parameter", "pvalue")
for (i in 1:length(ImmuneMeasure_ID)) {
  P <- combinedData[which(combinedData$Subject_ID %in% Subject_ID_P),ImmuneMeasure_ID[i]]
  NP <- combinedData[which(combinedData$Subject_ID %in% Subject_ID_NP),ImmuneMeasure_ID[i]]
  results[i,1] <- ImmuneMeasure_ID[i]
  results[i,2] <- comput.p(as.numeric(P), as.numeric(NP), 0.05)}
results$fdr <- p.adjust(results$pvalue, method = "BH")
results
VIR_complete <- results[which(results$pvalue < 0.05 & results$fdr < 0.2),]
VIR_complete


# Random Forest model

# Define features X and response variable Y
# with all sign factors from T1 (except MSD NANP because bias due to NA replacing)
subject_info<-subject_info[c(1:9,13:15),]

X <- matrixT1[,c(-19,-62)]
Y <- make.names(subject_info$Status)

X <- matrix[,c(79,95,25,61,7,89,78,94,75)]
Y <- make.names(subject_info$Status)

#Build random forest model, optimize hyper parameters using repeated cross validation
control <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 100,
  verboseIter = FALSE,
  returnData = FALSE,
  allowParallel = TRUE,
  #summaryFunction = twoClassSummary,
  classProbs = TRUE,
  savePredictions = T)

preProc <- preProcess(X)

set.seed(112358)
rfUnbalance <- caret::train(
  x = data.frame(predict(preProc, X)),
  y = Y,
  method = "rf", #random forest
  ntree = 1500,
  tuneLength = 10,
  metric = "Accuracy", #area under ROC curve
  trControl = control)

#Check 
rfUnbalance
rfUnbalance$finalModel

#Calculation of variable importance
varImportance <- varImp(rfUnbalance)$importance
varImportance[order(varImportance$Overall, decreasing = TRUE),,drop = FALSE]

# PCA 
X <- matrix[,c(1,79,95,25,61,89,78,94,75)]
X<-merge(X, subject_info,by="Subject_ID")

X[,10]<-c("P","NP","NP","NP","P","P","P","P","NP","P","P","P")
colnames(X)=c("Subject_ID","T CD4+TNFa+ SPZ.T1","T CD8+CCR6-CXCR3-TNFa+ CSP.T1", "transBC SPZ.T1",
              "T CD4+CXCR3+TNFa+ SPZ.T1","T CD8+CCR6+IFNg+ CSP.T1",                 
              "T CD4+TNFa+ SPZ.T0","T CD8+CCR6-CXCR3-TNFa+ CSP.T0","T CD4+TNFa+IFNg+ CSP.T1", "Status", "Cohort" )

PCA(X[,-c(1,10,11)], scale.unit = TRUE, ncp = 5, graph = TRUE)

res.pca <- PCA(X[,-c(1,10,11)], graph = FALSE)
print(res.pca)



fviz_pca_biplot (res.pca, geom="point", axes=c(1,2),
                 col.ind = X$Status, palette = c("red","blue"), mean.point=FALSE, pointsize=4,
                 addEllipses = FALSE, label = "var", 
                 col.var = "black", repel = TRUE,
                 legend.title = "Protection status") 
