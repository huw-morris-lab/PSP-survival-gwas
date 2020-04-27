# PSP-survival-gwas
PSP survival GWAS script written by Edwin Jabbari and Maryam Shoai.

library(data.table)
library(survival)
# Load .RData file which is the genetic dataset in RAWfile format merged with phenotype file (columns 1-9 in the following order: IID, Diseaseduration (in years), Event (yes or no), Richardson (yes or no), Gender (male or female), Ageonset (in years), PCA1 (value), PCA2 (value), PCA3 (value))
load("cxph_one.RData")

TABLE<- TABLE_one
coefficients<-as.data.frame(matrix(ncol= 9))

# Results output file will contain these columns
names(coefficients) <- c("SNP","Coeff", "se", "Pvalue", "Cox.zphPVal", "N", "ov.lik.ratio","logrank", "r2" )

# Run the Cox survival model for your dataset with disease duration as the outcome measure and death as the "Event" - using gender, age at onset, PCA 1-3 and PSP phenotype ("Richardson") as covariates
for (i in 10:ncol(TABLE)) {
  print(colnames(TABLE)[i])
  snp<- TABLE[,c(i,1:9)]
  model.cox<- coxph(Surv(snp$Diseaseduration, snp$Event) ~ snp[,1]+ snp$Ageonset + snp$Gender + snp$PCA1 + snp$PCA2 + snp$PCA3 + snp$Richardson, data=snp)  
  kmz<- cox.zph(model.cox, transform = "km")
  j= i-9
  coefficients[j,1]<- paste(colnames(TABLE)[i])
  coefficients[j,2]<- summary(model.cox)$coefficients[1,1] 
  coefficients[j,3]<- summary(model.cox)$coefficients[1,3] 
  coefficients[j,4]<- summary(model.cox)$coefficients[1,5] 
  coefficients[j,5]<- kmz$table[1,3]
  coefficients[j,6]<- model.cox$n
  coefficients[j,7]<- summary(model.cox)$logtest[[1]]
  coefficients[j,8]<- summary(model.cox)$sctest[[1]]
  coefficients[j,9]<- summary(model.cox)$rsq[[1]] # nagelkerke r square
  }
  
  # Write results
  fwrite(coefficients, "coefficients_results_1.txt", row.names=FALSE, sep="\t", quote= FALSE)
