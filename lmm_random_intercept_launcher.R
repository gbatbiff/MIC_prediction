rm(list = ls())
library(rstan)
library(dplyr)
library(LearnBayes)
library(bayesplot)
library(OneR)
library(scales)
###
args=commandArgs(trailingOnly = TRUE)

phen_file <- args[1] ### simulated phenotypes
matrix_file<-args[2] ### pres/abs matrix
lineages_file<-args[3] ### poppunk lineages
stan_model<-args[4]

dd<-read.delim(matrix_file,sep='\t') 
#colnames(dd)<-gsub("X","",colnames(dd))
#colnames(dd)<-gsub(".fna","",colnames(dd))
rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))
tdd[tdd>1]<-1
phen<-read.delim(phen_file,sep='',header = F)

#phen$V1<-NULL
#phen$V2<-gsub("X","",phen$V2)
#phen$V2<-gsub(".fna","",phen$V2)
colnames(phen)[1]<-"ID"
colnames(phen)[2]<-"MIC"
###
tdd$ID<-rownames(tdd)

mm<-merge(tdd,phen,by="ID")
rownames(mm)<-mm$ID
#mm$ID<-NULL

mm$sim_cont_mic<-rescale(mm$MIC,to=c(0.25,16))
mm$MIC<-NULL


st<-read.delim(lineages_file,header = F)
st$V2<-gsub("-","_",st$V2)
st$V1<-as.character(st$V1)
st$V2<-as.character(st$V2)
colnames(st)<-c("ID","PP")

mm_clust<-merge(mm,st,by="ID")

rownames(mm_clust)<-mm_clust$ID

mm_clust$ID<-NULL
snps<-grep("rs",colnames(mm_clust),value = T)
mm_clust$PP<-gsub("PP","",mm_clust$PP)


for (i in snps){
  
  y<-mm_clust$sim_cont_mic
  x<-mm_clust[,i] 
  PP<-as.integer(mm_clust$PP)
  J<-nlevels(as.factor(mm_clust$PP))
  N<-nrow(mm_clust)
  
  
  stanDat_mix <- list(y = y, x = x, N = N,PP=PP,
                      J = J))
  
  ranIntFit <- stan(file = stan_model, data = stanDat_mix,
                    iter = 2000, chains = 4)
  

  
  df_diag<-as.data.frame(summary(ranIntFit)$summary)
  df_diag$snp<-i
  
  write.table(df_diag,sprintf("model_diagn_snp_%s",i),quote = F,sep='\t')
  
}


# 
# traceplot(ranIntFit)
# 
# mcmc_areas(ranIntFit, pars = c("beta[1]", "beta[2]", "sigma_u"), prob = 0.89) +
#   scale_y_discrete(expand = c(0, 0))
# 
# samples <- as.data.frame(rstan::extract(ranIntFit, pars = c("beta[1]", "beta[2]"))) %>%
#   rename(intercept = beta.1., slope = beta.2.)
# 
# ggplot(samples) +
#   geom_abline(aes(intercept = intercept, slope = slope), color = "cadetblue4", alpha = 1) +
#   geom_abline(aes(intercept = mean(intercept), slope = mean(slope)), color = "black", size = 2) +
#   scale_x_continuous(limits = c(-1, 1)) +
#   scale_y_continuous(limits = c(-0.45, 0.45)) +
#   labs(x = "x",
#        y = "y",
#        title = "Samples from beta posteriors")
# 

