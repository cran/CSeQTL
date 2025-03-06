## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,comment = "#>",
	echo = TRUE,cache = TRUE,
	dev = "png")

## ----setup,warning = FALSE----------------------------------------------------
# Load libraries
req_packs = c("ggplot2","smarter","CSeQTL")
for(pack in req_packs){
  library(package = pack,character.only = TRUE)
}

# List package's exported functions
ls("package:CSeQTL")

# Fix seed
set.seed(2)

## ----fig.dim = c(8,6)---------------------------------------------------------
# sample size
NN = 3e2

# fold-change between q-th and 1st cell type
true_KAPPA  = c(1,3,1)

# eQTL effect size per cell type, 
#   fold change between B and A allele
true_ETA    = c(1,1,1)

# cis/trans effect size
true_ALPHA  = c(1,1,1)

# count number of cell types
QQ = length(true_KAPPA)

# TReC model overdispersion
true_PHI = 0.1

# ASReC model overdispersion
true_PSI = 0.05

# Simulate cell type proportions
tRHO = gen_true_RHO(wRHO = 1,NN = NN,QQ = QQ)
plot_RHO(RHO = tRHO)
boxplot(tRHO,xlab = "Cell Type",ylab = "Proportion")

# Simulate a data object for a gene and SNP
sim = CSeQTL_dataGen(NN = NN,MAF = 0.3,true_BETA0 = log(1000),
  true_KAPPA = true_KAPPA,true_ETA = true_ETA,true_PHI = true_PHI,
  true_PSI = true_PSI,prob_phased = 0.05,true_ALPHA = true_ALPHA,
  RHO = tRHO,cnfSNP = TRUE,show = FALSE)

names(sim)

# TReC, ASReC, haplotype 2 counts
sim$dat[1:3,]

# Batch covariates including the intercept
sim$XX[1:3,]

sim$dat$SNP = sim$true_SNP
sim$dat$SNP = factor(sim$dat$SNP,
  levels = sort(unique(sim$dat$SNP)),
  labels = c("AA","AB","BA","BB"))

# TReC vs SNP
ggplot(data = sim$dat,
  mapping = aes(x = SNP,y = log10(total + 1))) +
  geom_violin(aes(fill = SNP)) + geom_jitter(width = 0.25) +
  geom_boxplot(width = 0.1) +
  xlab("Phased Genotype") + ylab("log10(TReC + 1)") +
  theme(legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12))

# TReC vs ASReC
ggplot(data = sim$dat,
  mapping = aes(x = total,y = total_phased,color = SNP)) +
  geom_point() + xlab("Total Read Count (TReC)") + 
  ylab("Total Allele-specific Read Count (ASReC)") +
  theme(legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12))

## ----R.options = list(width = 200)--------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(sim$XX)

mRHO = matrix(1,NN,1)
colnames(mRHO) = "Bulk"
cistrans_thres = 0.01
res = c()

for(trec_only in c(TRUE,FALSE)){
for(neg_binom in c(TRUE,FALSE)){
for(beta_binom in c(TRUE,FALSE)){
  
  cat(sprintf("%s: trec_only = %s, neg_binom = %s, beta_binom = %s ...\n",
    date(),trec_only,neg_binom,beta_binom))
  
  PHASE = sim$dat$phased * ifelse(trec_only,0,1)
  upPHI = ifelse(neg_binom,1,0)
  upPSI = ifelse(beta_binom,1,0)
  
  sout = CSeQTL_smart(TREC = sim$dat$total,hap2 = sim$dat$hap2,
    ASREC = sim$dat$total_phased,PHASE = PHASE,
    SNP = sim$true_SNP,RHO = mRHO,XX = sim$XX,upPHI = upPHI,
    upKAPPA = 1,upETA = 1,upPSI = upPSI,upALPHA = 1,
    iFullModel = FALSE,trim = FALSE,thres_TRIM = 10,
    hypotest = TRUE,swap = FALSE,numAS = 5,numASn = 5,
    numAS_het = 5,cistrans_thres = cistrans_thres)
  # sout$HYPO
  
  res = rbind(res,smart_df(Model = ifelse(trec_only,"TReC-only","TReCASE"),
    TReC_Dist = ifelse(neg_binom,"Negative Binomial","Poisson"),
    ASReC_Dist = ifelse(beta_binom,"Beta-Binomial","Binomial"),
    sout$res))
  rm(sout)
  
}}}

num_vars = which(unlist(lapply(res,function(xx) class(xx))) == "numeric")
res[,num_vars] = apply(res[,num_vars],2,function(xx) smart_SN(x = xx,digits = 2))
res$Interpret = apply(res[,c("cistrans","PVAL_eQTL")],1,function(xx){
  ct = xx[1]
  pval = as.numeric(xx[2])
  ifelse(pval < 0.05,sprintf("%s eQTL",ct),"no eQTL")
})
res$Correct_Model = apply(res[,c("TReC_Dist","ASReC_Dist")],1,function(xx){
  chk_trec = (( true_PHI > 0 & xx[1] == "Negative Binomial" )
    | ( true_PHI == 0 & xx[1] == "Poisson" ))
  chk_trec
  
  chk_asrec = (( true_PSI > 0 & xx[2] == "Beta-Binomial" )
    | ( true_PSI == 0 & xx[2] == "Binomial" ))
  
  chk_trec = ifelse(chk_trec,"TReC-Yes","TReC-No")
  chk_asrec = ifelse(chk_asrec,"ASReC-Yes","ASReC-No")
  sprintf("%s;%s",chk_trec,chk_asrec)
})

# Simplified output
smart_rmcols(res,c("LRT_trec","LRT_trecase","LRT_cistrans","cistrans"))

## ----R.options = list(width = 200)--------------------------------------------------------------------------------------------------------------------------------------------------------------------
cistrans_thres = 0.01
res = c()

for(trec_only in c(TRUE,FALSE)){
for(neg_binom in c(TRUE,FALSE)){
for(beta_binom in c(TRUE,FALSE)){
  
  cat(sprintf("%s: trec_only = %s, neg_binom = %s, beta_binom = %s ...\n",
    date(),trec_only,neg_binom,beta_binom))
  
  PHASE = sim$dat$phased * ifelse(trec_only,0,1)
  upPHI = ifelse(neg_binom,1,0)
  upPSI = ifelse(beta_binom,1,0)
  upKAPPA = rep(1,QQ)
  upETA = upKAPPA
  upALPHA = upETA
  
  sout = CSeQTL_smart(TREC = sim$dat$total,hap2 = sim$dat$hap2,
    ASREC = sim$dat$total_phased,PHASE = PHASE,
    SNP = sim$true_SNP,RHO = sim$true_RHO,XX = sim$XX,upPHI = upPHI,
    upKAPPA = upKAPPA,upETA = upETA,upPSI = upPSI,upALPHA = upALPHA,
    iFullModel = FALSE,trim = FALSE,thres_TRIM = 10,
    hypotest = TRUE,swap = FALSE,numAS = 5,numASn = 5,
    numAS_het = 5,cistrans_thres = cistrans_thres)
  # sout$HYPO
  
  res = rbind(res,smart_df(Model = ifelse(trec_only,"TReC-only","TReCASE"),
    TReC_Dist = ifelse(neg_binom,"Negative Binomial","Poisson"),
    ASReC_Dist = ifelse(beta_binom,"Beta-Binomial","Binomial"),
    sout$res))
  rm(sout)
  
}}}

num_vars = which(unlist(lapply(res,function(xx) class(xx))) == "numeric")
res[,num_vars] = apply(res[,num_vars],2,function(xx) smart_SN(x = xx,digits = 2))
res$Interpret = apply(res[,c("cistrans","PVAL_eQTL")],1,function(xx){
  ct = xx[1]
  pval = as.numeric(xx[2])
  ifelse(pval < 0.05,sprintf("%s eQTL",ct),"no eQTL")
})
res$Correct_Model = apply(res[,c("TReC_Dist","ASReC_Dist")],1,function(xx){
  chk_trec = (( true_PHI > 0 & xx[1] == "Negative Binomial" )
    | ( true_PHI == 0 & xx[1] == "Poisson" ))
  chk_trec
  
  chk_asrec = (( true_PSI > 0 & xx[2] == "Beta-Binomial" )
    | ( true_PSI == 0 & xx[2] == "Binomial" ))
  
  chk_trec = ifelse(chk_trec,"TReC-Yes","TReC-No")
  chk_asrec = ifelse(chk_asrec,"ASReC-Yes","ASReC-No")
  sprintf("%s;%s",chk_trec,chk_asrec)
})

# Simplified output
smart_rmcols(res,c("LRT_trec","LRT_trecase","LRT_cistrans","cistrans"))


## ----ols_bulk-----------------------------------------------------------------
ols_out = CSeQTL_linearTest(input = sim$dat,XX = sim$XX,
  RHO = sim$true_RHO,SNP = sim$true_SNP,MARG = TRUE)
names(ols_out)

# Model estimates
summary(ols_out$lm_out)

# Effect sizes and hypothesis test
ols_out$out_df

## ----ols_cts,fig.dim = c(6,7)-------------------------------------------------
ols_out = CSeQTL_linearTest(input = sim$dat,XX = sim$XX,
  RHO = sim$true_RHO,SNP = sim$true_SNP,MARG = FALSE)
names(ols_out)

# Model estimates
summary(ols_out$lm_out)

# Effect sizes and hypothesis tests
ols_out$out_df

## -----------------------------------------------------------------------------
sessionInfo()

