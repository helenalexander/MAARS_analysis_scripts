# Show that microbes excluded by data filtering process improves the quality of our dataset
#1. Filtered dataset variance is unchanged from unfiltered dataset
#2. Filtered dataset diversity is unchanged from unfiltered dataset
#3. Filtered dataset is more ‘skin-like’ than unfiltered dataset (and not more like other microbial environments)

rm(list=ls())

#load filtered and unfiltered abundance data
load ("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/data-filtering/20160304met_corr_sample_filt_clades_abund.Rdata")
load ("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/data-filtering/20160304metaphlan_unfilt_clades.Rdata")
load ("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/data-filtering/20160304met_corr_filt_clades.Rdata")

# 1. Variance
Vunfilt = var(met_unfilt_clades)
VU = sum(abs(Vunfilt))
VU
Vcorfilt = var(met_corr_filt_clades)
VC = sum(abs(Vcorfilt))
VC
Vcorsamfilt = var(met_corr_sample_filt_clades_abund)
VCS = sum(abs(Vcorsamfilt))
VCS
boxplot (VU, VC, VCS, ylim=c(0, 226000), ylab = "variance", names = c("unfiltered", "step 1", "step 2"))

# 2. diversity

library("vegan")
dim(met_corr_filt_clades)

wSpecU = grep("^S_", (colnames(met_unfilt_clades)))
speciesU = met_unfilt_clades[,wSpecU]
wSpecCF  = grep("^S_", (colnames(met_corr_filt_clades)))
speciesCF = met_corr_filt_clades[,wSpecCF]
wSpecCSF  = grep("^S_", (colnames(met_corr_sample_filt_clades_abund)))
speciesCSF = met_corr_sample_filt_clades_abund[,wSpecCSF]
DivSpecUnfilt = diversity(speciesU, index = "shannon", MARGIN = 1, base = exp(1))
DivSpecCorFilt = diversity(speciesCF, index = "shannon", MARGIN = 1, base = exp(1))
DivSpecCorSamFilt = diversity(speciesCSF, index = "shannon", MARGIN = 1, base = exp(1)) 
boxplot(DivSpecUnfilt, DivSpecCorFilt, DivSpecCorSamFilt, ylab = "species_diversity", las = 2, names = c("prefilter", "filter 1", "filter 2"))

wGenU = grep("^G_", (colnames(met_unfilt_clades)))
GenU = met_unfilt_clades[,wGenU]
wGenCF  = grep("^G_", (colnames(met_corr_filt_clades)))
GenCF = met_corr_filt_clades[,wGenUCF]
wGenCSF  = grep("^G_", (colnames(met_corr_sample_filt_clades_abund)))
GenCSF = met_corr_sample_filt_clades_abund[,wGenCSF]
DivGenUnfilt = diversity(GenU, index = "shannon", MARGIN = 1, base = exp(1))
DivGenCorFilt = diversity(GenCF, index = "shannon", MARGIN = 1, base = exp(1))
DivGenCorSamFilt = diversity(GenCSF, index = "shannon", MARGIN = 1, base = exp(1)) 
boxplot(DivGenUnfilt, DivGenCorFilt, DivGenCorSamFilt, ylab = "genus_diversity", las = 2, names = c("prefilter", "filter 1", "filter 2"))

wFaU = grep("^F_", (colnames(met_unfilt_clades)))
FaU = met_unfilt_clades[,wFaU]
wFaCF  = grep("^F_", (colnames(met_corr_filt_clades)))
FaCF = met_corr_filt_clades[,wFaCF]
wFaCSF  = grep("^F_", (colnames(met_corr_sample_filt_clades_abund)))
FaCSF = met_corr_sample_filt_clades_abund[,wFaCSF]
DivFaUnfilt = diversity(FaU, index = "shannon", MARGIN = 1, base = exp(1))
DivFaCorFilt = diversity(FaCF, index = "shannon", MARGIN = 1, base = exp(1))
DivFaCorSamFilt = diversity(FaCSF, index = "shannon", MARGIN = 1, base = exp(1)) 
boxplot(DivFaUnfilt, DivFaCorFilt, DivFaCorSamFilt, ylab = "shannon_diversity_of_family", las = 2, names = c("unfiltered", "step 1", "step 2"))

wOrU = grep("^O_", (colnames(met_unfilt_clades)))
OrU = met_unfilt_clades[,wOrU]
wOrCF  = grep("^O_", (colnames(met_corr_filt_clades)))
OrCF = met_corr_filt_clades[,wOrCF]
wOrCSF  = grep("^O_", (colnames(met_corr_sample_filt_clades_abund)))
OrCSF = met_corr_sample_filt_clades_abund[,wOrCSF]
DivOrUnfilt = diversity(OrU, index = "shannon", MARGIN = 1, base = exp(1))
DivOrCorFilt = diversity(OrCF, index = "shannon", MARGIN = 1, base = exp(1))
DivOrCorSamFilt = diversity(OrCSF, index = "shannon", MARGIN = 1, base = exp(1)) 
boxplot(DivOrUnfilt, DivOrCorFilt, DivOrCorSamFilt, ylab = "shannon_diversity_of_orders", las = 2, names = c("unfiltered", "step 1", "step 2"))

wClU = grep("^C_", (colnames(met_unfilt_clades)))
ClU = met_unfilt_clades[,wClU]
wClCF  = grep("^C_", (colnames(met_corr_filt_clades)))
ClCF = met_corr_filt_clades[,wClCF]
wClCSF  = grep("^C_", (colnames(met_corr_sample_filt_clades_abund)))
ClCSF = met_corr_sample_filt_clades_abund[,wClCSF]
DivClUnfilt = diversity(ClU, index = "shannon", MARGIN = 1, base = exp(1))
DivClCorFilt = diversity(ClCF, index = "shannon", MARGIN = 1, base = exp(1))
DivClCorSamFilt = diversity(ClCSF, index = "shannon", MARGIN = 1, base = exp(1)) 
boxplot(DivClUnfilt, DivClCorFilt, DivClCorSamFilt, ylab = "shannon_diversity_of_class", las = 2, names = c("unfiltered", "step 1", "step 2"))

wPhU = grep("^P_", (colnames(met_unfilt_clades)))
PhU = met_unfilt_clades[,wPhU]
wPhCF  = grep("^P_", (colnames(met_corr_filt_clades)))
PhCF = met_corr_filt_clades[,wPhCF]
wPhCSF  = grep("^P_", (colnames(met_corr_sample_filt_clades_abund)))
PhCSF = met_corr_sample_filt_clades_abund[,wPhCSF]
DivPhUnfilt = diversity(PhU, index = "shannon", MARGIN = 1, base = exp(1))
DivPhCorFilt = diversity(PhCF, index = "shannon", MARGIN = 1, base = exp(1))
DivPhCorSamFilt = diversity(PhCSF, index = "shannon", MARGIN = 1, base = exp(1)) 
boxplot(DivPhUnfilt, DivPhCorFilt, DivPhCorSamFilt, ylab = "shannon_diversity_of_phylum", las = 2, names = c("unfiltered", "step 1", "step 2"))
