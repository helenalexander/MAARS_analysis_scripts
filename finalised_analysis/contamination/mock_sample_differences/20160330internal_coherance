rm(list=ls())

# load reduced abundance matrix
load ("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/data-filtering/20160304met_corr_sample_filt_clades_abund.Rdata")

head(met_corr_sample_filt_clades_abund)
dim(met_corr_sample_filt_clades_abund)

## mocks
# load metadata 
load ("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/20160301clin_seq_pt_metadata.Rdata")
colnames(clin_seq_pt_metadata)
clin_seq_pt_metadata[,"MOCK_ID"]

wmock = !is.na(clin_seq_pt_metadata[,"MOCK_ID"])
mocks = clin_seq_pt_metadata[wmock,"scilifeID"]
length(mocks)

MockSamples = c("P1411_1037", "P1411_1062" ,"P1411_1086", "P1411_1108" ,"P1709_1037", "P1761_1012", "P1761_1037", "P1761_1062" ,
                "P1761_1085", "P1761_1107" ,"P1821_1048", "P1821_1108" ,"P1896_1047", "P1896_1083", "P1896_1106",
                "P1896_1107", "P1900_1047", "P2029_1078", "P2029_1079", "P752_119",   "P752_120",   "P752_121"  ,
                "P752_122"  , "P752_123",   "P752_124"  )
MockAbund = met_corr_sample_filt_clades_abund[MockSamples,]
dim(MockAbund)
MockAbt = t(MockAbund)
dim(MockAbt)

corMocks = cor(MockAbt)
diag (corMocks) = 0
hist(corMocks, 30)
boxplot(corMocks)

#plot correlation matrix
library(corrplot)
dir.create("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/contamination/contamination_plots")
pdf("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/contamination/contamination_plots/corrplot_mocks.pdf", width = 20, height = 20)
corrplot(corMocks, method="circle", tl.pos="lt", type="upper",
         tl.col="black", tl.cex=0.9, tl.srt=45)
dev.off()


# CTRLs
CTRLSamples = c( "P1411_1001", "P1411_1002", "P1411_1007", "P1411_1008", "P1411_1013", "P1411_1014", "P1411_1019", "P1411_1020", "P1411_1023", "P1411_1024", "P1411_1025", "P1411_1026",
                 "P1411_1031", "P1411_1032", "P1411_1035", "P1411_1036", "P1411_1038", "P1411_1039", "P1411_1044", "P1411_1045", "P1411_1048", "P1411_1049", "P1411_1050", "P1411_1051",
                 "P1411_1056", "P1411_1057", "P1411_1060", "P1411_1061", "P1411_1069", "P1411_1070", "P1411_1075", "P1411_1076", "P1411_1081", "P1411_1082", "P1411_1085", "P1411_1087",
                 "P1411_1098", "P1411_1099", "P1411_1106", "P1411_1107", "P1709_1007", "P1709_1008", "P1709_1085", "P1709_1086", "P1709_1087", "P1709_1088", "P1709_1089", "P1709_1090",
                 "P1709_1091", "P1709_1092", "P1709_1093", "P1709_1094", "P1709_1095", "P1709_1096", "P1709_1097", "P1709_1098", "P1709_1099", "P1709_1100", "P1709_1101", "P1709_1102",
                 "P1709_1103", "P1709_1104", "P1709_1105", "P1761_1001", "P1761_1002", "P1761_1007", "P1761_1008", "P1761_1025", "P1761_1026", "P1761_1031", "P1761_1032", "P1761_1035",
                 "P1761_1036", "P1761_1038", "P1761_1039", "P1761_1044", "P1761_1045", "P1761_1050", "P1761_1051", "P1761_1056", "P1761_1057", "P1761_1060", "P1761_1061", "P1761_1069",
                 "P1761_1070", "P1761_1075", "P1761_1076", "P1761_1081", "P1761_1084", "P1761_1086", "P1761_1097", "P1761_1098", "P1761_1105", "P1761_1106", "P1821_1017","P1821_1073", 
                 "P1821_1074","P1821_1075", "P1821_1078", "P1821_1079", "P1821_1080", "P1821_1081", "P1821_1082", "P1821_1083", "P1821_1084", "P1821_1085", "P1821_1086",
                 "P1821_1087", "P1821_1088", "P1821_1089", "P1821_1090", "P1821_1091", "P1821_1092", "P1821_1097", "P1821_1098", "P1821_1103", "P1821_1104", "P1896_1030", "P1896_1031",
                 "P1896_1032", "P1896_1040", "P1896_1042", "P1896_1044", "P1896_1053", "P1896_1054", "P1896_1055", "P1896_1056", "P1896_1057",
                 "P1896_1058", "P1896_1059", "P1896_1060", "P1896_1062", "P1896_1063", "P1896_1064", "P1896_1065", "P1896_1067", "P1896_1068", "P1896_1069", "P1896_1072", "P1896_1073",
                 "P1896_1074", "P1896_1076", "P1896_1077", "P1896_1078", "P1896_1079", "P1896_1080", "P1896_1081", "P1896_1082", "P1896_1084", "P1896_1085", "P1896_1090", "P1896_1091",
                 "P1896_1096", "P1896_1097", "P1896_1102", "P1896_1103", "P1900_1001", "P1900_1003", "P1900_1004", "P1900_1005", "P1900_1010", "P1900_1013", "P1900_1014", "P1900_1015",
                 "P1900_1016", "P1900_1017", "P1900_1018", "P1900_1019", "P1900_1020", "P1900_1021", "P1900_1022", "P1900_1024", "P1900_1060", "P1900_1079", "P1900_1080", "P1900_1082",
                 "P1900_1085", "P1900_1086", "P1900_1087", "P1900_1088", "P1900_1089", "P1900_1090", "P1900_1091", "P1900_1092", "P1900_1093", "P1900_1096", "P1900_1097", "P1900_1102",
                 "P1900_1103", "P2029_1018", "P2029_1019", "P2029_1020", "P2029_1022", "P2029_1024", "P2029_1026", "P2029_1029", "P2029_1031", "P2029_1033", "P2029_1060", "P2029_1063",
                 "P2029_1076", "P2029_1077", "P752_101"  , "P752_102" ,  "P752_107",   "P752_108" ,  "P752_113",   "P752_114" )
#Abundances of our CTRL samples 
CTRLAbund = met_corr_sample_filt_clades_abund[CTRLSamples,]
dim(CTRLAbund)
CTRLAbt = t(CTRLAbund)
dim(CTRLAbt)

corCTRL = cor(CTRLAbt)
diag (corCTRL) = 0
hist(corCTRL, 30)
head(corCTRL)
boxplot(corCTRL)

#plot correlation matrix
pdf("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/contamination/contamination_plots/corrplot_CTRL.pdf", width = 20, height = 20)
corrplot(corCTRL, method="circle", tl.pos="lt", type="upper",
         tl.col="black", tl.cex=0.5, tl.srt=45)
dev.off()

### calculate mean correlation for each sample/mock with all other samples/mocks
MeanCorMocks = rowMeans(corMocks)
head(MeanCorMocks)
length(MeanCorMocks)

MeanCorCTRL = rowMeans(corCTRL)
head(MeanCorCTRL)
length(MeanCorCTRL)


# calculate correlation of samples with mocks
dim(CTRLAbt)
dim(MockAbt)
colnames(MockAbt) = paste ("M", colnames(MockAbt), sep="_")
colnames(CTRLAbt) = paste ("C", colnames(CTRLAbt), sep="_")

CTRLMockAbt = cbind(CTRLAbt, MockAbt)
dim(CTRLMockAbt)
head(CTRLMockAbt)

corCTRLMock = cor(CTRLMockAbt)
diag (corCTRLMock) = 0
hist(corCTRLMock, 30)
head(corCTRLMock)
colnames(corCTRLMock)

wMC = grep("^M", colnames(corCTRLMock))
MC = corCTRLMock[,wMC]
dim(corCTRLMock)
head(MC)
dim(MC)
wCM = grep("^C", rownames(MC))
CM = MC[wCM,]
head(CM)
dim(CM)

MeanCorCTRLMocks = rowMeans(CM)
head(MeanCorCTRLMocks)
length(MeanCorCTRLMocks)

boxplot(MeanCorCTRL, MeanCorMocks, MeanCorCTRLMocks, names = c("CTRL:CTRL", "Mock:Mock", "CTRL:Mock"), ylab = "correlation")
