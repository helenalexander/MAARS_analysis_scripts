rm(list=ls())

#load segre abundance data
segre_data <- read.csv("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/contamination/segre/20160308segredata.csv",
               header = TRUE,
               skip = 3,
               quote="\"",
               stringsAsFactors= TRUE,
               strip.white = TRUE)
head(segre_data)
dim(segre_data)
#make clade names rownames
rownames(segre_data) = segre_data[,1]
head(segre_data)
segre_data = segre_data[,-1]
head(segre_data)

# make abundances numeric
segre_data_num <- apply(segre_data,2,as.numeric) 
head(segre_data_num)
rownames(segre_data_num) = rownames(segre_data)
head(segre_data_num)
dim(segre_data_num)

# make clade names only include species name
cc= strsplit(x=(rownames(segre_data_num)), split = ";") 
auxFun = function(x){ l=length(x); return(x[l])}
Segspecies = unlist(lapply(FUN=auxFun, X = cc))
rownames(segre_data_num) = Segspecies
head(segre_data_num)
dim(segre_data_num)


###########################################
# 1. Compare which species have the top mean abundance in mocks and controls with Segre data
############################################
# Calculate mean abundance of each species
s_data_means = rowMeans(segre_data_num)
# Order the mean abundances so most abundant is first
SS = sort(-s_data_means)
length(SS)
# the top 20 most abundant species:
S20 = SS[1:20]
# the top 50:
S50 = SS[1:50]
# the top 10
S10 = SS[1:10]
length(S20)
names(S50)
length(S10)

# load reduced abundance matrix
load ("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/data-filtering/20160304met_corr_sample_filt_clades_abund.Rdata")

head(met_corr_sample_filt_clades_abund)
dim(met_corr_sample_filt_clades_abund)

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
# subset the species only
wspec = grep("^S", rownames(CTRLAbt))
SpCA = CTRLAbt[wspec,]
dim(SpCA)
head(SpCA)
# mean abundances of species in control samples
MeansSpCA = rowMeans(SpCA)
length(MeansSpCA)
# order these with most abundant first
SC = sort(-MeansSpCA)
length(SC)
CTRL20 = SC[1:20]
CTRL50 = SC[1:50]
CTRL10 = SC[1:10]


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

wspecm = grep("^S", rownames(MockAbt))
SpM = MockAbt[wspecm,]
head(SpM)
MeansSpM = rowMeans(SpM)
length(MeansSpM)

SM = sort(-MeansSpM)
length(SM)
Mock20 = SM[1:20]
Mock50 = SM[1:50]
Mock10 = SM[1:10]

names(S10) = paste ("S", names(S10), sep="_")
names(S20) = paste ("S", names(S20), sep="_")
names(S50) = paste ("S", names(S50), sep="_")
names(SS) = paste ("S", names(SS), sep="_")

# intersect between segre top 10 species and control top 10 species
SegCTRL10 = intersect(names(S10), names(CTRL10))
length(SegCTRL10)

SegMock10 = intersect (names(S10), names(Mock10))
length(SegMock10)

# overalap in top 20 species in Segre's data and our controls and mocks
SegCTRL20 = intersect (names(S20), names(CTRL20))
length(SegCTRL20)

SegMock20 = intersect (names(S20), names(Mock20))
length(SegMock20)

# overalap in top 50 species in Segre's data and our controls and mocks
SegCTRL50 = intersect (names(S50), names(CTRL50))
length(SegCTRL50)

SegMock50 = intersect (names(S50), names(Mock50))
length(SegMock50)

# overalap in all species in Segre's data and our controls and mocks
SegCTRLAll = intersect (names(SS), names(SC))
length(SegCTRLAll)

SegMockAll = intersect (names(SS), names(SM))
length(SegMockAll)

###########################################
# 2. Compare the intersection with segres samples with the abundance of all of the top 20 species in each sample and top 20 species in each mock
############################################
#Segre top 20 species
S20

# intersect for top 20 species in each control sample with top 20 species in segres data
dim(SpCA)
head(SpCA)

LIC = list()
for (i in 1:dim(SpCA)[2]){
  top20ab = sort(-(SpCA[,i]))[1:20]
  int = intersect(names(S20), names(top20ab))
  LengthIntersect = length(int)
  LIC[[(i)]] = LengthIntersect
  LengthIntersectCTRL = as.numeric(LIC)
}
LengthIntersectCTRL

# intersect for top 20 species in each mock sample with top 20 species in segres data
LIM = list()
for (i in 1:dim(SpM)[2]){
  top20ab = sort(-(SpM[,i]))[1:20]
  int = intersect(names(S20), names(top20ab))
  LengthIntersect = length(int)
  LIM[[(i)]] = LengthIntersect
  LengthIntersectMOCK = as.numeric(LIM)
}
LengthIntersectMOCK

# plot the disctributions
boxplot(LengthIntersectCTRL, LengthIntersectMOCK, names=c("Control_samples", "Mock_samples"), ylab = "length_of_intersect_with_skin_reference")




