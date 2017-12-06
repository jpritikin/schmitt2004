##################################################################################
# TEMPLATE.R
# MASTER R FILE TO RUN GENETICALLY INFORMATIVE, DYNAMIC LATENT GROWTH CURVE MODELS
#
# J. Eric Schmitt and Michael C. Neale
# 
# From Schmitt, Neale, Fassassi, Perez, Lenroot, Wells, and Giedd
# The dynamic role of genetics on cortical patterning during childhood and adolescence
#
# Required Files:
# 1) TEMPLATE.R     (this file)
# 2) voxAC.R        (Open Mx Model Script)
# 3) concatdata.R   (place model parameters in a vector)
# 4) submodels.R    (submodels to test statistical significance of key parameters)
# 5) your data
# 6) R libraries OpenMx, plyr, reshape

voxseg<-1:5000 #SELECTS A SUBSAMPLE OF BRAIN VERTICES FOR PARALLEL COMPUTING

load("../../Giedd.RData") #LOAD DATA


rm(dt_r) # REMOVES CONTRALATERAL VERTEX-LEVEL DATA FROM WORKSPACE TO FREE MEMORY

#SELECT DEMOGRAPHIC VARIABLES
dem<-c("FAMILYID","GROUP","MRINUM","PERSONID","SEX","AGESCAN","ordernum","ZYG","TWIN")


###### REFORMATTING VARIABLES NAMES SO THAT MX LIKES THEM

selVarSuf1<-c("_1_1","_2_1","_3_1","_4_1","_5_1","_6_1","_7_1","_8_1")
selVarSuf2<-c("_1_2","_2_2","_3_2","_4_2","_5_2","_6_2","_7_2","_8_2")
selVarSuf3<-c("_1_3","_2_3","_3_3","_4_3","_5_3","_6_3","_7_3","_8_3")
selVarSuf4<-c("_1_4","_2_4","_3_4","_4_4","_5_4","_6_4","_7_4","_8_4")
selVarSuf5<-c("_1_5","_2_5","_3_5","_4_5","_5_5","_6_5","_7_5","_8_5")

ageVars1<-paste("data.AGESCAN",selVarSuf1,sep="")
ageVars2<-paste("data.AGESCAN",selVarSuf2,sep="")
ageVars3<-paste("data.AGESCAN",selVarSuf3,sep="")
ageVars4<-paste("data.AGESCAN",selVarSuf4,sep="")
ageVars5<-paste("data.AGESCAN",selVarSuf5,sep="")

nageVars1<-paste("AGESCAN",selVarSuf1,sep="")
nageVars2<-paste("AGESCAN",selVarSuf2,sep="")
nageVars3<-paste("AGESCAN",selVarSuf3,sep="")
nageVars4<-paste("AGESCAN",selVarSuf4,sep="")
nageVars5<-paste("AGESCAN",selVarSuf5,sep="")


#Concatenate age variables
ageVars<-c(nageVars1,nageVars2,nageVars3,nageVars4,nageVars5)


defVars<-c("ZYG_1_1","TWIN_1_1") # DEFINITIION VARIABLES (ZYGOSITY OF TWIN, YES/NO TWIN)


selVars1<-paste("CTvec",selVarSuf1,sep="") #Tell Mx which variables are observed in the model.                               
selVars2<-paste("CTvec",selVarSuf2,sep="")
selVars3<-paste("CTvec",selVarSuf3,sep="")
selVars4<-paste("CTvec",selVarSuf4,sep="")
selVars5<-paste("CTvec",selVarSuf5,sep="")


selVars<-c(selVars1,selVars2,selVars3,selVars4,selVars5) # SELECTED VARIABLES FOR MODEL

allVars<-c(selVars,ageVars,defVars)



## LOAD REQUIRED R LIBRARIES
library(plyr,lib.loc="~/R/lib/")
library(reshape,lib.loc="~/R/lib/")
library(OpenMx)


# BEGIN LOOP FOR EVERY VERTEX

for(i in voxseg){

#GENERATE TEMPORARY OUTPUT VECTORS FOR VERTEX i

outmatrix<-matrix(NA,nrow=1,ncol=48) # PARAMETER ESTIMATES
submodels<-matrix(NA,nrow=1,ncol=31) # SUBMODELS



#SUBSET DATA TO ONLY INCLUDE VARIABLES FOR THIS ANALYSIS

CTvec<-dt_l[i,]*10 #SELECT CORTICAL THICKNESS MEASURES FOR iTH VERTEX

modvars<-c(dem) #SELECT SUBSET VARIABLES TO PUSH TO OMX

subdat<-subset(gf_l, select=modvars) #SUBSET DEMOGRAPHIC/SCAN DATA 

a<-as.data.frame(cbind(subdat,CTvec)) # COMBINE SUBSETTED DEMOGRAPHIC/SCAN DATA WITH iTH VERTEX

rm(subdat) # FREE MEMORY

m1<-sapply(subset(a,select=CTvec),mean) #CALCULATE MEAN STARTING VALUE
V1<-sapply(subset(a,select=CTvec),var)  #CALCULATE ROUGH VARIANCE FOR STARTING VALUE

regP<-lm(a$CTvec~a$AGESCAN) # LINEAR REGRESSION TO PROVIDE REASONABLE STARTING VALUES
intS<-regP$coefficients[1]
slopeS<-regP$coefficients[2]



#RESHAPE DATA FROM LONG FORMAT TO WIDE FORMAT


#FIRST PLACE ALL DATA FROM A SINGLE INDIVIDUAL ON A SINGLE LINE
	a<-subset(a,MRINUM<9)
	a<-as.data.frame(sapply(a[,],as.numeric))
	a<-a[order(a$MRINUM),] #sort because reshape treats timevar as a factor

	a$AGESCAN<-a$AGESCAN-3.27 #RECENTER AGE TO EARLIEST SCAN AS TIME ZERO


	ind<-reshape(a, idvar=c("PERSONID","FAMILYID","ordernum"),timevar="MRINUM", direction="wide") # INDIVIDUAL-WISE RECORDS

#THEN CONVERT RECORDS TO BE FAMILY-WISE

    ind<-ind[order(ind$ordernum),]
    fam<-reshape(ind,idvar="FAMILYID",timevar="ordernum",direction="wide")


rm(ind,a) # FREE MEMORY


#RENAME VARIABLES SO THAT OMX LIKES THEM

	names(fam)<-gsub("[.]","_",names(fam))


##NECESSARY EVIL SINCE DEFINITION VARIABLES CAN'T  BE NA

fam$AGESCAN_1_1[is.na(fam$AGESCAN_1_1)]<- 0
fam$AGESCAN_2_1[is.na(fam$AGESCAN_2_1)]<- 0
fam$AGESCAN_3_1[is.na(fam$AGESCAN_3_1)]<- 0
fam$AGESCAN_4_1[is.na(fam$AGESCAN_4_1)]<- 0
fam$AGESCAN_5_1[is.na(fam$AGESCAN_5_1)]<- 0
fam$AGESCAN_6_1[is.na(fam$AGESCAN_6_1)]<- 0
fam$AGESCAN_7_1[is.na(fam$AGESCAN_7_1)]<- 0
fam$AGESCAN_8_1[is.na(fam$AGESCAN_8_1)]<- 0

fam$AGESCAN_1_2[is.na(fam$AGESCAN_1_2)]<- 0
fam$AGESCAN_2_2[is.na(fam$AGESCAN_2_2)]<- 0
fam$AGESCAN_3_2[is.na(fam$AGESCAN_3_2)]<- 0
fam$AGESCAN_4_2[is.na(fam$AGESCAN_4_2)]<- 0
fam$AGESCAN_5_2[is.na(fam$AGESCAN_5_2)]<- 0
fam$AGESCAN_6_2[is.na(fam$AGESCAN_6_2)]<- 0
fam$AGESCAN_7_2[is.na(fam$AGESCAN_7_2)]<- 0
fam$AGESCAN_8_2[is.na(fam$AGESCAN_8_2)]<- 0

fam$AGESCAN_1_3[is.na(fam$AGESCAN_1_3)]<- 0
fam$AGESCAN_2_3[is.na(fam$AGESCAN_2_3)]<- 0
fam$AGESCAN_3_3[is.na(fam$AGESCAN_3_3)]<- 0
fam$AGESCAN_4_3[is.na(fam$AGESCAN_4_3)]<- 0
fam$AGESCAN_5_3[is.na(fam$AGESCAN_5_3)]<- 0
fam$AGESCAN_6_3[is.na(fam$AGESCAN_6_3)]<- 0
fam$AGESCAN_7_3[is.na(fam$AGESCAN_7_3)]<- 0
fam$AGESCAN_8_3[is.na(fam$AGESCAN_8_3)]<- 0

fam$AGESCAN_1_4[is.na(fam$AGESCAN_1_4)]<- 0
fam$AGESCAN_2_4[is.na(fam$AGESCAN_2_4)]<- 0
fam$AGESCAN_3_4[is.na(fam$AGESCAN_3_4)]<- 0
fam$AGESCAN_4_4[is.na(fam$AGESCAN_4_4)]<- 0
fam$AGESCAN_5_4[is.na(fam$AGESCAN_5_4)]<- 0
fam$AGESCAN_6_4[is.na(fam$AGESCAN_6_4)]<- 0
fam$AGESCAN_7_4[is.na(fam$AGESCAN_7_4)]<- 0
fam$AGESCAN_8_4[is.na(fam$AGESCAN_8_4)]<- 0

fam$AGESCAN_1_5[is.na(fam$AGESCAN_1_5)]<- 0
fam$AGESCAN_2_5[is.na(fam$AGESCAN_2_5)]<- 0
fam$AGESCAN_3_5[is.na(fam$AGESCAN_3_5)]<- 0
fam$AGESCAN_4_5[is.na(fam$AGESCAN_4_5)]<- 0
fam$AGESCAN_5_5[is.na(fam$AGESCAN_5_5)]<- 0
fam$AGESCAN_6_5[is.na(fam$AGESCAN_6_5)]<- 0
fam$AGESCAN_7_5[is.na(fam$AGESCAN_7_5)]<- 0
fam$AGESCAN_8_5[is.na(fam$AGESCAN_8_5)]<- 0

fam$ZYG_1_1[is.na(fam$ZYG_1_1)]<- 0
fam$TWIN_1_1[is.na(fam$TWIN_1_1)]<- 0


fam$SEX_1_1[is.na(fam$SEX_1_1)]<- 0
fam$SEX_1_2[is.na(fam$SEX_1_2)]<- 0
fam$SEX_1_3[is.na(fam$SEX_1_3)]<- 0
fam$SEX_1_4[is.na(fam$SEX_1_4)]<- 0
fam$SEX_1_5[is.na(fam$SEX_1_5)]<- 0


source("voxAC.R") # RUN MX MODEL

source("concatdata.R") # GENERATE OUTPUT VECTOR

outmatrix<-ROIout #ASSIGN OUTPUT VECTOR TO MATRIX (REDUNDANT FOR ERROR CHECKING)

source("submodels.R") # RUN SUBMODELS

submodels<-submodelstats
names(submodels)[1]<-"i2"


f<-as.data.frame(cbind(t(outmatrix),t(submodels))) # COMBINE RESULTS FROM INITIAL AND SUBMODELS

write.table(f,file="LGCout.txt",col.names=F,row.names=F,na=".",sep=" ",append=TRUE) #RIGHT COMBINED OUTPUT TO TXT FILE

# LOOP TO i+1 VERTEX

} 

##ASSIGN COLUMN NAMES TO OUTPUT FILES

parnam<-c("a11","a21","a31","a22","a32","a33",
"c11","c21","c31","c22","c32","c33",
"e11","e21","e31","e22","e32","e33",
"t11","t21","t31","t22","t32","t33",
"v11","v21","v31","v22","v32","v33")

vcnam <-c("A11","A21","A31","A22","A32","A33",
"C11","C21","C31","C22","C32","C33",
"E11","E21","E31","E22","E32","E33",
"T11","T21","T31","T22","T32","T33",
"V11","V21","V31","V22","V32","V33")

stdnam<-c("a2_11","a2_21","a2_31","a2_22","a2_32","a2_33",
"c2_11","c2_21","c2_31","c2_22","c2_32","c2_33",
"e2_11","e2_21","e2_31","e2_22","e2_32","e2_33",
"t2_11","t2_21","t2_31","t2_22","t2_32","t2_33")

# SUBMODEL PARAMETERS (MANY NOT UTILIZED IN THE PRESENT STUDY)

nsubmodel_no <-c("LLnoA","X2noA","pNoA","LLnoE","X2noE","pNoE")
nsubmodel_sq <-c("LLnoAsq","X2noAsq","pNoAsq","LLnoEsq","X2noEsq","pNoEsq")
nsubmodel_del<-c("LLnoAdel","X2noAdel","pNoAdel","LLnoEdel","X2noEdel","pNoEdel")
nsubmodel_cor<-c("LLnoAcor","X2noAcor","pNoAcor","LLnoEcor","X2noEcor","pNoEcor")
nsubmodel_tot<-c("LLnosquare","X2nosquare","pNosquare","LLnoslope","X2noslope","pNoslope")
statnam<-c("i2",nsubmodel_no,nsubmodel_sq,nsubmodel_del,nsubmodel_cor,nsubmodel_tot)

outnam<-c("i","ROI",parnam,vcnam,stdnam,"int","age","age2","sex","err","LL")

write.table(outnam,file="outnam.txt",row.names=F,col.names=F)

