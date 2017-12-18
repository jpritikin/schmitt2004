############################################################################################
# voxAC.R
# GENETICALLY INFORMATIVE LATENT GROWTH CURVE MODEL
# J. Eric Schmitt and Michael C. Neale
#
# From Schmitt, Neale, Fassassi, Perez, Lenroot, Wells, and Giedd
# 
# The dynamic role of genetics on cortical patterning during childhood and adolescence
# 
#
#


require(OpenMx) #Load Open Mx 


# Define dimensions of parameter space
nl<-3       #Number of latent common factors (Intercept, Linear Slope, Quadratic Slope)
nvars<-8    #Maximum number of observed variables per individual
nsibs<-5    #Maximum number of siblings (including MZ and DZ twins)
allvars<-40 #Total number of observed variables (nvars*nsibs)

#STARTING VALUES
startVC<-V1**.5
startM<-m1
slopeS2<-0
sexS<-0


####### DEFINE LATENT FACTOR STRUCTURE FOR VARIANCE DECOMPOSITION

# Genetic factor (Lower triangular matrix)/Cholesky Decomposition

lgcACE<- mxModel("lgcACE",

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=TRUE,
		values=c(startVC/10,0,0,0,0,0),
		labels=c("a11","a21","a31","a22","a32","a33"),
		name="Al"
	),

# Shared Environmental Factor Option (Fixed to 0 in the current analysis) (Lower Triangular/Cholesky)
# Note that freeing these parameters will significantly increase computation times.
	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=FALSE,
		values=c(0,0,0,0,0,0),
		labels=c("c11","c21","c31","c22","c32","c33"),
		name="Cl"
	),


#Unique Environmental COMMON factor (Lower Triangular/Cholesky)
	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=TRUE,
		values=c(startVC/10,0,0,0,0,0),
		labels=c("e11","e21","e31","e22","e32","e33"),
		name="El"
	),


#Twin Specific Environmental Factor Option (Fixed to 0)
# Note that freeing these parameters will significantly increase computation times.

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=FALSE,
		values=c(0,0,0,0,0,0),
		labels=c("t11","t21","t31","t22","t32","t33"),
		name="Tl"
	),
												

### AGE Definition Variables 

#EACH VARIABLE IS A VECTOR OF AGE AT TIME OF SCAN FOR THE ith INDIVIDUAL OF A FAMILY

	mxMatrix(
		type="Full",
		nrow=nvars,
		ncol=1,
		free=FALSE,
		labels=ageVars1,
		name = "AGE1"
	),

	mxMatrix(
		type="Full",
		nrow=nvars,
		ncol=1,
		free=FALSE,
		labels=ageVars2,
		name = "AGE2"
	),

	mxMatrix(
		type="Full",
		nrow=nvars,
		ncol=1,
		free=FALSE,
		labels=ageVars3,
		name = "AGE3"
	),

	mxMatrix(
		type="Full",
		nrow=nvars,
		ncol=1,
		free=FALSE,
		labels=ageVars4,
		name = "AGE4"
	),       	

	mxMatrix(
		type="Full",
		nrow=nvars,
		ncol=1,
		free=FALSE,
		labels=ageVars5,
		name = "AGE5"
	),



#ZYGOSITY DEFINITION VARIABLES FOR TWIN PAIRS {1 MZ, 0.5 DZ or no Twins in family}

	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.ZYG_1_1"),
		name = "ZYG"
	),


#TWIN DEFINITION VARIABLES {1 if twin, 0 if not twin}; NOT USED IN VERTEX-LEVEL MODELING
#FOR MODELING TWIN SPECIFIC ENVIRONMENT 

	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.TWIN_1_1"),
		name = "TWIN"
	),


#SEX DEFINITION VARIABLES FOR THE ith INDIVIDUAL IN EACH FAMILY (up to 5)

	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.SEX_1_1"),
		name = "SEX1"
	),


	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.SEX_1_2"),
		name = "SEX2"
	),

	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.SEX_1_3"),
		name = "SEX3"
	),


	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.SEX_1_4"),
		name = "SEX4"
	),


	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=FALSE,
		labels=c("data.SEX_1_5"),
		name = "SEX5"
	),




####           GENERATE ACET VARIANCE COMPONENTS (C AND T FIXED TO 0)


### ADDITIVE GENETIC VARIANCE
mxAlgebra(
	expression=Al %*% t(Al),
	name="A"
),

### SHARED ENVIRONMENTAL VARIANCE OPTION
mxAlgebra(
	expression=Cl %*% t(Cl),
	name="C"
),

### UNIQUE ENVIRONMENTAL VARIANCE
mxAlgebra(
	expression=El %*% t(El),
	name="E"
),

### TWIN SPECIFIC ENVIRONMENTAL VARIANCE OPTION
mxAlgebra(
	expression=Tl %*% t(Tl),
	name="T"
),
	


#### EXPECTED LATENT COVARIANCE MATRIX


mxAlgebra(
	expression= rbind(
		cbind(A+C+T+E,  		ZYG%x%A+C+TWIN%x%T,	0.5%x%A+C,	0.5%x%A+C,	0.5%x%A+C),
		cbind(ZYG%x%A+C+TWIN%x%T,	A+C+T+E,		0.5%x%A+C,	0.5%x%A+C,	0.5%x%A+C),
		cbind(0.5%x%A+C,		0.5%x%A+C,		A+C+T+E,	0.5%x%A+C,	0.5%x%A+C),
		cbind(0.5%x%A+C,		0.5%x%A+C,		0.5%x%A+C,	A+C+T+E,	0.5%x%A+C),
		cbind(0.5%x%A+C,		0.5%x%A+C,		0.5%x%A+C,	0.5%x%A+C,	A+C+T+E)
	),
	name="GenCov"
),




### BUILDING BLOCKS FOR LGC COMPONENT OF THE MODEL


### IDENTITY MATRIX
	mxMatrix(
		type="Iden",
		nrow=allvars,
		ncol=allvars,
		free=FALSE,
		name = "Ivars"
	),

### ZERO MATRIX
	mxMatrix(
		type="Zero",
		nrow=nvars,
		ncol=nl,
		free=FALSE,
		name="zero",
	),

### UNIT MATRIX
	mxMatrix(
		type="Unit",
		nrow=nvars,
		ncol=1,
		free=FALSE,
		name="IntLoad"
	),


### MEANS VECTOR (regressing CT on age, age^2, age^3, and gender)
	mxMatrix(
		type="Full",
		nrow=4,
		ncol=1,
		free=TRUE,
		value=c(intS,slopeS,slopeS2,sexS),
		labels=c("bint","bslope","bslope2","sex"),
		name="betas"
	),

### ERROR/RESIDUAL LATENT FACTOR. 
	
	mxMatrix(
		type="Full",
		nrow=1,
		ncol=1,
		free=TRUE,
		value=startVC/2,
		labels=("error"),
		name="err",
	),



### CREATE DYNAMIC LATENT GROWTH CURVE MODEL WITH FACTOR LOADINGS BASED ON AGE AT TIME OF SCAN DEFINITION VARIABLES

mxAlgebra(
	expression=rbind(
			cbind(IntLoad,AGE1,AGE1*AGE1,zero,zero,zero,zero),
			cbind(zero, IntLoad,AGE2,AGE2*AGE2,zero,zero,zero),
			cbind(zero,zero,IntLoad,AGE3,AGE3*AGE3,zero,zero),
			cbind(zero,zero,zero,IntLoad,AGE4,AGE4*AGE4,zero),
			cbind(zero,zero,zero,zero,IntLoad,AGE5,AGE5*AGE5)
			),
			name="LGC"
),		



### COMBINE VARIANCE DECOMPOSITION WITH LGC TO GENERATE AN EXPECTED COVARIANCE MATRIX

mxAlgebra( 
	expression=((LGC %*% GenCov %*% t(LGC))+Ivars%x%(err**2)),
	name="CovMat"
	),





### MEANS VECTOR (for regressing CT on age, age^2, age^3, and gender)
        mxMatrix(
                type="Full",
                nrow=4,
                ncol=1,
                free=TRUE,
                value=c(intS,slopeS,slopeS2,sexS),
                labels=c("bint","bslope","bslope2","sex"),
                name="betas"
        ),


### GENERATE EXPECTED MEANS MATRIX

mxAlgebra(
	expression=cbind(
			t(cbind(IntLoad,AGE1,AGE1*AGE1,SEX1%x%IntLoad) %*% betas),
			t(cbind(IntLoad,AGE2,AGE2*AGE2,SEX2%x%IntLoad) %*% betas),
			t(cbind(IntLoad,AGE3,AGE3*AGE3,SEX3%x%IntLoad) %*% betas),
			t(cbind(IntLoad,AGE4,AGE4*AGE4,SEX4%x%IntLoad) %*% betas),
			t(cbind(IntLoad,AGE5,AGE5*AGE5,SEX5%x%IntLoad) %*% betas)),
			name="MeansMat"
),




### LINK OBSERVED DATA WITH SCRIPT AND DEFINE COVARIANCE AND MEANS MATRICES 

	mxData(
		observed=fam,
		type="raw",
	),
	mxFIMLObjective(

		covariance="CovMat",
		means="MeansMat",
		dimnames=selVars
	)


)

### RUN THE MODEL

fit<-mxRun(lgcACE)


