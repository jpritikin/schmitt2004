#############################################################
#
# submodels.R
#
# HYPOTHESIS TESTING LATENT GROWTH CURVE MODEL COMPONENTS
#
# J. Eric Schmitt and Michael C. Neale
#
# From Schmitt, Neale, Fassassi, Perez, Lenroot, Wells, and Giedd
#
# The dynamic role of genetics on cortical patterning during childhood and adolescence. 
#


##########################################
# 
# SIGNIFICANCE OF ALL GENETIC EFFECTS (environmental factors remain free)

lgcNoA<-mxRename(lgcACE,"noA")

lgcNoA$Al<-

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=FALSE,
		values=c(0,0,0,0,0,0),
		labels=c("a11","a21","a31","a22","a32","a33"),
		name="Al"
	)

FnoA<-mxRun(lgcNoA)

LLnoA<-mxEval(objective,FnoA)
X2noA<-LLnoA-LL
pNoA<-pchisq(X2noA,6,lower.tail=F)

##########################################
#
# SIGNIFICANCE OF ENVIRONMENTAL COMMON EFFECTS (genetic factors and error remain free)

lgcNoE<-mxRename(lgcACE,"noE")

lgcNoE$El<-

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=FALSE,
		values=c(0,0,0,0,0,0),
		labels=c("e11","e21","e31","e22","e32","e33"),
		name="El"
	)

FnoE<-mxRun(lgcNoE)

LLnoE<-mxEval(objective,FnoE)
X2noE<-LLnoE-LL
pNoE<-pchisq(X2noE,6,lower.tail=F)

###############################################


##########################################
#
# SIGNIFICANCE OF GENETIC CHANGE (allowing for a main effect)

lgcNoAdel<-mxRename(lgcACE,"noAdel")

lgcNoAdel$Al<-

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
		values=c(startVC/40,0,0,0,0,0),
		labels=c("a11","a21","a31","a22","a32","a33"),
		name="Al"
	)

FnoAdel<-mxRun(lgcNoAdel)

LLnoAdel<-mxEval(objective,FnoAdel)
X2noAdel<-LLnoAdel-LL
pNoAdel<-pchisq(X2noAdel,5,lower.tail=F)


##########################################
#
# SIGNIFICANCE OF ENVIRONMENTAL CHANGE (allowing for a main effect)

lgcNoEdel<-mxRename(lgcACE,"noEdel")

lgcNoEdel$El<-

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
		values=c(startVC/40,0,0,0,0,0),
		labels=c("e11","e21","e31","e22","e32","e33"),
		name="El"
	)

FnoEdel<-mxRun(lgcNoEdel)

LLnoEdel<-mxEval(objective,FnoEdel)
X2noEdel<-LLnoEdel-LL
pNoEdel<-pchisq(X2noEdel,5,lower.tail=F)

###############################################


##########################################
#
#
# SIGNIFICNANCE OF ALL CHANGE (BOTH GENETIC AND ENVIRONMENTAL)
# ALLOWS FOR GENETIC AND NONGENETIC MAIN EFFECTS AND ERROR

lgcNoslope<-mxRename(lgcACE,"noSLOPE")

lgcNoslope$Al<-

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
		values=c(startVC/40,0,0,0,0,0),
		labels=c("a11","a21","a31","a22","a32","a33"),
		name="Al"
	)

lgcNoslope$El<-

	mxMatrix(
		type="Lower",
		nrow=nl,
		ncol=nl,
		free=c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
		values=c(startVC/40,0,0,0,0,0),
		labels=c("e11","e21","e31","e22","e32","e33"),
		name="El"
	)

Fnoslope<-mxRun(lgcNoslope)

LLnoslope<-mxEval(objective,Fnoslope)
X2noslope<-LLnoslope-LL
pNoslope<-pchisq(X2noslope,10,lower.tail=F)

###############################################



### CONCATENATE STATISTICS TO FACILITATE MANAGEMENT

submodel_no <-c(LLnoA,X2noA,pNoA,LLnoE,X2noE,pNoE)
submodel_del<-c(LLnoAdel,X2noAdel,pNoAdel,LLnoEdel,X2noEdel,pNoEdel)
submodel_tot<-c(LLnoslope,X2noslope,pNoslope)

#LINK SUMMARY STATS TO VERTEX NUMBER TO ENSURE PROPER SYNCHRONIZATION OF OUTPUT

submodelstats<-c(i,submodel_no,submodel_del,submodel_tot)



