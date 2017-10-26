#-----------------------------------------------------------------------#
# Program: Bivariate Saturated and ACE Model				                    #
#  Author: Kelly Spoon                                                  #
#    Date: 07 18 2011                                                   #
#Modified: 09 17 2015                                                   #
#                                                                       #
# Description: Complete bivariate analysis (Saturated, ACE model, 	    #
#		   & Nested Models							                                    #
#     	   							                                                #
#       Input: Raw data file (CSV) with the following columns	 	        #
#  		   IN ORDER: subid, case, twin, zygosity, phenotype 1,2,3...      #
#		   names of columns do not matter					                          #
#		   *twin variable should be defined A or B			                    #
#		   *zygosity should be defined as 1 for MZ and 2 for DZ	            #
#												                                                #
#	   Data: Continuous								                                    #
#												                                                #
#      Output: All bivariate statistics                                 #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user       #
#-----------------------------------------------------------------------#

## Loading Required/Useful Libraries ##
require(OpenMx)		# Loads OpemMx
require(psych)		# Loads Psych package
source("http://www.vipbg.vcu.edu/~vipbg/Tc24/GenEpiHelperFunctions.R") # Loads GenEpi functions from VCU

### Reading in Original Data ###
#*** Set working directory - location of data file and desired output
	setwd("C:/Users/mspanizzon/Desktop/Jeremy Mx Training/Bivariate ACE Model")

#*** Read in data file
	twins <- read.csv("V1_HeightAFQT.csv",header=T)

#*** Quick check that data file is properly read in
	dim(twins)[1] 		# should be number of subjects
	dim(twins)[2] 		# should be number of variables in data file
	names(twins)[1:4] 	# should be subid, case, twin, zygosity
	names(twins)[7:8] 	# should be the phenotypes of interest


#*** Defining number of variables
nv <- 2				#Number of phenotypes (NVAR)
ntv <- nv*2

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg10"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
	names(newtwins) 		# Gives the variable names for the new paired dataset

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg10==1,c(7,8,13,14)))
	DZdata <- as.data.frame(subset(newtwins,zyg10==2,c(7,8,13,14)))

selvars <- names(MZdata)

# Print Descriptive Statistics
describe(MZdata)
describe(DZdata)

cor(MZdata,use="complete")
cor(DZdata,use="complete")

# Getting Number of Twins Used in Analysis
	cpMZ <- sum(complete.cases(MZdata)) # Number of complete MZ pairs
	sMZ <- dim(MZdata)[1]-cpMZ 		# Number of MZ singletons
	cpDZ <- sum(complete.cases(DZdata)) # Number of complete DZ pairs
	sDZ <- dim(DZdata)[1]-cpDZ 		# Number of DZ singletons

#---------------------------------#
## Fit Bivariate Saturated Model ##
#---------------------------------#

meanSV <- c(0,0,0,0)  				# Start values for the means
cholSV <- c(.9,0,0,0,.9,0,0,.9,0,.9) 	# Start values for the lower var/cov matrix 

# The following lines are optional, but are useful for labeling elements of the variance/covariance matrix #
repM <- rep("M",10)
repD <- rep("D",10)
row <- c(1,2,3,4,2,3,4,3,4,4)
col <- c(1,1,1,1,2,2,2,3,3,4)
labM <- paste(repM,row,col,sep="")
meanM <- paste("meanM",c(1,2),c("A","A","B","B"),sep="")
meanD <- paste("meanD",c(1,2),c("A","A","B","B"),sep="")
labD <- paste(repD,row,col,sep="")

multiTwinSatModel <- mxModel("multiTwinSat",
    mxModel("MZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=cholSV, lbound=-4, ubound=4, labels=labM, name="CholMZ" ),
        mxAlgebra( expression=CholMZ %*% t(CholMZ), name="expCovMZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanSV, labels=meanM, name="expMeanMZ" ),
        mxData( observed=MZdata, type="raw" ),
	  mxExpectationNormal( covariance="MZ.expCovMZ", means="MZ.expMeanMZ", dimnames=selvars),  
	  mxFitFunctionML(),

    # Algebra's needed for equality constraints    
        mxAlgebra( expression=t(diag2vec(expCovMZ)), name="expVarMZ"),
        mxAlgebra( expression=expVarMZ[1,1:nv], name="expVarMZtA"),
        mxAlgebra( expression=expVarMZ[1,(nv+1):ntv], name="expVarMZtB"),
        mxAlgebra( expression=expMeanMZ[1,1:nv], name="expMeanMZtA"),
        mxAlgebra( expression=expMeanMZ[1,(nv+1):ntv], name="expMeanMZtB")
    ),
    mxModel("DZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=cholSV, lbound=-4, ubound=4, labels=labD, name="CholDZ" ),
        mxAlgebra( expression=CholDZ %*% t(CholDZ), name="expCovDZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanSV, labels=meanD, name="expMeanDZ" ),
        mxData( observed=DZdata, type="raw" ),
	  mxExpectationNormal( covariance="DZ.expCovDZ", means="DZ.expMeanDZ", dimnames=selvars),  
	  mxFitFunctionML(),

    # Algebra's needed for equality constraints    
        mxAlgebra( expression=t(diag2vec(expCovDZ)), name="expVarDZ"),
        mxAlgebra( expression=expVarDZ[1,1:nv], name="expVarDZtA"),
        mxAlgebra( expression=expVarDZ[1,(nv+1):ntv], name="expVarDZtB"),
        mxAlgebra( expression=expMeanDZ[1,1:nv], name="expMeanDZtA"),
        mxAlgebra( expression=expMeanDZ[1,(nv+1):ntv], name="expMeanDZtB")
    ),
    mxAlgebra( MZ.objective + DZ.objective, name="neg2sumll" ),
    mxFitFunctionAlgebra("neg2sumll")	# Calculates the -2LL
)

multiTwinSatFit <- mxRun(multiTwinSatModel)
multiTwinSatSumm <- summary(multiTwinSatFit)
multiTwinSatSumm

# Generate Saturated Output
parameterSpecifications(multiTwinSatFit)
expectedMeansCovariances(multiTwinSatFit)
tableFitStatistics(multiTwinSatFit)

#-------------------------#
# Fit Bivariate ACE Model #
#-------------------------#
meanSVnv <- rep(0,nv)
cholASVnv <- rep(.3,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholCSVnv <- rep(.1,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholESVnv <- rep(.3,nv*(nv+1)/2)  #Establishes starting values based on the data 

rows <- c()
cols <- c()

# Optional code for creating labels #
for (i in 1:nv){
row <- c(i:nv)
col <- rep(i,(nv+1-i))
rows <- c(rows,row)
cols <- c(cols,col)
}

AFac <- paste("A",rows,cols,sep="")
CFac <- paste("C",rows,cols,sep="")
EFac <- paste("E",rows,cols,sep="")

mean <- as.character(c("mean1","mean2")) # Optional: Provides labels for the means

multiCholACEModel <- mxModel("ACE",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, labels=AFac, name="a", lbound=-3, ubound=3 ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, labels=CFac, name="c", lbound=-3, ubound=3 ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, labels=EFac, name="e", lbound=-3, ubound=3 ),

    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=a %*% t(a), name="A" ),
        mxAlgebra( expression=c %*% t(c), name="C" ),
        mxAlgebra( expression=e %*% t(e), name="E" ),

    # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=A+C+E, name="V" ),
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),

    # Confidence intervals portion for covariance matrices
	  mxAlgebra( expression=A/V,name="stndVCA"),
	  mxAlgebra( expression=C/V,name="stndVCC"),
	  mxAlgebra( expression=E/V,name="stndVCE"),
	  mxCI(c("stndVCA","stndVCC","stndVCE")),

    # Confidence intervals for Phenotypic correlation pieces
	  mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[1,2],name="Vcorr"),
	  mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[1,2],name="Acorr"),
	  mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C)))[1,2],name="Ccorr"),
	  mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[1,2],name="Ecorr"),
	  mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),

## Note that the rest of the mxModel statements do not change for bi/multivariate case
    # Matrix & Algebra for expected means vector
        mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=meanSVnv, labels=mean, name="Mean" ),
        mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
    # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
    # Algebra for expected variance/covariance matrix in DZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ),
    mxModel("MZ",
        mxData( observed=MZdata, type="raw" ),
	  mxExpectationNormal( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selvars),  
	  mxFitFunctionML()
	),

    mxModel("DZ", 
        mxData( observed=DZdata, type="raw" ),
	  mxExpectationNormal( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selvars),  
	  mxFitFunctionML()
	),
    mxAlgebra( expression=MZ.objective + DZ.objective, name="neg2sumll" ),
    mxFitFunctionAlgebra("neg2sumll")
)
multiCholACEFit <- mxRun(multiCholACEModel,intervals=FALSE)
multiCholACESumm <- summary(multiCholACEFit)
multiCholACESumm

# Generate Bivariate Cholesky ACE Output
parameterSpecifications(multiCholACEFit)
expectedMeansCovariances(multiCholACEFit)
tableFitStatistics(multiTwinSatFit, multiCholACEFit)

# Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices

Vars <- c("var")

# Cholesky Output (Path Estimates)
ACEpathMatrices <- c("ACE.a","ACE.c","ACE.e","ACE.iSD","ACE.iSD %*% ACE.a","ACE.iSD %*% ACE.c","ACE.iSD %*% ACE.e")
ACEpathLabels <- c("pathEst_a","pathEst_c","pathEst_e","sd","stPathEst_a","stPathEst_c","stPathEst_e")
formatOutputMatrices(multiCholACEFit,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
ACEcovMatrices <- c("ACE.A","ACE.C","ACE.E","ACE.V","ACE.A/ACE.V","ACE.C/ACE.V","ACE.E/ACE.V")
ACEcovLabels <- c("covComp_A","covComp_C","covComp_E","Var","stCovComp_A","stCovComp_C","stCovComp_E")
formatOutputMatrices(multiCholACEFit,ACEcovMatrices,ACEcovLabels,Vars,4)

#-----------------#
## Nested Models ##
#-----------------#

# Set C = 0, Bivariate AE Model
# -----------------------------------------------------------------------

multiCholAEModel <- mxRename(multiCholACEModel, "AE")
multiCholAEModel$AE.c <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=FALSE, values=0, name="c" )
multiCholAEModel$AE.Ccorr <- mxAlgebra(0,name="Ccorr")
multiCholAEFit <- mxRun(multiCholAEModel,intervals=T)
multiCholAESumm <- summary(multiCholAEFit)
multiCholAESumm

tableFitStatistics(multiCholACEFit,multiCholAEFit)

# Cholesky Output (Path Estimates)
AEpathMatrices <- c("AE.a","AE.c","AE.e","AE.iSD","AE.iSD %*% AE.a","AE.iSD %*% AE.c","AE.iSD %*% AE.e")
formatOutputMatrices(multiCholAEFit,AEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
AEcovMatrices <- c("AE.A","AE.C","AE.E","AE.V","AE.A/AE.V","AE.C/AE.V","AE.E/AE.V")
formatOutputMatrices(multiCholAEFit,AEcovMatrices,ACEcovLabels,Vars,4)


# Set r_C = 0 
# -----------------

multiCholACEModel_noCcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noCcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c" )
multiCholACEFit_noCcorr <- mxRun(multiCholACEModel_noCcorr)
multiCholACESumm_noCcorr <- summary(multiCholACEFit_noCcorr)
multiCholACESumm_noCcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noCcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noCcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noCcorr,ACEcovMatrices,ACEcovLabels,Vars,4)


# Set r_A = 0 
# -----------------------------------------------------------------------

multiCholACEModel_noAcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noAcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
multiCholACEFit_noAcorr <- mxRun(multiCholACEModel_noAcorr)
multiCholACESumm_noAcorr <- summary(multiCholACEFit_noAcorr)
multiCholACESumm_noAcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noAcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noAcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noAcorr,ACEcovMatrices,ACEcovLabels,Vars,4)


# Set r_E = 0 
# -----------------------------------------------------------------------

multiCholACEModel_noEcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
multiCholACEFit_noEcorr <- mxRun(multiCholACEModel_noEcorr)
multiCholACESumm_noEcorr <- summary(multiCholACEFit_noEcorr)
multiCholACESumm_noEcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noEcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noEcorr,ACEcovMatrices,ACEcovLabels,Vars,4)

# Set r_A & r_C = 0
# -----------------------------------------------------------------------

multiCholACEModel_noACcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noACcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
multiCholACEModel_noACcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c")
multiCholACEFit_noACcorr <- mxRun(multiCholACEModel_noACcorr)
multiCholACESumm_noACcorr <- summary(multiCholACEFit_noACcorr)
multiCholACESumm_noACcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noACcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noACcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noACcorr,ACEcovMatrices,ACEcovLabels,Vars,4)

# Set r_A & r_E = 0
# -----------------------------------------------------------------------

multiCholACEModel_noAEcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noAEcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
multiCholACEModel_noAEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
multiCholACEFit_noAEcorr <- mxRun(multiCholACEModel_noACcorr)
multiCholACESumm_noAEcorr <- summary(multiCholACEFit_noACcorr)
multiCholACESumm_noAEcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noAEcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noAEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noACcorr,ACEcovMatrices,ACEcovLabels,Vars,4)

# Set r_E & r_C = 0
# -----------------------------------------------------------------------

multiCholACEModel_noCEcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noCEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
multiCholACEModel_noCEcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c")
multiCholACEFit_noCEcorr <- mxRun(multiCholACEModel_noCEcorr)
multiCholACESumm_noCEcorr <- summary(multiCholACEFit_noCEcorr)
multiCholACESumm_noCEcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noCEcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noCEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noCEcorr,ACEcovMatrices,ACEcovLabels,Vars,4)

# Test overall phenotypic correlation = 0 (i.e, r_A = r_C = r_E = 0)
# -----------------------------------------------------------------------

multiCholACEModel_nocorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_nocorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
multiCholACEModel_nocorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c")
multiCholACEModel_nocorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
multiCholACEFit_nocorr <- mxRun(multiCholACEModel_nocorr)
multiCholACESumm_nocorr <- summary(multiCholACEFit_nocorr)
multiCholACESumm_nocorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_nocorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_nocorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_nocorr,ACEcovMatrices,ACEcovLabels,Vars,4)




#-------------------------------------------------------------------#
### OUTPUT SUMMARY ###
#-------------------------------------------------------------------#

sat2ll <- multiTwinSatSumm$Minus2LogLikelihood
satdf <- multiTwinSatFit$degreesOfFreedom

#ACE
IFAIL<- multiCholACEFit@output$status[[1]]
ACEfit<- multiCholACEFit@output$Minus2LogLikelihood
AICace<- multiCholACESumm$AIC
x2acefit<- ACEfit-sat2ll
pacefit<- pchisq(ACEfit-sat2ll,lower.tail=F,17)

# Test of no A corr
ACEfit_noAcorr<- multiCholACEFit_noAcorr@output$Minus2LogLikelihood
x2noAcorr<- ACEfit_noAcorr-ACEfit
pnoAcorr<- pchisq(ACEfit_noAcorr-ACEfit,lower.tail=F,1)
AICace_noAcorr<- multiCholACESumm_noAcorr$AIC

# Test of no C corr
ACEfit_noCcorr<-multiCholACEFit_noCcorr@output$Minus2LogLikelihood
x2noCcorr<- ACEfit_noCcorr-ACEfit
pnoCcorr<- pchisq(ACEfit_noCcorr-ACEfit,lower.tail=F,1)
AICace_noCcorr<- multiCholACESumm_noCcorr$AIC

#AC (test of no E corr)
ACEfit_noEcorr<-multiCholACEFit_noEcorr@output$Minus2LogLikelihood
x2noEcorr<- ACEfit_noEcorr-ACEfit
pnoEcorr<- pchisq(ACEfit_noEcorr-ACEfit,lower.tail=F,1)
AICace_noEcorr<- multiCholACESumm_noEcorr$AIC

#E Only (test of no A or C correlation)
ACEfit_noACcorr<-multiCholACEFit_noACcorr@output$Minus2LogLikelihood
x2noACcorr<- ACEfit_noACcorr-ACEfit
pnoACcorr<- pchisq(ACEfit_noACcorr-ACEfit,lower.tail=F,2)
AICace_noACcorr<- multiCholACESumm_noACcorr$AIC

#C Only (test of no A or E correlation)
ACEfit_noAEcorr<-multiCholACEFit_noAEcorr@output$Minus2LogLikelihood
x2noAEcorr<- ACEfit_noAEcorr-ACEfit
pnoAEcorr<- pchisq(ACEfit_noAEcorr-ACEfit,lower.tail=F,2)
AICace_noAEcorr<- multiCholACESumm_noAEcorr$AIC

#A Only (test of no C or E correlation)
ACEfit_noCEcorr<-multiCholACEFit_noCEcorr@output$Minus2LogLikelihood
x2noCEcorr<- ACEfit_noCEcorr-ACEfit
pnoCEcorr<- pchisq(ACEfit_noCEcorr-ACEfit,lower.tail=F,2)
AICace_noCEcorr<- multiCholACESumm_noCEcorr$AIC

#No Phenotypic Corr
ACEfit_nocorr<-multiCholACEFit_nocorr@output$Minus2LogLikelihood
x2nocorr<- ACEfit_nocorr-ACEfit
pnocorr<- pchisq(ACEfit_nocorr-ACEfit,lower.tail=F,2)
AICace_nocorr<- multiCholACESumm_nocorr$AIC

##Combine Output
bivariate<-c(sat2ll,ACEfit,x2acefit,pacefit,AICace,ACEfit_noAcorr,x2noAcorr,pnoAcorr,AICace_noAcorr,
ACEfit_noCcorr,x2noCcorr,pnoCcorr,AICace_noCcorr,ACEfit_noEcorr,x2noEcorr,pnoEcorr,AICace_noEcorr,
ACEfit_noACcorr,x2noACcorr,pnoACcorr,AICace_noACcorr,
ACEfit_noAEcorr,x2noAEcorr,pnoAEcorr,AICace_noAEcorr,
ACEfit_noCEcorr,x2noCEcorr,pnoCEcorr,AICace_noCEcorr,
ACEfit_nocorr,x2nocorr,pnocorr,AICace_nocorr)

names(bivariate)<-c("satfit","ACEfit","x2acefit","pacefit","AICace",
"ACEfit_noAcorr","x2noAcorr","pnoAcorr","AICace_noAcorr",
"ACEfit_noCcorr","x2noCcorr","pnoCcorr","AICace_noCcorr",
"ACEfit_noEcorr","x2noEcorr","pnoEcorr","AICeace_noEcorr",
"ACEfit_noACcorr","x2noACcorr","pnoACcorr","AICace_noACcorr",
"ACEfit_noAEcorr","x2noAEcorr","pnoAEcorr","AICace_noAEcorr",
"ACEfit_noCEcorr","x2noCEcorr","pnoCEcorr","AICeace_noCEcorr",
"ACEfit_nocorr","x2nocorr","pnocorr","AICace_nocorr")

bivariate # Provides all model fit results for the entire script



