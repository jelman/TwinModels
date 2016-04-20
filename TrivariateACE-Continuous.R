#-----------------------------------------------------------------------#
# Program: Trivariate Saturated and ACE Model				                    #
#  Author: Kelly Spoon                                                  #
#    Date: 07 18 2011                                                   #
#Modified: 10 20 2015 by Jeremy Elman                                   #
#                                                                       #
# Description: Complete trivariate analysis (Saturated, ACE model, 	    #
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
#      Output: All trivariate statistics                                #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user       #
#-----------------------------------------------------------------------#

## Loading Required/Useful Libraries ##

# Load OpenMX
require(OpenMx)   	#Loads OpenMx
require(psych)   		#Loads Psych package
library(dplyr)      #Loads dplyr package
library(tidyr)      #Loads tidyr package
library(stringr)    #Loads stringr package

# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL")  

# Load additional GenEpi tools
source("K:/code/TwinModels/GenEpiHelperFunctions.R")

#*** Set working directory - location of data file and desired output
setwd("K:/Projects/Cortical_DTI/data/Merged/")

#*** Read in data file
	twins <- read.csv("V2_MD_CT_WMMD_GlobalAgeSiteAdj.csv",header=T)

#*** Quick check that data file is properly read in
	dim(twins)[1] 		# should be number of subjects
	dim(twins)[2] 		# should be number of variables in data file
	names(twins)[1:4] 	# should be subid, case, twin, zygosity
	names(twins)[7:8] 	# should be the phenotypes of interest


#*** Defining number of variables
nv <- 3				#Number of phenotypes (NVAR)
ntv <- nv*2

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg14"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
	names(newtwins) 		# Gives the variable names for the new paired dataset

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg14==1,c(5,6,7,10,11,12)))
	DZdata <- as.data.frame(subset(newtwins,zyg14==2,c(5,6,7,10,11,12)))

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
## Fit Trivariate Saturated Model ##
#---------------------------------#

meanSV <- rep(0, 6)  				# Start values for the means
cholSV <- diag(x=.9, 6,6)[lower.tri(diag(6), diag=TRUE)] 	# Start values for the lower var/cov matrix 

# The following lines are optional, but are useful for labeling elements of the variance/covariance matrix #
repM <- rep("M",length(cholSV))
repD <- rep("D",length(cholSV))
rows <- c()
cols <- c()

# Optional code for creating labels #
for (i in 1:ntv){
  row <- c(i:ntv)
  col <- rep(i,(ntv+1-i))
  rows <- c(rows,row)
  cols <- c(cols,col)
}
labM <- paste(repM,rows,cols,sep="")
meanM <- paste("meanM",c(1,2,3),c("A","A","A","B","B","B"),sep="")
meanD <- paste("meanD",c(1,2,3),c("A","A","A","B","B","B"),sep="")
labD <- paste(repD,rows,cols,sep="")

multiTwinSatModel <- mxModel("multiTwinSat",
    mxModel("MZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=cholSV, lbound=-4, ubound=4, 
                  labels=labM, name="CholMZ" ),
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
# Fit Trivariate ACE Model #
#-------------------------#
meanSVnv <- rep(0,nv)
cholASVnv <- diag(.3, nrow=nv, ncol=nv)  #Establishes starting values based on the data 
cholCSVnv <- diag(.1, nrow=nv, ncol=nv)  #Establishes starting values based on the data 
cholESVnv <- diag(.3, nrow=nv, ncol=nv)  #Establishes starting values based on the data 

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

mean <- as.character(paste0("mean",seq(1:nv))) # Optional: Provides labels for the means

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
    mxAlgebra( expression=(A/V)[1,1],name="stndVCA_V1"),
    mxAlgebra( expression=(A/V)[2,2],name="stndVCA_V2"),
    mxAlgebra( expression=(A/V)[3,3],name="stndVCA_V3"),
    mxAlgebra( expression=(C/V)[1,1],name="stndVCC_V1"),
    mxAlgebra( expression=(C/V)[2,2],name="stndVCC_V2"),
    mxAlgebra( expression=(C/V)[3,3],name="stndVCC_V3"),    
    mxAlgebra( expression=(E/V)[1,1],name="stndVCE_V1"),
    mxAlgebra( expression=(E/V)[2,2],name="stndVCE_V2"),
    mxAlgebra( expression=(E/V)[3,3],name="stndVCE_V3"),
    mxCI(c("stndVCA_V1","stndVCA_V2","stndVCA_V3",
           "stndVCC_V1","stndVCC_V2","stndVCC_V3",
           "stndVCE_V1","stndVCE_V2","stndVCE_V3")),
    # Confidence intervals for Phenotypic correlation pieces
    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[2,1],name="Vcorr_V1V2"),
    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[3,1],name="Vcorr_V1V3"),
    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[3,2],name="Vcorr_V2V3"),
    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[2,1],name="Acorr_V1V2"),
    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[3,1],name="Acorr_V1V3"),
    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[3,2],name="Acorr_V2V3"),
    mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C)))[2,1],name="Ccorr_V1V2"),
    mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C)))[3,1],name="Ccorr_V1V3"),
    mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C)))[3,2],name="Ccorr_V2V3"),
    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[2,1],name="Ecorr_V1V2"),
    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[3,1],name="Ecorr_V1V3"),
    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[3,2],name="Ecorr_V2V3"),
    mxCI(c("Vcorr_V1V2","Vcorr_V1V3","Vcorr_V2V3",
           "Acorr_V1V2","Acorr_V1V3","Acorr_V2V3",
           "Ccorr_V1V2","Ccorr_V1V3","Ccorr_V2V3",
           "Ecorr_V1V2","Ecorr_V1V3","Ecorr_V2V3")),
    
    # Confidence intervals for paths
    mxCI(c("a","c","e")),


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
multiCholACEFit <- mxRun(multiCholACEModel,intervals=TRUE)
multiCholACESumm <- summary(multiCholACEFit)
multiCholACESumm

# Generate Trivariate Cholesky ACE Output
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


# Set C = 0, Trivariate AE Model
# -----------------------------------------------------------------------

# multiCholAEModel <- mxRename(multiCholACEModel, "AE")
# multiCholAEModel$AE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, name="c" )
# multiCholAEModel$AE.Ccorr <- mxAlgebra(0,name="Ccorr")
# multiCholAEFit <- mxRun(multiCholAEModel,intervals=TRUE)
# multiCholAESumm <- summary(multiCholAEFit, verbose=T)
# multiCholAESumm

multiCholAEModel <- mxModel("AE",
    # Matrices a, c, and e to store a, c, and e path coefficients
    mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, labels=AFac, name="a", lbound=-3, ubound=3 ),
    mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, labels=CFac, name="c", lbound=-3, ubound=3 ),
    mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, labels=EFac, name="e", lbound=-3, ubound=3 ),
    
    # Matrices A, and E compute variance components
    mxAlgebra( expression=a %*% t(a), name="A" ),
    mxAlgebra( expression=e %*% t(e), name="E" ),
    
    # Algebra to compute total variances and standard deviations (diagonal only)
    mxAlgebra( expression=A+E, name="V" ),
    mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
    mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
    
    # Confidence intervals portion for covariance matrices
    mxAlgebra( expression=(A/V)[1,1],name="stndVCA_V1"),
    mxAlgebra( expression=(A/V)[2,2],name="stndVCA_V2"),
    mxAlgebra( expression=(A/V)[3,3],name="stndVCA_V3"),
    mxAlgebra( expression=(E/V)[1,1],name="stndVCE_V1"),
    mxAlgebra( expression=(E/V)[2,2],name="stndVCE_V2"),
    mxAlgebra( expression=(E/V)[3,3],name="stndVCE_V3"),
    mxCI(c("stndVCA_V1","stndVCA_V2","stndVCA_V3",
           "stndVCE_V1","stndVCE_V2","stndVCE_V3")),
    # Confidence intervals for Phenotypic correlation pieces
    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[2,1],name="Vcorr_V1V2"),
    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[3,1],name="Vcorr_V1V3"),
    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[3,2],name="Vcorr_V2V3"),
    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[2,1],name="Acorr_V1V2"),
    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[3,1],name="Acorr_V1V3"),
    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[3,2],name="Acorr_V2V3"),
    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[2,1],name="Ecorr_V1V2"),
    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[3,1],name="Ecorr_V1V3"),
    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[3,2],name="Ecorr_V2V3"),
    mxCI(c("Vcorr_V1V2","Vcorr_V1V3","Vcorr_V2V3",
           "Acorr_V1V2","Acorr_V1V3","Acorr_V2V3",
           "Ecorr_V1V2","Ecorr_V1V3","Ecorr_V2V3")),
    
    # Confidence intervals for paths
    mxCI(c("a","e")),
    
    
    ## Note that the rest of the mxModel statements do not change for bi/multivariate case
    # Matrix & Algebra for expected means vector
    mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=meanSVnv, labels=mean, name="Mean" ),
    mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
    # Algebra for expected variance/covariance matrix in MZ
    mxAlgebra( expression= rbind  ( cbind(A+E , A),
                                    cbind(A   , A+E)), name="expCovMZ" ),
    # Algebra for expected variance/covariance matrix in DZ
    mxAlgebra( expression= rbind  ( cbind(A+E     , 0.5%x%A),
                                    cbind(0.5%x%A , A+E)),  name="expCovDZ" ),
    mxModel("MZ",
            mxData( observed=MZdata, type="raw" ),
            mxExpectationNormal( covariance="AE.expCovMZ", means="AE.expMean", dimnames=selvars),  
            mxFitFunctionML()
    ),
    
    mxModel("DZ", 
            mxData( observed=DZdata, type="raw" ),
            mxExpectationNormal( covariance="AE.expCovDZ", means="AE.expMean", dimnames=selvars),  
            mxFitFunctionML()
    ),
    mxAlgebra( expression=MZ.objective + DZ.objective, name="neg2sumll" ),
    mxFitFunctionAlgebra("neg2sumll")
)
multiCholAEFit <- mxRun(multiCholAEModel,intervals=TRUE)
multiCholAESumm <- summary(multiCholAEFit)
multiCholAESumm

tableFitStatistics(multiCholACEFit,multiCholAEFit)

# Cholesky Output (Path Estimates)
AEpathMatrices <- c("AE.a","AE.e","AE.iSD","AE.iSD %*% AE.a","AE.iSD %*% AE.e")
formatOutputMatrices(multiCholAEFit,AEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
AEcovMatrices <- c("AE.A","AE.E","AE.V","AE.A/AE.V","AE.E/AE.V")
formatOutputMatrices(multiCholAEFit,AEcovMatrices,ACEcovLabels,Vars,4)


# Set r_C = 0 
# -----------------

multiCholACEModel_noCcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noCcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                             free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                             values=cholCSVnv, labels=CFac, name="c" )
multiCholACEFit_noCcorr <- mxRun(multiCholACEModel_noCcorr)
multiCholACESumm_noCcorr <- summary(multiCholACEFit_noCcorr, intervals=FALSE)
multiCholACESumm_noCcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noCcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noCcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noCcorr,ACEcovMatrices,ACEcovLabels,Vars,4)


# Set r_A = 0 
# -----------------------------------------------------------------------

multiCholACEModel_noAcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noAcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                             free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                             values=cholASVnv, labels=AFac, name="a")
multiCholACEFit_noAcorr <- mxRun(multiCholACEModel_noAcorr, intervals=FALSE)
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
multiCholACEModel_noEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                             free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                             values=cholESVnv, labels=EFac, name="e")
multiCholACEFit_noEcorr <- mxRun(multiCholACEModel_noEcorr)
multiCholACESumm_noEcorr <- summary(multiCholACEFit_noEcorr, intervals=FALSE)
multiCholACESumm_noEcorr
tableFitStatistics(multiCholACEFit,multiCholACEFit_noEcorr)

# Cholesky Output (Path Estimates)
formatOutputMatrices(multiCholACEFit_noEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)

# Correlated Factors Output (Cov and Corr Matrices)
formatOutputMatrices(multiCholACEFit_noEcorr,ACEcovMatrices,ACEcovLabels,Vars,4)

# Set r_A & r_C = 0
# -----------------------------------------------------------------------

multiCholACEModel_noACcorr <- mxRename(multiCholACEModel, "ACE")
multiCholACEModel_noACcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                              free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                              values=cholASVnv, labels=AFac, name="a")
multiCholACEModel_noACcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                              free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                              values=cholCSVnv, labels=CFac, name="c")
multiCholACEFit_noACcorr <- mxRun(multiCholACEModel_noACcorr, intervals=FALSE)
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
multiCholACEModel_noAEcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                              free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                              values=cholASVnv, labels=AFac, name="a")
multiCholACEModel_noAEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                              free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                              values=cholESVnv, labels=EFac, name="e")
multiCholACEFit_noAEcorr <- mxRun(multiCholACEModel_noACcorr, intervals=FALSE)
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
multiCholACEModel_noCEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                              free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                              values=cholESVnv, labels=EFac, name="e")
multiCholACEModel_noCEcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                              free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                              values=cholCSVnv, labels=CFac, name="c")
multiCholACEFit_noCEcorr <- mxRun(multiCholACEModel_noCEcorr, intervals=FALSE)
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
multiCholACEModel_nocorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                            free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                            values=cholASVnv, labels=AFac, name="a")
multiCholACEModel_nocorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                            free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                            values=cholCSVnv, labels=CFac, name="c")
multiCholACEModel_nocorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, 
                                            free=c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE), 
                                            values=cholESVnv, labels=EFac, name="e")
multiCholACEFit_nocorr <- mxRun(multiCholACEModel_nocorr, intervals=FALSE)
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
satdf <- multiTwinSatSumm$degreesOfFreedom

#ACE
IFAIL<- multiCholACEFit@output$status[[1]]
ACEfit<- multiCholACEFit@output$Minus2LogLikelihood
ACEdf<- multiCholACESumm$degreesOfFreedom
AICace<- multiCholACESumm$AIC
x2acefit<- ACEfit-sat2ll
x2df_ACE<- ACEdf - satdf
pacefit<- pchisq(ACEfit-sat2ll,lower.tail=F,x2df_ACE)

# Test of no A corr
ACEfit_noAcorr<- multiCholACEFit_noAcorr@output$Minus2LogLikelihood
x2noAcorr<- ACEfit_noAcorr-ACEfit
ACEdf_noAcorr<- multiCholACESumm_noAcorr$degreesOfFreedom
x2df_noAcorr<- ACEdf_noAcorr - ACEdf
pnoAcorr<- pchisq(ACEfit_noAcorr-ACEfit,lower.tail=F,x2df_noAcorr)
AICace_noAcorr<- multiCholACESumm_noAcorr$AIC

# Test of no C corr
ACEfit_noCcorr<-multiCholACEFit_noCcorr@output$Minus2LogLikelihood
x2noCcorr<- ACEfit_noCcorr-ACEfit
ACEdf_noCcorr<- multiCholACESumm_noCcorr$degreesOfFreedom
x2df_noCcorr<- ACEdf_noCcorr - ACEdf
pnoCcorr<- pchisq(ACEfit_noCcorr-ACEfit,lower.tail=F,x2df_noCcorr)
AICace_noCcorr<- multiCholACESumm_noCcorr$AIC

#AC (test of no E corr)
ACEfit_noEcorr<-multiCholACEFit_noEcorr@output$Minus2LogLikelihood
x2noEcorr<- ACEfit_noEcorr-ACEfit
ACEdf_noEcorr<- multiCholACESumm_noEcorr$degreesOfFreedom
x2df_noEcorr<- ACEdf_noEcorr - ACEdf
pnoEcorr<- pchisq(ACEfit_noEcorr-ACEfit,lower.tail=F,x2df_noEcorr)
AICace_noEcorr<- multiCholACESumm_noEcorr$AIC

#E Only (test of no A or C correlation)
ACEfit_noACcorr<-multiCholACEFit_noACcorr@output$Minus2LogLikelihood
x2noACcorr<- ACEfit_noACcorr-ACEfit
ACEdf_noACcorr<- multiCholACESumm_noACcorr$degreesOfFreedom
x2df_noACcorr<- ACEdf_noACcorr - ACEdf
pnoACcorr<- pchisq(ACEfit_noACcorr-ACEfit,lower.tail=F,x2df_noACcorr)
AICace_noACcorr<- multiCholACESumm_noACcorr$AIC

#C Only (test of no A or E correlation)
ACEfit_noAEcorr<-multiCholACEFit_noAEcorr@output$Minus2LogLikelihood
x2noAEcorr<- ACEfit_noAEcorr-ACEfit
ACEdf_noAEcorr<- multiCholACESumm_noAEcorr$degreesOfFreedom
x2df_noAEcorr<- ACEdf_noAEcorr - ACEdf
pnoAEcorr<- pchisq(ACEfit_noAEcorr-ACEfit,lower.tail=F,x2df_noAEcorr)
AICace_noAEcorr<- multiCholACESumm_noAEcorr$AIC

#A Only (test of no C or E correlation)
ACEfit_noCEcorr<-multiCholACEFit_noCEcorr@output$Minus2LogLikelihood
x2noCEcorr<- ACEfit_noCEcorr-ACEfit
ACEdf_noCEcorr<- multiCholACESumm_noCEcorr$degreesOfFreedom
x2df_noCEcorr<- ACEdf_noCEcorr - ACEdf
pnoCEcorr<- pchisq(ACEfit_noCEcorr-ACEfit,lower.tail=F,x2df_noCEcorr)
AICace_noCEcorr<- multiCholACESumm_noCEcorr$AIC

#No Phenotypic Corr
ACEfit_nocorr<-multiCholACEFit_nocorr@output$Minus2LogLikelihood
x2nocorr<- ACEfit_nocorr-ACEfit
ACEdf_nocorr<- multiCholACESumm_nocorr$degreesOfFreedom
x2df_nocorr<- ACEdf_nocorr - ACEdf
pnocorr<- pchisq(ACEfit_nocorr-ACEfit,lower.tail=F,x2df_nocorr)
AICace_nocorr<- multiCholACESumm_nocorr$AIC

##Combine Output
Trivariate<-c(sat2ll,ACEfit,x2acefit,pacefit,AICace,ACEfit_noAcorr,x2noAcorr,pnoAcorr,AICace_noAcorr,
ACEfit_noCcorr,x2noCcorr,pnoCcorr,AICace_noCcorr,ACEfit_noEcorr,x2noEcorr,pnoEcorr,AICace_noEcorr,
ACEfit_noACcorr,x2noACcorr,pnoACcorr,AICace_noACcorr,
ACEfit_noAEcorr,x2noAEcorr,pnoAEcorr,AICace_noAEcorr,
ACEfit_noCEcorr,x2noCEcorr,pnoCEcorr,AICace_noCEcorr,
ACEfit_nocorr,x2nocorr,pnocorr,AICace_nocorr)

names(Trivariate)<-c("satfit","ACEfit","x2acefit","pacefit","AICace",
"ACEfit_noAcorr","x2noAcorr","pnoAcorr","AICace_noAcorr",
"ACEfit_noCcorr","x2noCcorr","pnoCcorr","AICace_noCcorr",
"ACEfit_noEcorr","x2noEcorr","pnoEcorr","AICeace_noEcorr",
"ACEfit_noACcorr","x2noACcorr","pnoACcorr","AICace_noACcorr",
"ACEfit_noAEcorr","x2noAEcorr","pnoAEcorr","AICace_noAEcorr",
"ACEfit_noCEcorr","x2noCEcorr","pnoCEcorr","AICeace_noCEcorr",
"ACEfit_nocorr","x2nocorr","pnocorr","AICace_nocorr")

Trivariate # Provides all model fit results for the entire script



