#-----------------------------------------------------------------------#
# Program: Multivariate ACE Model                           		        #
# Contents: Saturated Model, Cholesky Model, and Common Pathways       	#
# Created: December 5, 2016			  		                                  #
# Edited: July 25, 2017 by Jeremy Elman		                              #
# Data Type: Continuous						                                      #
#											                                                  #
#-----------------------------------------------------------------------#

### Reading in Original Data ###

#*** Set working directory - location of data file and desired output
setwd("~/netshare/K/Projects/HippoSubfields/")

## Loading Required/Useful Libraries ##
require(OpenMx)		# Loads OpenMX
require(psych)		# Loads Psych package
require(dplyr)
source("http://www.vipbg.vcu.edu/~vipbg/Tc24/GenEpiHelperFunctions.R") #Loads GenEpi package from VCU

# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL") 

### Reading in Original Data ###
#*** Read in data file
twins <- read.csv("./data/HippoSubfieldsAdjusted_07252017.txt",header=T)



#*** Recode twin from 1/2 to A/B
twins$twin = gsub("1","A",twins$twin)
twins$twin = gsub("2","B",twins$twin)

# Quick check that data file is properly read in
	names(twins) 	  # list variables in data file
	dim(twins)[1] 	# number of subjects
	dim(twins)[2] 	# number of variables in the file


## Preparing the data for OpenMx ##

# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zygos"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
	names(newtwins)   # this is what your paired dataset now looks like
	dim(newtwins)[1] 	# number of pairs (complete and incomplete)
	dim(newtwins)[2] 	# number of variables in the file

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zygos==1))
	DZdata <- as.data.frame(subset(newtwins,zygos==2))
	dim(MZdata)[1]
	dim(DZdata)[1]


# Defining number of phenotype for the model
nv <- 12		# Number of phenotypes 
ntv <- nv*2		# Number of variables (phenotypes x 2)

names(newtwins)

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zygos==1,c(5:16,19:30)))		
	DZdata <- as.data.frame(subset(newtwins,zygos==2,c(5:16,19:30)))		

selvars <- names(MZdata)
selvars

sfVars = grep("_A", selvars, value = T)
sfVars = gsub("_A","",sfVars)
sfVars = gsub("bilat_","",sfVars)

# Print Descriptive Statistics
describe(MZdata)		 
describe(DZdata)

cor(MZdata,use="complete")
cor(DZdata,use="complete")


# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL") 

#-----------------------------------------------------------------------#

## Fit Multivariate Saturated Model ##

#-----------------------------------------------------------------------#


multiSatModel <- mxModel("multiSat",
    mxModel("MZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=.4, lbound=-1, ubound=1, name="CholMZ" ),
        mxAlgebra( expression=CholMZ %*% t(CholMZ), name="expCovMZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=0, name="expMeanMZ" ),
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
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=.4, lbound=-1, ubound=1, name="CholDZ" ),
        mxAlgebra( expression=CholDZ %*% t(CholDZ), name="expCovDZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=0, name="expMeanDZ" ),
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

multiSatFit <- mxRun(multiSatModel)
multiSatSumm <- summary(multiSatFit)
multiSatSumm

tableFitStatistics(multiSatFit)

#-----------------------------------------------------------------------#

## Fit Multivariate ACE Model ##

#-----------------------------------------------------------------------#


meanSVnv <- rep(0.0,nv)
cholASVnv <- rep(.4,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholCSVnv <- rep(.1,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholESVnv <- rep(.4,nv*(nv+1)/2)  #Establishes starting values based on the data 

multiCholACEModel <- mxModel("ACE",
    # Lower matrices for path coefficients / parameter estimates
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, name="a", lbound=-1, ubound=1 ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, name="c", lbound=-1, ubound=1 ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, name="e", lbound=-1, ubound=1 ),

    # Full matrices for genetic and environmental variance
        mxAlgebra( expression=a %*% t(a), name="A" ),
        mxAlgebra( expression=c %*% t(c), name="C" ),
        mxAlgebra( expression=e %*% t(e), name="E" ),

    # Algebra to compute total variances and standard deviations 
        mxAlgebra( expression=A+C+E, name="V" ),
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),

    # Confidence intervals portion for covariance matrices
	  mxAlgebra( expression=A/V,name="stndVCA"),
	  mxAlgebra( expression=C/V,name="stndVCC"),
	  mxAlgebra( expression=E/V,name="stndVCE"),
	  #mxCI(c("stndVCA","stndVCC","stndVCE")),
    
    # Algebra to take diagonal of stndVCA, i.e., heritability
    mxAlgebra( expression= diag2vec(stndVCA), name="diagVCA" ),
    mxCI(c("diagVCA")),

    # Confidence intervals for Phenotypic correlation pieces
	  mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V))),name="Vcorr"),
	  mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A))),name="Acorr"),
	  mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C))),name="Ccorr"),
	  mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E))),name="Ecorr"),
	  #mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),

## Note that the rest of the mxModel statements do not change for bi/multivariate case
    # Matrix & Algebra for expected means vector
        mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=meanSVnv, name="Mean" ),
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

tableFitStatistics(multiCholACEFit)

round(multiCholACEFit$ACE.A@result,4)
round(multiCholACEFit$ACE.C@result,4)
round(multiCholACEFit$ACE.E@result,4)

round(multiCholACEFit$ACE.stndVCA@result,4)
round(multiCholACEFit$ACE.stndVCC@result,4)
round(multiCholACEFit$ACE.stndVCE@result,4)

Amat = multiCholACEFit$A$result
rownames(Amat) = sfVars
colnames(Amat) = sfVars
fa.parallel(Amat, fm = 'ml', fa = 'fa', n.obs=130)
fa.diagram(fa(Amat,nfactors = 3,rotate = "varimax",fm="ml", covar = TRUE), cut=.1)


#-----------------------------------------------------------------------#

## Fit Multivariate ACE Model, no A covariance ##

#-----------------------------------------------------------------------#


meanSVnv <- rep(0.0,nv)
cholASVnv <- .4  #Establishes starting values based on the data 
cholCSVnv <- rep(.1,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholESVnv <- rep(.4,nv*(nv+1)/2)  #Establishes starting values based on the data 

multiCholACEModel_NoA <- mxModel("ACE",
                             # Lower matrices for path coefficients / parameter estimates
                             mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, name="a", lbound=-1, ubound=1 ),
                             mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, name="c", lbound=-1, ubound=1 ),
                             mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, name="e", lbound=-1, ubound=1 ),
                             
                             # Full matrices for genetic and environmental variance
                             mxAlgebra( expression=a %*% t(a), name="A" ),
                             mxAlgebra( expression=c %*% t(c), name="C" ),
                             mxAlgebra( expression=e %*% t(e), name="E" ),
                             
                             # Algebra to compute total variances and standard deviations 
                             mxAlgebra( expression=A+C+E, name="V" ),
                             mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
                             mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
                             
                             # Confidence intervals portion for covariance matrices
                             mxAlgebra( expression=A/V,name="stndVCA"),
                             mxAlgebra( expression=C/V,name="stndVCC"),
                             mxAlgebra( expression=E/V,name="stndVCE"),
                             #mxCI(c("stndVCA","stndVCC","stndVCE")),
                             
                             # Algebra to take diagonal of stndVCA, i.e., heritability
                             mxAlgebra( expression= diag2vec(stndVCA), name="diagVCA" ),
                             mxCI(c("diagVCA")),
                             
                             # Confidence intervals for Phenotypic correlation pieces
                             mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V))),name="Vcorr"),
                             mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A))),name="Acorr"),
                             mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C))),name="Ccorr"),
                             mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E))),name="Ecorr"),
                             #mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),
                             
                             ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                             # Matrix & Algebra for expected means vector
                             mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=meanSVnv, name="Mean" ),
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

multiCholACEFit_NoA <- mxRun(multiCholACEModel_NoA,intervals=FALSE)
tableFitStatistics(multiCholACEFit, multiCholACEFit_NoA)
multiCholACESumm_NoA = summary(multiCholACEFit_NoA)


#-----------------------------------------------------------------------#

## Fit Multivariate ACE Model, no C covariance ##

#-----------------------------------------------------------------------#


meanSVnv <- rep(0.0,nv)
cholASVnv <- rep(.4,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholCSVnv <- .1  #Establishes starting values based on the data 
cholESVnv <- rep(.4,nv*(nv+1)/2)  #Establishes starting values based on the data 

multiCholACEModel_NoC <- mxModel("ACE",
                                 # Lower matrices for path coefficients / parameter estimates
                                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, name="a", lbound=-1, ubound=1 ),
                                 mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, name="c", lbound=-1, ubound=1 ),
                                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, name="e", lbound=-1, ubound=1 ),
                                 
                                 # Full matrices for genetic and environmental variance
                                 mxAlgebra( expression=a %*% t(a), name="A" ),
                                 mxAlgebra( expression=c %*% t(c), name="C" ),
                                 mxAlgebra( expression=e %*% t(e), name="E" ),
                                 
                                 # Algebra to compute total variances and standard deviations 
                                 mxAlgebra( expression=A+C+E, name="V" ),
                                 mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
                                 mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
                                 
                                 # Confidence intervals portion for covariance matrices
                                 mxAlgebra( expression=A/V,name="stndVCA"),
                                 mxAlgebra( expression=C/V,name="stndVCC"),
                                 mxAlgebra( expression=E/V,name="stndVCE"),
                                 #mxCI(c("stndVCA","stndVCC","stndVCE")),
                                 
                                 # Algebra to take diagonal of stndVCA, i.e., heritability
                                 mxAlgebra( expression= diag2vec(stndVCA), name="diagVCA" ),
                                 mxCI(c("diagVCA")),
                                 
                                 # Confidence intervals for Phenotypic correlation pieces
                                 mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V))),name="Vcorr"),
                                 mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A))),name="Acorr"),
                                 mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C))),name="Ccorr"),
                                 mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E))),name="Ecorr"),
                                 #mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),
                                 
                                 ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                                 # Matrix & Algebra for expected means vector
                                 mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=meanSVnv, name="Mean" ),
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

multiCholACEFit_NoC <- mxRun(multiCholACEModel_NoC,intervals=FALSE)
tableFitStatistics(multiCholACEFit, multiCholACEFit_NoC)
multiCholACESumm_NoC = summary(multiCholACEFit_NoC)


#-----------------------------------------------------------------------#

## Fit Multivariate ACE Model, no E covariance ##

#-----------------------------------------------------------------------#


meanSVnv <- rep(0.0,nv)
cholASVnv <- rep(.4,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholCSVnv <- rep(.1,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholESVnv <- .4  #Establishes starting values based on the data 

multiCholACEModel_NoE <- mxModel("ACE",
                                 # Lower matrices for path coefficients / parameter estimates
                                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, name="a", lbound=-1, ubound=1 ),
                                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, name="c", lbound=-1, ubound=1 ),
                                 mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, name="e", lbound=-1, ubound=1 ),
                                 
                                 # Full matrices for genetic and environmental variance
                                 mxAlgebra( expression=a %*% t(a), name="A" ),
                                 mxAlgebra( expression=c %*% t(c), name="C" ),
                                 mxAlgebra( expression=e %*% t(e), name="E" ),
                                 
                                 # Algebra to compute total variances and standard deviations 
                                 mxAlgebra( expression=A+C+E, name="V" ),
                                 mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
                                 mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
                                 
                                 # Confidence intervals portion for covariance matrices
                                 mxAlgebra( expression=A/V,name="stndVCA"),
                                 mxAlgebra( expression=C/V,name="stndVCC"),
                                 mxAlgebra( expression=E/V,name="stndVCE"),
                                 #mxCI(c("stndVCA","stndVCC","stndVCE")),
                                 
                                 # Algebra to take diagonal of stndVCA, i.e., heritability
                                 mxAlgebra( expression= diag2vec(stndVCA), name="diagVCA" ),
                                 mxCI(c("diagVCA")),
                                 
                                 # Confidence intervals for Phenotypic correlation pieces
                                 mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V))),name="Vcorr"),
                                 mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A))),name="Acorr"),
                                 mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C))),name="Ccorr"),
                                 mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E))),name="Ecorr"),
                                 #mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),
                                 
                                 ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                                 # Matrix & Algebra for expected means vector
                                 mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=meanSVnv, name="Mean" ),
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

multiCholACEFit_NoE <- mxRun(multiCholACEModel_NoE,intervals=FALSE)
tableFitStatistics(multiCholACEFit, multiCholACEFit_NoE)
multiCholACESumm_NoE = summary(multiCholACEFit_NoE)


# -----------------------------------------------------------------------
# Single Factor Independent Pathway Model - A only
# -----------------------------------------------------------------------

nf <- 1

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T)
fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP1FAmodel <- mxModel("IndPath1FA",
  mxModel("ACE",
          # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
          mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.3, name="al" ),
          # Algebra for variances components of latent phenotype
          mxAlgebra( expression= al %*% t(al), name="Al" ),
          # Matrix and Algebra for constraint on variance of latent phenotype
          mxAlgebra( expression= diag2vec(Al), name="VarLP" ),
          mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
          mxConstraint(VarLP == Unit),
          # Matrix f for factor loadings on latent phenotype
          mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
          mxCI(c("fl")),
          
          # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
          mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="as" ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="es" ),
          
          # Algebra for residual variance components, with CIs
          mxAlgebra( expression= as %*% t(as), name="As" ),
          mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
          mxAlgebra( expression= es %*% t(es), name="Es" ),
          mxAlgebra( expression= As + Cs + Es, name="Vs"), 
          mxCI(c("Vs","As","Cs","Es")),
          
          # Matrices A, C, and E compute overall variance components
          mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
          mxAlgebra( expression=cs %*% t(cs), name="C" ),
          mxAlgebra( expression=es %*% t(es), name="E" ),
          # Algebra to compute total variances and standard deviations (diagonal only)
          mxAlgebra( expression=A+C+E, name="V" ),
          ## Note that the rest of the mxModel statements do not change for bi/multivariate case
          # Matrix & Algebra for expected means vector
          mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
          mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
          # Algebra for expected variance/covariance matrix in MZ
          mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                          cbind(A+C   , A+C+E)), name="expCovMZ" ),
          # Algebra for expected variance/covariance matrix in DZ
          mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                          cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
  ),
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

IP1FAFit <- mxRun(IP1FAmodel,intervals=FALSE)
IP1FASumm <- summary(IP1FAFit)
IP1FASumm

tableFitStatistics(multiCholACEFit,IP1FAFit)


# -----------------------------------------------------------------------
# Two Factor Independent Pathway Model - A only
# -----------------------------------------------------------------------

nf <- 2

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T,
           F,T,T,T,T,T,T,T,T,T,T,T)

fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP2FAmodel <- mxModel("IndPath2FA",
  mxModel("ACE",
          # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
          mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.3, name="al" ),
          # Algebra for variances components of latent phenotype
          mxAlgebra( expression= al %*% t(al), name="Al" ),
          # Matrix and Algebra for constraint on variance of latent phenotype
          mxAlgebra( expression= diag2vec(Al), name="VarLP" ),
          mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
          mxConstraint(VarLP == Unit),
          # Matrix f for factor loadings on latent phenotype
          mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
          mxCI(c("fl")),
          
          # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
          mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="as" ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="es" ),
          
          # Algebra for residual variance components, with CIs
          mxAlgebra( expression= as %*% t(as), name="As" ),
          mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
          mxAlgebra( expression= es %*% t(es), name="Es" ),
          mxAlgebra( expression= As + Cs + Es, name="Vs"), 
          mxCI(c("Vs","As","Cs","Es")),
          
          # Matrices A, C, and E compute overall variance components
          mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
          mxAlgebra( expression=cs %*% t(cs), name="C" ),
          mxAlgebra( expression=es %*% t(es), name="E" ),
          # Algebra to compute total variances and standard deviations (diagonal only)
          mxAlgebra( expression=A+C+E, name="V" ),
          ## Note that the rest of the mxModel statements do not change for bi/multivariate case
          # Matrix & Algebra for expected means vector
          mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
          mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
          # Algebra for expected variance/covariance matrix in MZ
          mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                          cbind(A+C   , A+C+E)), name="expCovMZ" ),
          # Algebra for expected variance/covariance matrix in DZ
          mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                          cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
  ),
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

IP2FAFit <- mxRun(IP2FAmodel,intervals=FALSE)
IP2FASumm <- summary(IP2FAFit)
IP2FASumm

tableFitStatistics(multiCholACEFit,c(IP1FAFit,IP2FAFit))
tableFitStatistics(IP2FAFit,IP1FAFit)

twoFAfl = IP2FAFit$ACE$fl$values
rownames(twoFAfl) = sfVars
twoFAfl
varimax(twoFAfl)
fa.diagram(varimax(twoFAfl)$loadings, cut=.1)



# -----------------------------------------------------------------------
# Three Factor Independent Pathway Model - A only
# -----------------------------------------------------------------------

nf <- 3

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T,
           F,T,T,T,T,T,T,T,T,T,T,T,
           F,F,T,T,T,T,T,T,T,T,T,T)

fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP3FAmodel <- mxModel("IndPath3FA",
                      mxModel("ACE",
                              # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
                              mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.3, name="al" ),
                              # Algebra for variances components of latent phenotype
                              mxAlgebra( expression= al %*% t(al), name="Al" ),
                              # Matrix and Algebra for constraint on variance of latent phenotype
                              mxAlgebra( expression= diag2vec(Al), name="VarLP" ),
                              mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                              mxConstraint(VarLP == Unit),
                              # Matrix f for factor loadings on latent phenotype
                              mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
                              mxCI(c("fl")),
                              
                              # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                              mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="as" ),
                              mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
                              mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="es" ),
                              
                              # Algebra for residual variance components, with CIs
                              mxAlgebra( expression= as %*% t(as), name="As" ),
                              mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
                              mxAlgebra( expression= es %*% t(es), name="Es" ),
                              mxAlgebra( expression= As + Cs + Es, name="Vs"), 
                              #mxCI(c("Vs","As","Cs","Es")),
                              
                              # Algebra to get diagonal of residual A matrix
                              mxAlgebra( expression= diag2vec(As), name="ResidA" ),
                              mxCI(c("ResidA")),
                              
                              # Matrices A, C, and E compute overall variance components
                              mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
                              mxAlgebra( expression=cs %*% t(cs), name="C" ),
                              mxAlgebra( expression=es %*% t(es), name="E" ),
                              # Algebra to compute total variances and standard deviations (diagonal only)
                              mxAlgebra( expression=A+C+E, name="V" ),
                              ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                              # Matrix & Algebra for expected means vector
                              mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
                              mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
                              # Algebra for expected variance/covariance matrix in MZ
                              mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                                              cbind(A+C   , A+C+E)), name="expCovMZ" ),
                              # Algebra for expected variance/covariance matrix in DZ
                              mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                              cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
                      ),
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

IP3FAFit <- mxRun(IP3FAmodel,intervals=FALSE)
IP3FASumm <- summary(IP3FAFit)
IP3FASumm

tableFitStatistics(multiCholACEFit,c(IP1FAFit,IP2FAFit,IP3FAFit))
tableFitStatistics(IP3FAFit,IP2FAFit)

threeFAfl = IP3FAFit$ACE$fl$values
rownames(threeFAfl) = sfVars
threeFAfl
varimax(threeFAfl)
fa.diagram(varimax(threeFAfl)$loadings, cut=.1)


# -----------------------------------------------------------------------
# Three Factor Independent Pathway Model - A only, no residual A
# -----------------------------------------------------------------------

nf <- 3

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T,
           F,T,T,T,T,T,T,T,T,T,T,T,
           F,F,T,T,T,T,T,T,T,T,T,T)

fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP3FAmodel_noResid <- mxModel("IndPath3FA_noResid",
                              mxModel("ACE",
                                      # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
                                      mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.3, name="al" ),
                                      # Algebra for variances components of latent phenotype
                                      mxAlgebra( expression= al %*% t(al), name="Al" ),
                                      # Matrix and Algebra for constraint on variance of latent phenotype
                                      mxAlgebra( expression= diag2vec(Al), name="VarLP" ),
                                      mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                                      mxConstraint(VarLP == Unit),
                                      # Matrix f for factor loadings on latent phenotype
                                      mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
                                      mxCI(c("fl")),
                                      
                                      # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                                      mxMatrix( type="Diag", nrow=nv, ncol=nv, free=FALSE, values=0, name="as" ),
                                      mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
                                      mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="es" ),
                                      
                                      # Algebra for residual variance components, with CIs
                                      mxAlgebra( expression= as %*% t(as), name="As" ),
                                      mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
                                      mxAlgebra( expression= es %*% t(es), name="Es" ),
                                      mxAlgebra( expression= As + Cs + Es, name="Vs"), 
                                      mxCI(c("Vs","As","Cs","Es")),
                                      
                                      # Matrices A, C, and E compute overall variance components
                                      mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
                                      mxAlgebra( expression=cs %*% t(cs), name="C" ),
                                      mxAlgebra( expression=es %*% t(es), name="E" ),
                                      # Algebra to compute total variances and standard deviations (diagonal only)
                                      mxAlgebra( expression=A+C+E, name="V" ),
                                      ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                                      # Matrix & Algebra for expected means vector
                                      mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
                                      mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
                                      # Algebra for expected variance/covariance matrix in MZ
                                      mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                                                      cbind(A+C   , A+C+E)), name="expCovMZ" ),
                                      # Algebra for expected variance/covariance matrix in DZ
                                      mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                                      cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
                              ),
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

IP3FAFit_noResid <- mxRun(IP3FAmodel_noResid,intervals=FALSE)
IP3FASumm_noResid <- summary(IP3FAFit_noResid)
IP3FASumm_noResid

tableFitStatistics(multiCholACEFit,c(IP1FAFit,IP2FAFit,IP3FAFit, IP3FAFit_noResid))
tableFitStatistics(IP3FAFit,IP3FAFit_noResid)

threeFAfl_noResid = IP3FAFit_noResid$ACE$fl$values
rownames(threeFAfl_noResid) = sfVars
threeFAfl_noResid
varimax(threeFAfl_noResid)
fa.diagram(varimax(threeFAfl_noResid)$loadings, cut=.1)

# -----------------------------------------------------------------------
# Three Factor Independent Pathway Model - A only, partial residual A
# -----------------------------------------------------------------------

nf <- 3

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T,
           F,T,T,T,T,T,T,T,T,T,T,T,
           F,F,T,T,T,T,T,T,T,T,T,T)

fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)


AsTrue = c(F,F,F,F,F,T,F,F,F,F,F,F)
AsLoad = sapply(AsTrue, function(x) ifelse(x==TRUE, .2, 0))


IP3FAmodel_partResid <- mxModel("IndPath3FA_partResid",
                                mxModel("ACE",
                                        # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
                                        mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.3, name="al" ),
                                        # Algebra for variances components of latent phenotype
                                        mxAlgebra( expression= al %*% t(al), name="Al" ),
                                        # Matrix and Algebra for constraint on variance of latent phenotype
                                        mxAlgebra( expression= diag2vec(Al), name="VarLP" ),
                                        mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                                        mxConstraint(VarLP == Unit),
                                        # Matrix f for factor loadings on latent phenotype
                                        mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
                                        mxCI(c("fl")),
                                        
                                        # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                                        mxMatrix( type="Diag", nrow=nv, ncol=nv, free=AsTrue, values=AsLoad, name="as" ),
                                        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
                                        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="es" ),
                                        
                                        # Algebra for residual variance components, with CIs
                                        mxAlgebra( expression= as %*% t(as), name="As" ),
                                        mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
                                        mxAlgebra( expression= es %*% t(es), name="Es" ),
                                        mxAlgebra( expression= As + Cs + Es, name="Vs"), 
                                        #mxCI(c("Vs","As","Cs","Es")),
                                        
                                        # Matrices A, C, and E compute overall variance components
                                        mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
                                        mxAlgebra( expression=cs %*% t(cs), name="C" ),
                                        mxAlgebra( expression=es %*% t(es), name="E" ),
                                        # Algebra to compute total variances and standard deviations (diagonal only)
                                        mxAlgebra( expression=A+C+E, name="V" ),
                                        ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                                        # Matrix & Algebra for expected means vector
                                        mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
                                        mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
                                        # Algebra for expected variance/covariance matrix in MZ
                                        mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
                                        # Algebra for expected variance/covariance matrix in DZ
                                        mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
                                ),
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

IP3FAFit_partResid <- mxRun(IP3FAmodel_partResid,intervals=FALSE)
IP3FASumm_partResid <- summary(IP3FAFit_partResid)
IP3FASumm_partResid

tableFitStatistics(multiCholACEFit,c(IP1FAFit,IP2FAFit,IP3FAFit,IP3FAFit_noResid,IP3FAFit_partResid))
tableFitStatistics(IP3FAFit,IP3FAFit_partResid)
tableFitStatistics(IP3FAFit_partResid,IP3FAFit_noResid)


threeFAfl_partResid = IP3FAFit_partResid$ACE$fl$values
rownames(threeFAfl_partResid) = sfVars
threeFAfl_partResid
fa.diagram(varimax(threeFAfl_partResid)$loadings, cut=.1)

varimax(threeFAfl_partResid)
factor2cluster(varimax(threeFAfl_partResid)$loadings)
IP3FASumm_partResid$CI


# --------------------------------------------------------------------------------------
# Three Factor Independent Pathway Model - A only, partial residual A, correlated factors simple structure model 1
# --------------------------------------------------------------------------------------

promax(threeFAfl_partResid)
factor2cluster(promax(threeFAfl_partResid)$loadings)

nf <- 3

ftrue <- c(T,F,F,F,F,F,F,F,F,F,F,F,
           F,T,F,F,T,T,F,F,F,F,T,F,
           F,F,T,T,F,F,T,T,T,T,F,T)

fload <- sapply(ftrue, function(x) ifelse(x==TRUE, .5, 0))

AsTrue = c(F,F,F,F,F,T,F,F,F,F,F,F)
AsLoad = sapply(AsTrue, function(x) ifelse(x==TRUE, .2, 0))

IP3FAmodel_simpleCorrResid <- mxModel("IndPath3FA_simpleCorrResid",
                                  mxModel("ACE",
                                          # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
                                          mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.3, name="al" ),
                                          mxMatrix( type="Stand", nrow=nf, ncol=nf, free=TRUE, values=0.7, lbound=-1, ubound=1, name="aR" ),
                                          # Algebra for variances components of latent phenotype
                                          mxAlgebra( expression= al %*% aR %*% t(al), name="Al" ),
                                          # Matrix and Algebra for constraint on variance of latent phenotype
                                          mxAlgebra( expression= diag2vec(Al), name="VarLP" ),
                                          mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                                          mxConstraint(VarLP == Unit),
                                          # Matrix f for factor loadings on latent phenotype
                                          mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
                                          mxCI(c("fl")),
                                          
                                          # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                                          mxMatrix( type="Diag", nrow=nv, ncol=nv, free=AsTrue, values=AsLoad, name="as" ),
                                          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
                                          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="es" ),
                                          
                                          # Algebra for residual variance components, with CIs
                                          mxAlgebra( expression= as %*% t(as), name="As" ),
                                          mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
                                          mxAlgebra( expression= es %*% t(es), name="Es" ),
                                          mxAlgebra( expression= As + Cs + Es, name="Vs"), 
                                          #mxCI(c("Vs","As","Cs","Es")),
                                          
                                          # Matrices A, C, and E compute overall variance components
                                          mxAlgebra( expression=(fl %*% al %*% aR %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
                                          mxAlgebra( expression=cs %*% t(cs), name="C" ),
                                          mxAlgebra( expression=es %*% t(es), name="E" ),
                                          # Algebra to compute total variances and standard deviations (diagonal only)
                                          mxAlgebra( expression=A+C+E, name="V" ),
                                          ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                                          # Matrix & Algebra for expected means vector
                                          mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
                                          mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
                                          # Algebra for expected variance/covariance matrix in MZ
                                          mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                                                          cbind(A+C   , A+C+E)), name="expCovMZ" ),
                                          # Algebra for expected variance/covariance matrix in DZ
                                          mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                                          cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
                                  ),
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

IP3FAmodel_simpleCorrResid <- mxOption(IP3FAmodel_simpleCorrResid , "Major iterations", 50000)

IP3FAFit_simpleCorrResid <- mxRun(IP3FAmodel_simpleCorrResid, intervals=FALSE)
IP3FASumm_simpleCorrResid <- summary(IP3FAFit_simpleCorrResid)
IP3FASumm_simpleCorrResid

#tableFitStatistics(multiCholACEFit,c(IP1FAFit_noResid,IP2FAFit_noResid,IP3FAFit_noResid,IP3FAFit_partResid,IP3FAFit_simpleResid))
tableFitStatistics(IP3FAFit_partResid,IP3FAFit_simpleCorrResid)
#tableFitStatistics(IP3FAFit_noResid,IP3FAFit_simpleResid)

tableFitStatistics(multiCholACEFit,c(IP3FAFit,IP3FAFit_partResid,IP3FAFit_simpleCorrResid))


threeFAfl_simpleCorrResid = IP3FAFit_simpleCorrResid$ACE$fl$values
rownames(threeFAfl_simpleCorrResid) = sfVars
threeFAfl_simpleCorrResid
fa.diagram(threeFAfl_simpleCorrResid, simple=FALSE, sort=TRUE)
IP3FASumm_simpleCorrResid$CI

###########################
#     Find E structure    #
###########################

# -----------------------------------------------------------------------
# Single Factor Independent Pathway Model - E only
# -----------------------------------------------------------------------

nf <- 1

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T)
fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP1FEmodel <- mxModel("IndPath1FE",
mxModel("ACE",
        # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
        mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.5, name="el" ),
        # Algebra for variances components of latent phenotype
        mxAlgebra( expression= el %*% t(el), name="El" ),
        # Matrix and Algebra for constraint on variance of latent phenotype
        mxAlgebra( expression= diag2vec(El), name="VarLP" ),
        mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
        mxConstraint(VarLP == Unit),
        # Matrix f for factor loadings on latent phenotype
        mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
        mxCI(c("fl")),
        
        # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.4, name="as" ),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
        mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="es" ),
        
        # Algebra for residual variance components, with CIs
        mxAlgebra( expression= as %*% t(as), name="As" ),
        mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
        mxAlgebra( expression= es %*% t(es), name="Es" ),
        mxAlgebra( expression= As + Cs + Es, name="Vs"), 
        mxCI(c("Vs","As","Cs","Es")),
        
        # Matrices A, C, and E compute overall variance components
        mxAlgebra( expression=(fl %*% el %*% t(el) %*% t(fl)) + (es %*% t(es)), name="E" ),
        mxAlgebra( expression=cs %*% t(cs), name="C" ),
        mxAlgebra( expression=as %*% t(as), name="A" ),
        # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=A+C+E, name="V" ),
        ## Note that the rest of the mxModel statements do not change for bi/multivariate case
        # Matrix & Algebra for expected means vector
        mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
        mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
        # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
        # Algebra for expected variance/covariance matrix in DZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
),
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

IP1FEFit <- mxRun(IP1FEmodel,intervals=FALSE)
IP1FESumm <- summary(IP1FEFit)
IP1FESumm

tableFitStatistics(multiCholACEFit,IP1FEFit)


# -----------------------------------------------------------------------
# Two Factor Independent Pathway Model - E only
# -----------------------------------------------------------------------

nf <- 2

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T,
           F,T,T,T,T,T,T,T,T,T,T,T)

fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP2FEmodel <- mxModel("IndPath2FE",
                      mxModel("ACE",
                              # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
                              mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.5, name="el" ),
                              # Algebra for variances components of latent phenotype
                              mxAlgebra( expression= el %*% t(el), name="El" ),
                              # Matrix and Algebra for constraint on variance of latent phenotype
                              mxAlgebra( expression= diag2vec(El), name="VarLP" ),
                              mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                              mxConstraint(VarLP == Unit),
                              # Matrix f for factor loadings on latent phenotype
                              mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
                              mxCI(c("fl")),
                              
                              # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                              mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.4, name="as" ),
                              mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
                              mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="es" ),
                              
                              # Algebra for residual variance components, with CIs
                              mxAlgebra( expression= as %*% t(as), name="As" ),
                              mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
                              mxAlgebra( expression= es %*% t(es), name="Es" ),
                              mxAlgebra( expression= As + Cs + Es, name="Vs"), 
                              mxCI(c("Vs","As","Cs","Es")),
                              
                              # Matrices A, C, and E compute overall variance components
                              mxAlgebra( expression=(fl %*% el %*% t(el) %*% t(fl)) + (es %*% t(es)), name="E" ),
                              mxAlgebra( expression=cs %*% t(cs), name="C" ),
                              mxAlgebra( expression=as %*% t(as), name="A" ),
                              # Algebra to compute total variances and standard deviations (diagonal only)
                              mxAlgebra( expression=A+C+E, name="V" ),
                              ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                              # Matrix & Algebra for expected means vector
                              mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
                              mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
                              # Algebra for expected variance/covariance matrix in MZ
                              mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                                              cbind(A+C   , A+C+E)), name="expCovMZ" ),
                              # Algebra for expected variance/covariance matrix in DZ
                              mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                              cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
                      ),
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

IP2FEFit <- mxRun(IP2FEmodel,intervals=FALSE)
IP2FESumm <- summary(IP2FEFit)
IP2FESumm

tableFitStatistics(multiCholACEFit,c(IP1FEFit,IP2FEFit))
tableFitStatistics(IP2FEFit,IP1FEFit)

# -----------------------------------------------------------------------
# Three Factor Independent Pathway Model - E only
# -----------------------------------------------------------------------

nf <- 3

ftrue <- c(T,T,T,T,T,T,T,T,T,T,T,T,
           F,T,T,T,T,T,T,T,T,T,T,T,
           F,F,T,T,T,T,T,T,T,T,T,T)

fload <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,
           0,0,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)

IP3FEmodel <- mxModel("IndPath3FE",
                      mxModel("ACE",
                              # Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
                              mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=.5, name="el" ),
                              # Algebra for variances components of latent phenotype
                              mxAlgebra( expression= el %*% t(el), name="El" ),
                              # Matrix and Algebra for constraint on variance of latent phenotype
                              mxAlgebra( expression= diag2vec(El), name="VarLP" ),
                              mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                              mxConstraint(VarLP == Unit),
                              # Matrix f for factor loadings on latent phenotype
                              mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-3, ubound=3),
                              mxCI(c("fl")),
                              
                              # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                              mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.4, name="as" ),
                              mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs" ),
                              mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="es" ),
                              
                              # Algebra for residual variance components, with CIs
                              mxAlgebra( expression= as %*% t(as), name="As" ),
                              mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
                              mxAlgebra( expression= es %*% t(es), name="Es" ),
                              mxAlgebra( expression= As + Cs + Es, name="Vs"), 
                              mxCI(c("Vs","As","Cs","Es")),
                              
                              # Matrices A, C, and E compute overall variance components
                              mxAlgebra( expression=(fl %*% el %*% t(el) %*% t(fl)) + (es %*% t(es)), name="E" ),
                              mxAlgebra( expression=cs %*% t(cs), name="C" ),
                              mxAlgebra( expression=as %*% t(as), name="A" ),
                              # Algebra to compute total variances and standard deviations (diagonal only)
                              mxAlgebra( expression=A+C+E, name="V" ),
                              ## Note that the rest of the mxModel statements do not change for bi/multivariate case
                              # Matrix & Algebra for expected means vector
                              mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
                              mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
                              # Algebra for expected variance/covariance matrix in MZ
                              mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                                              cbind(A+C   , A+C+E)), name="expCovMZ" ),
                              # Algebra for expected variance/covariance matrix in DZ
                              mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                              cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )
                      ),
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

IP3FEFit <- mxRun(IP3FEmodel,intervals=FALSE)
IP3FESumm <- summary(IP3FEFit)
IP3FESumm

tableFitStatistics(multiCholACEFit,c(IP1FEFit,IP2FEFit,IP3FEFit))
tableFitStatistics(IP3FEFit,IP2FEFit)


#############################
# Model Identification Test #
#############################
n    <- 20            # How many permutations of the model do you want to run.
test <- IP3FAmodel_simpleCorrResid     # Model name - what you specify in your mxRun statement
lab <- names(omxGetParameters(test))

resCP1  <- matrix(NA, n, 2*length(lab)+2)
cm   <- 1e10
for (i in 1:n){
  param <- runif(length(lab), -1, 1)    # Range of start values that you want to sample from
  test  <- omxSetParameters(test,
                            labels=lab,
                            values=param
  )
  tr <- mxRun(test, intervals=FALSE)
  resCP1[i,] <- c(param,
                  omxGetParameters(tr),
                  tr@output$Minus2LogLikelihood,
                  tr@output$status[[1]])
  if (tr@output$Minus2LogLikelihood<cm){
    best <- tr
    cm <- tr@output$Minus2LogLikelihood
  }
  print(i)
}

resCP1 <- data.frame(resCP1)
names(resCP1) <- c(
  paste("start", lab, sep=""),
  paste("est",   lab, sep=""),
  "M2LL",
  "status"
)

write.csv(resCP1, "~/Desktop/IP3FAmodel_simpleCorrResid_IDCheck.csv")

