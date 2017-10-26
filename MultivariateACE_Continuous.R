#-----------------------------------------------------------------------#
# Program: Multivariate ACE Model                           			      #
# Contents: Saturated Model, Cholesky Model, and Common Pathways       	#
# Created: December 5, 2016			  					                            #
#												                                                #
# Data Type: Continuous									                                #
#												                                                #
#-----------------------------------------------------------------------#

### Reading in Original Data ###

# Set working directory - location of data file and desired output
setwd("C:/Users/...")

## Loading Required/Useful Libraries ##

# Load OpenMX
require(OpenMx)   #Loads OpenMx
require(psych)    #Loads Psych package

# Load additional GenEpi tools
source("GenEpiHelperFunctions.R")

# Read in data file
twins <- read.csv("V1Memory.csv",header=T)

# Quick check that data file is properly read in
	names(twins) 	  # list variables in data file
	dim(twins)[1] 	# number of subjects
	dim(twins)[2] 	# number of variables in the file


## Preparing the data for OpenMx ##

# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("CASE","ZYG14"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
	names(newtwins)   # this is what your paired dataset now looks like
	dim(newtwins)[1] 	# number of pairs (complete and incomplete)
	dim(newtwins)[2] 	# number of variables in the file

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,ZYG14==1))
	DZdata <- as.data.frame(subset(newtwins,ZYG14==2))
	dim(MZdata)[1]
	dim(DZdata)[1]


# Defining number of phenotype for the model
nv <- 6		# Number of phenotypes 
ntv <- nv*2		# Number of variables (phenotypes x 2)

names(newtwins)

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,ZYG14==1,c(5,6,7,8,9,10,13,14,15,16,17,18)))		
	DZdata <- as.data.frame(subset(newtwins,ZYG14==2,c(5,6,7,8,9,10,13,14,15,16,17,18)))		

selvars <- names(MZdata)
selvars

# Print Descriptive Statistics
describe(MZdata)		 
describe(DZdata)

cor(MZdata,use="complete")
cor(DZdata,use="complete")



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
cholASVnv <- rep(.6,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholCSVnv <- rep(.2,nv*(nv+1)/2)  #Establishes starting values based on the data 
cholESVnv <- rep(.3,nv*(nv+1)/2)  #Establishes starting values based on the data 

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
	  mxCI(c("stndVCA","stndVCC","stndVCE")),

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

tableFitStatistics(multiSatFit, multiCholACEFit)


# --------------------------------------------------------------#

## Two Factor Common Pathways Model ##

# --------------------------------------------------------------#

# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL")  

selvars

nf <- 2

ftrue <- c(T,T,T,T,T,T,
           F,T,F,T,F,T)

fload <- c(.4,.4,.4,.4,.4,.4,
           0,.4,0,.4,0,.4)

restrue <- c(T,T,F,F,F,F,
             T,F,F,F,F,
             T,T,F,F,
             T,F,F,
             T,T,
             T)

resload <- c(.1,.1,0,0,0,0,
             .1,0,0,0,0,
             .1,.1,0,0,
             .1,0,0,
             .1,.1,
             .1)

ModelCP1 <- mxModel("BiFactorCP",
  mxModel("ACE",
          # Matrices to store a, c, and e path coefficients for latent phenotype(s)
          mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=c(.4,.1), name="al", lbound=-3, ubound=3),
          mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=c(.1,.01), name="cl", lbound=-3, ubound=3 ),
          mxMatrix( type="Diag", nrow=nf, ncol=nf, free=TRUE, values=c(.5,.1), name="el", lbound=-3, ubound=3 ),
          
          # Algebra to calculate Variances for latent phenotype(s)
          mxAlgebra( al %*% t(al), name="Al"),
          mxAlgebra( cl %*% t(cl), name="Cl"),
          mxAlgebra( el %*% t(el), name="El"),
          mxAlgebra( Al + Cl + El, name="Vl"),
          mxCI(c("Al","Cl","El")),
          
          # Constraining the Variance to 1.0 
          mxAlgebra( diag2vec(Vl), name="latvar" ),
          mxMatrix( type="Full", nrow=nf, ncol=1, free=FALSE, values=1, name="unit"),
          mxConstraint( latvar == unit, name="varcons"),
          
          # Matrix f for factor loadings on latent phenotype
          mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-1, ubound=1),
          mxCI("fl"),
          
          # Matrices to store a, c, and e residual components for specific phenotypes
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=restrue, values=resload, name="as", lbound=-1, ubound=1 ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=restrue, values=resload, name="cs", lbound=-1, ubound=1 ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=restrue, values=resload, name="es", lbound=-1, ubound=1 ),

          # Algebra for residual variance components, with CIs
          mxAlgebra( expression= as %*% t(as), name="As" ),
          mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
          mxAlgebra( expression= es %*% t(es), name="Es" ),
          #mxCI(c("As","Cs","Es")),

          # Adding it all together
          mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
          mxAlgebra( expression=(fl %*% cl %*% t(cl) %*% t(fl)) + (cs %*% t(cs)), name="C" ),
          mxAlgebra( expression=(fl %*% el %*% t(el) %*% t(fl)) + (es %*% t(es)), name="E" ),
          
          # Algebra to compute total variances and standard deviations (diagonal only)
          mxAlgebra( expression=A+C+E, name="V" ),

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

ModelCP1 <- mxOption(ModelCP1, "Major iterations", 40000)

FitCP1 <- mxRun(ModelCP1,intervals=TRUE)
SummCP1 <- summary(FitCP1)
SummCP1 

tableFitStatistics(multiCholACEFit,FitCP1 )


#############################
# Model Identification Test #
#############################
n    <- 50            # How many permutations of the model do you want to run.
test <- ModelCP1     # Model name - what you specify in your mxRun statement
lab <- names(omxGetParameters(test))

resCP1  <- matrix(NA, n, 2*length(lab)+2)
cm   <- 1e10
for (i in 1:n){
  param <- runif(length(lab), -1, 1)    # Range of start values that you want to sample from
  test  <- omxSetParameters(test,
                            labels=lab,
                            values=param
  )
  tr <- mxRun(test)
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

write.csv(resCP1, "IDCheck.csv")

##------------------------##

nf <- 1

fload <- c(.4,.4,.4)

ftrue <- c(TRUE,TRUE,TRUE)

ModelCP1r <- mxModel("OneFactorCP",
  mxModel("ACE",
          # Matrices to store a, c, and e path coefficients for latent phenotype(s)
          mxMatrix( type="Lower", nrow=nf, ncol=nf, free=TRUE, values=.4, name="al", lbound=-1, ubound=1),
          mxMatrix( type="Lower", nrow=nf, ncol=nf, free=TRUE, values=.01, name="cl", lbound=-1, ubound=1 ),
          mxMatrix( type="Lower", nrow=nf, ncol=nf, free=TRUE, values=.5, name="el", lbound=-1, ubound=1 ),
          
          # Algebra to calculate Variances for latent phenotype(s)
          mxAlgebra( al %*% t(al), name="Al"),
          mxAlgebra( cl %*% t(cl), name="Cl"),
          mxAlgebra( el %*% t(el), name="El"),
          mxAlgebra( Al + Cl + El, name="Vl"),
          mxCI(c("Al","Cl","El")),
          
          # Constraining the Variance to 1.0 
          mxAlgebra( diag2vec(Vl), name="latvar" ),
          mxMatrix( type="Full", nrow=nf, ncol=1, free=FALSE, values=1, name="unit"),
          mxConstraint( latvar == unit, name="varcons"),
          
          # Matrix f for factor loadings on latent phenotype
          mxMatrix( type="Full", nrow=nv, ncol=nf, free=ftrue, values=fload, name="fl", lbound=-1, ubound=1),
          mxCI("fl"),
          
          # Matrices to store a, c, and e residual components for specific phenotypes
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(T,F,F,T,T,T), values=c(.1,0,0,.1,.1,.1), name="as", lbound=-1, ubound=1 ),
          mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.1, name="cs", lbound=-1, ubound=1 ),
          mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE, values=.2, name="es", lbound=-1, ubound=1 ),
          
          # Algebra for residual variance components, with CIs
          mxAlgebra( expression= as %*% t(as), name="As" ),
          mxAlgebra( expression= cs %*% t(cs), name="Cs" ),
          mxAlgebra( expression= es %*% t(es), name="Es" ),
          mxCI(c("As","Cs","Es")),
          
          # Adding it all together
          mxAlgebra( expression=(fl %*% al %*% t(al) %*% t(fl)) + (as %*% t(as)), name="A" ),
          mxAlgebra( expression=(fl %*% cl %*% t(cl) %*% t(fl)) + (cs %*% t(cs)), name="C" ),
          mxAlgebra( expression=(fl %*% el %*% t(el) %*% t(fl)) + (es %*% t(es)), name="E" ),
          
          # Algebra to compute total variances and standard deviations (diagonal only)
          mxAlgebra( expression=A+C+E, name="V" ),
          
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

ModelCP1r <- mxOption(ModelCP1r, "Major iterations", 40000)

FitCP1r <- mxRun(ModelCP1r,intervals=FALSE)
SummCP1r <- summary(FitCP1r)
SummCP1r 


tableFitStatistics(multiCholACEFit,c(FitCP1,FitCP1r ))
tableFitStatistics(FitCP1r,FitCP1)


#############################
# Model Identification Test #
#############################
n    <- 50            # How many permutations of the model do you want to run.
test <- ModelCP1r     # Model name - what you specify in your mxRun statement
lab <- names(omxGetParameters(test))

resCP1  <- matrix(NA, n, 2*length(lab)+2)
cm   <- 1e10
for (i in 1:n){
  param <- runif(length(lab), -1, 1)    # Range of start values that you want to sample from
  test  <- omxSetParameters(test,
                            labels=lab,
                            values=param
  )
  tr <- mxRun(test)
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

write.csv(resCP1, "IDCheck2.csv")

##------------------------##