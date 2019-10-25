#-----------------------------------------------------------------------#
# Program: Univariate Saturated and ACE Models, Ordinal Data	        	#
#  Author: Kelly Spoon                                                  #
#    Date: 08 16 2011                                                   #
#                                                                       #
# Description: Complete univariate analyses of Saturated and ACE models #
#		   for ordinal data							                                    #
#-----------------------------------------------------------------------#

#*** Read in data file
	twins <- read.csv("VETSA_Alcohol_Updated.csv",header=T)
	names(twins)
	
#*** Quick check that data file is properly read in
	dim(twins)[1] #*** should be number of subjects 
	names(twins)[1:4] #*** should be subid, case, twin, zygosity
	names(twins)[25] #*** should be phenotype
	
#*** Defining number of variables
	nv <- 1
	ntv <- 2*nv

## Creating Starting Values
	twins[,26] <- as.ordered(twins[,26]) # Defining variable as ordinal in R
	
##*** Note, this will sort your levels in alphabetical or numeric order, assuming
##*** i.e. 0 < 1 < 2 or A < B < C
##*** If you want to specify another order, it can be done adding a factor statement
##*** i.e. twins[,5] <- factor(twins[,5],levels=c("B","A"),ordered=T)
	nlev <- length(levels(twins[,26])) # Defining number of levels 
	nthr <- nlev - 1 # Defining number of thresholds
	levels(twins[,26]) #*** These are listed from least to greatest
	thresh <- c()
	val <- 0
	for (i in 1:nthr){
		val <- val + table(twins[,26])[i]/sum(table(twins[,26]))
		thresh[i] <- round(qnorm(val),2)
	}
	thlab <- paste("th",c(1:nthr),c("_A","_B"),sep="")
	thname <- paste("th",c(1:nthr),sep="")
	thlab <- paste("th",c(1:nthr),sep="")
	thlabd <- paste(thlab,"d",sep="")
	thnamed <- paste(thname,"d",sep="")	

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]
# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg2019"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg2019==1))
	DZdata <- as.data.frame(subset(newtwins,zyg2019==2))

	names(MZdata)
	
# Defining Variables of Interest
	selVars <- names(MZdata)[c(26,52)]
	selVars
	
	nthr
	
## Loading Required/Useful Libraries ##

# Load OpenMX
	require(OpenMx)   #Loads OpenMx
	require(psych)   #Loads Psych package
# Load additional GenEpi tools directly from VCU
	source("GenEpiHelperFunctions.R")


# Specify and Run Saturated Model (Tetrachoric correlations) with RawData and Matrix-style Input
# -----------------------------------------------------------------------
	
# Set default Mx optimizer to NPSOL
#mxOption(NULL, "Default optimizer", "NPSOL")  
	
twinSatModel <- mxModel("twinSat",
mxModel("MZ",
  # Matrix & Algebra for expected means vector (SND), Thresholds and correlation
	mxMatrix( type="Zero", nrow=1, ncol=nv, name="M" ),
	mxAlgebra( expression= cbind(M,M), name="expMeanMZ" ),
	mxMatrix(type="Full", nrow=nthr, ncol=ntv, free=TRUE, values=thresh, name="expThreMZ", dimnames=list(thname,selVars) ),
	mxMatrix(type="Stand", nrow=ntv, ncol=ntv, free=TRUE, values=.3, lbound=-1, ubound=1, name="expCorMZ"), 
  mxCI(c("expCorMZ")),
 	mxData(MZdata, type="raw"),
 	mxFIMLObjective( covariance="expCorMZ", means="expMeanMZ", dimnames=selVars, thresholds="expThreMZ" )),
 	
mxModel("DZ",
  # Matrix & Algebra for expected means vector (SND), Thresholds and correlation
	mxMatrix( type="Zero", nrow=1, ncol=nv, name="M" ),
	mxAlgebra( expression= cbind(M,M), name="expMeanDZ" ),
	mxMatrix( type="Full", nrow=nthr, ncol=ntv, free=TRUE, values=thresh, name="expThreDZ", dimnames=list(thnamed,selVars)),
	mxMatrix(type="Stand", nrow=ntv, ncol=ntv, free=TRUE, values=.03, lbound=-1, ubound=1, name="expCorDZ"), 
  mxCI(c("expCorDZ")),
 	mxData(DZdata, type="raw"),
 	mxFIMLObjective( covariance="expCorDZ", means="expMeanDZ", dimnames=selVars, thresholds="expThreDZ" )),


 	mxAlgebra( MZ.objective + DZ.objective , name="min2sumll" ),
 	mxAlgebraObjective("min2sumll") 
 
 ) 
twinSatFit <- mxTryHard(twinSatModel,intervals=TRUE, extraTries=30)
twinSatFit <- mxRun(twinSatModel,intervals=TRUE)
twinSatSumm <- summary(twinSatFit)
twinSatSumm


# Fit ACE Model with RawData and Matrices Input, ONE overall Threshold
# ---------------------------------------------------------------------

# Set default Mx optimizer to NPSOL
#mxOption(NULL, "Default optimizer", "NPSOL")  

univACEOrdModel <- mxModel("univACEOrd",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, name="a", lbound = -1, ubound = 1 ),
        mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.1, name="c", lbound = -1, ubound = 1 ),
        mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, name="e", lbound = -1, ubound = 1 ),

    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=a %*% t(a), name="A" ),
        mxAlgebra( expression=c %*% t(c), name="C" ),
        mxAlgebra( expression=e %*% t(e), name="E" ),
    
    # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=A+C+E, name="V" ),
    
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(sqrt(I*V)), name="sd"),
    
	  mxAlgebra( expression=cbind(A/V,C/V,E/V),name="stndVCs"),
    
    # Calculate 95% CIs here
	  mxCI(c("stndVCs")),
    
    # Constraint on variance of ordinal variables
        mxConstraint(V == I, name="Var"),
    
    # Matrix & Algebra for expected means vector
        mxMatrix( type="Zero", nrow=1, ncol=nv, name="M" ),
        mxAlgebra( expression= cbind(M,M), name="expMean" ),
    
        mxMatrix( type="Full", nrow=nthr, ncol=nv, free=TRUE, values=thresh, name="Thre" ),
        mxAlgebra( expression= cbind(Thre,Thre), name="expThre" ),
    
  
    # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
    # Algebra for expected variance/covariance matrix in DZ, note use of 0.5, converted to 1*1 matrix
        mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ),
   
    mxModel("MZ",
        mxData( observed=MZdata, type="raw" ),
        mxFIMLObjective( covariance="univACEOrd.expCovMZ", means="univACEOrd.expMean", dimnames=selVars, thresholds="univACEOrd.expThre" )),
    mxModel("DZ", 
        mxData( observed=DZdata, type="raw" ),
        mxFIMLObjective( covariance="univACEOrd.expCovDZ", means="univACEOrd.expMean", dimnames=selVars, thresholds="univACEOrd.expThre" )),

    mxAlgebra( expression=MZ.objective + DZ.objective, name="min2sumll" ),
    mxAlgebraObjective("min2sumll")
)

univACEOrdFit <- mxTryHard(univACEOrdModel,intervals=TRUE, extraTries=30)

univACEOrdFit <- mxRun(univACEOrdModel,intervals=T)
univACEOrdSumm <- summary(univACEOrdFit)
univACEOrdSumm

tableFitStatistics(twinSatFit,univACEOrdFit)

univACEOrdFit$univACEOrd.A@result
univACEOrdFit$univACEOrd.C@result
univACEOrdFit$univACEOrd.E@result


# Fit AE Model with RawData and Matrices Input, ONE overall Threshold
# ---------------------------------------------------------------------

# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL")  

univAEOrdModel <- mxModel("univAEOrd",
   # Matrices a, c, and e to store a, c, and e path coefficients
   mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label="a11", name="a", lbound = -1, ubound = 1 ),
   mxMatrix( type="Full", nrow=nv, ncol=nv, free=FALSE, values=.0, label="c11", name="c", lbound = -1, ubound = 1 ),
   mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label="e11", name="e", lbound = -1, ubound = 1 ),
   # Matrices A, C, and E compute variance components
   mxAlgebra( expression=a %*% t(a), name="A" ),
   mxAlgebra( expression=c %*% t(c), name="C" ),
   mxAlgebra( expression=e %*% t(e), name="E" ),
   # Algebra to compute total variances and standard deviations (diagonal only)
   mxAlgebra( expression=A+C+E, name="V" ),
   mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
   mxAlgebra( expression=solve(sqrt(I*V)), name="sd"),
   mxAlgebra( expression=cbind(A/VP,C/VP,E/VP),name="stndVCs"),
   # Calculate 95% CIs here
   mxAlgebra(A+C+E,name="VP"),
   mxCI(c("stndVCs")),
   # Constraint on variance of ordinal variables
   mxConstraint(V == I, name="Var1"),
   # Matrix & Algebra for expected means vector
   mxMatrix( type="Zero", nrow=1, ncol=nv, name="M" ),
   mxAlgebra( expression= cbind(M,M), name="expMean" ),
   mxMatrix( type="Full", nrow=nthr, ncol=ntv, free=TRUE, values=thresh, labels=thname, name="expThre" ),
   # Algebra for expected variance/covariance matrix in MZ
   mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                   cbind(A+C   , A+C+E)), name="expCovMZ" ),
   # Algebra for expected variance/covariance matrix in DZ, note use of 0.5, converted to 1*1 matrix
   mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                   cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ),
   mxModel("MZ",
           mxData( observed=MZdata, type="raw" ),
           mxFIMLObjective( covariance="univAEOrd.expCovMZ", means="univAEOrd.expMean", dimnames=selVars, thresholds="univAEOrd.expThre" )),
   mxModel("DZ", 
           mxData( observed=DZdata, type="raw" ),
           mxFIMLObjective( covariance="univAEOrd.expCovDZ", means="univAEOrd.expMean", dimnames=selVars, thresholds="univAEOrd.expThre" )),
   mxAlgebra( expression=MZ.objective + DZ.objective, name="min2sumll" ),
   mxAlgebraObjective("min2sumll")
)

univAEOrdFit <- mxRun(univAEOrdModel,intervals=T)
univAEOrdSumm <- summary(univAEOrdFit)
univAEOrdSumm

tableFitStatistics(twinSatFit,c(univACEOrdFit,univAEOrdFit))


# Generate ACE Output
# -----------------------------------------------------------------------
parameterSpecifications(univACEOrdFit)
expectedMeansCovariances(univACEOrdFit)
tableFitStatistics(univACEOrdFit)




# Fit ACE Model with RawData and Matrices Input, ONE overall Threshold
# ---------------------------------------------------------------------

# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL")  

univADEOrdModel <- mxModel("univADEOrd",
   # Matrices a, c, and e to store a, c, and e path coefficients
   mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.0, label="a11", name="a", lbound = -1, ubound = 1 ),
   mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.5, label="d11", name="d", lbound = -1, ubound = 1 ),
   mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.6, label="e11", name="e", lbound = -1, ubound = 1 ),
   # Matrices A, C, and E compute variance components
   mxAlgebra( expression=a %*% t(a), name="A" ),
   mxAlgebra( expression=d %*% t(d), name="D" ),
   mxAlgebra( expression=e %*% t(e), name="E" ),
   # Algebra to compute total variances and standard deviations (diagonal only)
   mxAlgebra( expression=A+D+E, name="V" ),
   mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
   mxAlgebra( expression=solve(sqrt(I*V)), name="sd"),
   mxAlgebra( expression=cbind(A/VP,D/VP,E/VP),name="stndVCs"),
   # Calculate 95% CIs here
   mxAlgebra(A+D+E,name="VP"),
   mxCI(c("stndVCs")),
   # Constraint on variance of ordinal variables
   mxConstraint(V == I, name="Var1"),
   # Matrix & Algebra for expected means vector
   mxMatrix( type="Zero", nrow=1, ncol=nv, name="M" ),
   mxAlgebra( expression= cbind(M,M), name="expMean" ),
   mxMatrix( type="Full", nrow=nthr, ncol=ntv, free=TRUE, values=thresh, labels=thname, name="expThre" ),
   # Algebra for expected variance/covariance matrix in MZ
   mxAlgebra( expression= rbind  ( cbind(A+D+E , A+D),
                                   cbind(A+D   , A+D+E)), name="expCovMZ" ),
   # Algebra for expected variance/covariance matrix in DZ, note use of 0.5, converted to 1*1 matrix
   mxAlgebra( expression= rbind  ( cbind(A+D+E     , 0.5%x%A+0.25%x%D),
                                   cbind(0.5%x%A+0.25%x%D , A+D+E)),  name="expCovDZ" ),
   mxModel("MZ",
           mxData( observed=MZdata, type="raw" ),
           mxFIMLObjective( covariance="univADEOrd.expCovMZ", means="univADEOrd.expMean", dimnames=selVars, thresholds="univADEOrd.expThre" )),
   mxModel("DZ", 
           mxData( observed=DZdata, type="raw" ),
           mxFIMLObjective( covariance="univADEOrd.expCovDZ", means="univADEOrd.expMean", dimnames=selVars, thresholds="univADEOrd.expThre" )),
   mxAlgebra( expression=MZ.objective + DZ.objective, name="min2sumll" ),
   mxAlgebraObjective("min2sumll")
)

univADEOrdFit <- mxRun(univADEOrdModel,intervals=T)
univADEOrdSumm <- summary(univADEOrdFit)
univADEOrdSumm
LL_ADE <- mxEval(objective, univADEOrdFit)
LL_ADE

tableFitStatistics(twinSatFit,c(univACEOrdFit,univADEOrdFit))


# Fit AE model
# ---------------------------------------------------------------------

univAEOrdModel <- mxRename(univADEOrdModel, "univAEOrd")
univAEOrdModel$univAEOrd.d  <-  mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, values=0, label="d11", name="d" ) # drop c at 0
univAEOrdFit <- mxRun(univAEOrdModel,intervals=T)
univAEOrdSumm <- summary(univAEOrdFit)
univAEOrdSumm
LL_AE <- mxEval(objective, univAEOrdFit)
LL_AE


tableFitStatistics(twinSatFit,c(univACEOrdFit,univADEOrdFit,univAEOrdFit))


# Fit CE model
# --------------------------------------------------------------------
univCEOrdModel <- mxRename(univACEOrdModel, "univCEOrd")
univCEOrdModel$univCEOrd.a  <-  mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, values=0, label="a11", name="a" ) # drop a at 0
univCEOrdFit <- mxRun(univCEOrdModel)
univCEOrdSumm <- summary(univCEOrdFit)
univCEOrdSumm

# Fit E model
# ---------------------------------------------------------------------
univEOrdModel <- mxRename(univAEOrdFit, "univEOrd")
univEOrdModel$univEOrd.a  <-  mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, values=0, label="a11", name="a" ) # drop a at 0
univEOrdFit <- mxRun(univEOrdModel)
univEOrdSumm <- summary(univEOrdFit)
univEOrdSumm

tableFitStatistics(twinSatFit,c(univACEOrdFit,univADEOrdFit,univAEOrdFit,univEOrdFit))

# Print Comparative Fit Statistics ACE models
# ---------------------------------------------------------------------
univACEOrdNested <- list(univAEOrdFit, univCEOrdFit, univEOrdFit)
tableFitStatistics(univACEOrdFit,univACEOrdNested)


