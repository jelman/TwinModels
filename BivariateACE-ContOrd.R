#-----------------------------------------------------------------------#
# Program: Bivariate ACE Model, Continuous and Ordinal Phenotype        #
#  Author: Kelly Spoon                                                  #
#    Date: 08 10 2012                                                   #
#                                                                       #
# Description: Bivariate ACE model					                            #
#		   							                                                    #
#       Input: Raw data file (CSV) with the following columns	 	        #
#  	 IN ORDER: subid, case, twin, zyg, cont phenotype, ord pheno	      #
#	 names of columns do not matter					                              #
#	 *twin variable should be defined A or B			                        #
#	 *zygosity should be defined as 1 for MZ and 2 for DZ		              #
#									                                                      #
#      Output: Bivariate ACE Summary		                                #
#                                                                       #
#	Notes: All comments with *** denote necessary input from user         #
#-----------------------------------------------------------------------#

## Reading in Original Data ##
#-----------------------------------------------------------------------#

#*** Read in data file
	twins <- read.csv("IGEMS_OATS_Delta.csv",header=T)
  names(twins)
  
# Telling R that ordinal variable is ordinal
	twins[,5] <- as.ordered(twins[,5])
##*** Note, this will sort your levels in alphabetical or numeric order, assuming
##*** i.e. 0 < 1 < 2 or A < B < C
##*** If you want to specify another order, it can be done adding a factor statement
##*** i.e. twins[,5] <- factor(twins[,5],levels=c("B","A"),ordered=T)
	levels(twins[,5]) #*** These are listed from least to greatest

## Creating MZ and DZ data sets ##
#-----------------------------------------------------------------------#
# Dividing up twins
	twinA <- twins[twins$twin=="1",]
	twinB <- twins[twins$twin=="2",]
# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("pair","zyg"),all.y=TRUE,suffixes=c("_A","_B"))
# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg==1))
	DZdata <- as.data.frame(subset(newtwins,zyg==2))

	names(newtwins)
	
## Loading Required/Useful Libraries ##
#-----------------------------------------------------------------------#
# Load OpenMX
	require(OpenMx)   #Loads OpenMx
	require(psych)   #Loads Psych package
# Load additional GenEpi tools directly from VCU
	source("GenEpiHelperFunctions.R")

## Defining Variables / Data for use with OpenMX ##
#-----------------------------------------------------------------------#
# Defining number of variables as 2 for OpenMx
	nv <- 2 #Number of phenotypes
	ntv <- nv*2
# Making Sparse Data Frame with only phenotype 
	mzData <- MZdata[,c(8,5,15,12)]
	dzData <- DZdata[,c(8,5,15,12)]
# Define Variable to use in OpenMX
	selVars <- c(names(mzData))
	Vars <- 'var'

# Defining start values for ordinal variable, data driven
	nlev <- length(levels(twins[,5])) # Defining number of levels 
	nthr <- nlev - 1 # Defining number of thresholds
# Start values for Thresholds by determining proportion of values and converting proportions 
# to values from standard normal distribution
	thresh <- c()
	all <- sum(table(twins[,5]))
	val <- 0
	for (i in 1:nthr){
		group <- table(twins[,5])[i]
		val <- val + group
		thresh[i] <- round(qnorm(val/all),2)
	}
# Names and labels for threshold
	th <- paste("th",c(1:nthr),sep="")

# Defining start values for continuous variable, may require tweaking
	avals <- c(.4,.3,.01)
	cvals <- c(.1,.1,.01)
	evals <- c(.4,.2,.2)

# Fit Bivariate ACE Model with Raw Data Input
# -----------------------------------------------------------------------#

bivACEModel <- mxModel("ACE",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=avals, name="a", lbound = -30, ubound = 30),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cvals, name="c", lbound = -30, ubound = 30),
        mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=evals, name="e", lbound = -30, ubound = 30),
    # Matrices A, C, and E compute variance components
        mxAlgebra( expression=a %*% t(a), name="A"),
        mxAlgebra( expression=c %*% t(c), name="C"),
        mxAlgebra( expression=e %*% t(e), name="E"),
    # Algebra to compute total variances and standard deviations (diagonal only)
        mxAlgebra( expression=A+C+E, name="V" ),
        mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
        mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
    
        mxAlgebra( expression= V[2,2], name="V22" ),
        mxConstraint(V22 == 1),
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
    # Matrix & Algebra for expected means vector, set as 0 for ordinal variable
        mxMatrix( type="Full", nrow=1, ncol=nv, free=c(TRUE,FALSE), values=0, labels=c("mean1","mean2"), name="Mean"),
        mxAlgebra( expression= cbind(Mean,Mean), name="expMean"),
    ## This is the line that allows OpenMx to deal with the ordinal data, it treats it as a normal variable, centered at 0	
        mxMatrix(type="Full", nrow=nthr, ncol=nv, free=TRUE, values=thresh, name="expThre", labels=th,dimnames=list(th,selVars[c(2,4)])),
    # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
                                        cbind(A+C   , A+C+E)), name="expCovMZ" ),
    # Algebra for expected variance/covariance matrix in DZ
        mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ),
    mxModel("MZ",
        mxData( observed=mzData, type="raw" ),
        mxFIMLObjective( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selVars, thresholds="ACE.expThre", threshnames=selVars[c(2,4)])
    ),
    mxModel("DZ", 
        mxData( observed=dzData, type="raw" ),
        mxFIMLObjective( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selVars, thresholds="ACE.expThre", threshnames=selVars[c(2,4)]) 
    ),
    mxAlgebra( expression=MZ.objective + DZ.objective, name="neg2sumll" ),
    mxAlgebraObjective("neg2sumll")
)
#bivACEFit <- mxTryHard(bivACEModel,intervals=TRUE, extraTries=50)
bivACEFit <- mxRun(bivACEModel,intervals=TRUE)
bivACESumm <- summary(bivACEFit)
bivACESumm

##In your output, the confidence intervals portion will likely be of the most interest to you
# ACE.stndVCA[1,1] will give you the estimate for the genetic portion for cont test measure, while
# ACE.stndVCA[2,2] will give you genetic portion for your ordinal measure
# ACE.Acorr will give you the genetic correlation between the two measures

#############################
# Model Identification Test #
#############################
n    <- 50
test <- bivACEModel 
lab <- names(omxGetParameters(test))

resCP1  <- matrix(NA, n, 2*length(lab)+2)
cm   <- 1e10
for (i in 1:n){
  param <- runif(length(lab), -10, 10)
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

