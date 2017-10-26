#-----------------------------------------------------------------------#
# Program: Univariate Saturated and ACE Model                           #
#  Author: Kelly Spoon                                                  #
#    Date: 07 07 2011                                                   #
#Modified: 08 31 2015 			                                    #
#                                                                       #
# Description: Complete univariate analysis (Saturated, ACE models,  	#
#		   & Nested Models							#
#									                  #
#       Input: Raw data file (CSV) with the following columns	 	#
#  		   IN ORDER: subid, case, twin, zygosity, phenotype		#
#		   names of columns do not matter					#
#		   *twin variable should be defined A or B			#
#		   *zygosity should be defined as 1 for MZ and 2 for DZ	#
#												#
#	   Data: Continuous								#
#												#
#      Output: All univariate statistics                                #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user   #
#-----------------------------------------------------------------------#

### Reading in Original Data ###

#*** Set working directory - location of data file and desired output
	setwd("C:/Users/mspanizzon/Desktop/Univariate ACE Model")

#*** Read in data file
	twins <- read.csv("V1_HeightAFQT.csv",header=T)

#*** Quick check that data file is properly read in
	dim(twins)[1] 	# should be number of subjects 
	dim(twins)[2] 	# should be 5
	names(twins) 	# should be subid, case, twin, zygocity, phenotype

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg10"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
	names(newtwins) # this is what your paired dataset now looks like

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg10==1))
	DZdata <- as.data.frame(subset(newtwins,zyg10==2))

## Loading Required/Useful Libraries ##

# Load OpenMX
	require(OpenMx)   #Loads OpenMx
	require(psych)   #Loads Psych package

# Load additional GenEpi tools directly from VCU
	source("http://www.vipbg.vcu.edu/~vipbg/Tc24/GenEpiHelperFunctions.R")

## Defining Variables / Data for use with OpenMX ##

	names(MZdata)

# Making Sparse Data Frame with only phenotype of interest
	MZdata <- MZdata[,c(8,14)] 	#*** This will grab the 5th and 14th columns of your data frame
	DZdata <- DZdata[,c(8,14)]

# Getting Number of Twins Used in Analysis
	cpMZ <- sum(complete.cases(MZdata)) # Number of complete MZ pairs
	sMZ <- dim(MZdata)[1]-cpMZ 		# Number of MZ singletons
	cpDZ <- sum(complete.cases(DZdata)) # Number of complete DZ pairs
	sDZ <- dim(DZdata)[1]-cpDZ 		# Number of DZ singletons

# Define Variable to use in OpenMX
	selvars <- c(names(MZdata))		# Specifying the variable names (MZ or DZ doesn't matter)
	Vars <- 'var'
	nv <- 1 					# Defines the number of phenotypes
	ntv <- nv*2


# Covariance between MZ twins 
	MZcov<-cov(MZdata,use="complete")
	MZcov

# Covariance between DZ twins 
	DZcov<-cov(DZdata,use="complete")
	DZcov


# Correlation between MZ twins 
	MZcor<-cor(MZdata,use="complete")
	MZcor

	cor.test(MZdata[,1],MZdata[,2]) # Tests the significance of the correlation

# Correlation between DZ twins 
	DZcor<-cor(DZdata,use="complete")
	DZcor

	cor.test(DZdata[,1],DZdata[,2]) # Tests the significance of the correlation

#-----------------------------------------------------------------------#

### Begin Univariate Analysis ###

#-----------------------------------------------------------------------#

## Fit Univariate Saturated Model ##

univTwinSatModel <- mxModel("univTwinSat",
    mxModel("MZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=.6, name="CholMZ" ),  #*** Starting value is 10 for all matrix elements
        mxAlgebra( expression=CholMZ %*% t(CholMZ), name="expCovMZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=0, name="expMeanMZ" ), #*** Starting value is 100 for all matrix elements (means)
        mxData( observed=MZdata, type="raw" ),
	  mxExpectationNormal( covariance="MZ.expCovMZ", means="MZ.expMeanMZ", dimnames=selvars),  
	  mxFitFunctionML(),

    # Algebra's needed for equality constraints
        mxAlgebra( expression=expMeanMZ[1,1:nv], name="expMeanMZt1"),
        mxAlgebra( expression=expMeanMZ[1,(nv+1):ntv], name="expMeanMZt2"),
        mxAlgebra( expression=t(diag2vec(expCovMZ)), name="expVarMZ"),
        mxAlgebra( expression=expVarMZ[1,1:nv], name="expVarMZt1"),
        mxAlgebra( expression=expVarMZ[1,(nv+1):ntv], name="expVarMZt2")
    ),

    mxModel("DZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=.6, name="CholDZ" ),  #*** Starting value is 10 for all matrix elements
        mxAlgebra( expression=CholDZ %*% t(CholDZ), name="expCovDZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=0, name="expMeanDZ" ),	  #*** Starting value is 100 for all matrix elements (means)
        mxData( observed=DZdata, type="raw" ),
	  mxExpectationNormal( covariance="DZ.expCovDZ", means="DZ.expMeanDZ", dimnames=selvars),  
	  mxFitFunctionML(),

    # Algebra's needed for equality constraints
        mxAlgebra( expression=expMeanDZ[1,1:nv], name="expMeanDZt1"),
        mxAlgebra( expression=expMeanDZ[1,(nv+1):ntv], name="expMeanDZt2"),
        mxAlgebra( expression=t(diag2vec(expCovDZ)), name="expVarDZ"),
        mxAlgebra( expression=expVarDZ[1,1:nv], name="expVarDZt1"),
        mxAlgebra( expression=expVarDZ[1,(nv+1):ntv], name="expVarDZt2")
    ),

    mxAlgebra( MZ.objective + DZ.objective, name="neg2sumll" ),
    mxFitFunctionAlgebra("neg2sumll")	# Calculates the -2LL

)
	univTwinSatFit <- mxRun(univTwinSatModel)
	univTwinSatSumm <- summary(univTwinSatFit)
	univTwinSatSumm

univTwinSatFit$MZ.CholMZ
univTwinSatFit$MZ.expCovMZ
univTwinSatFit$MZ.expMeanMZ
univTwinSatFit$DZ.CholDZ
univTwinSatFit$DZ.expCovDZ
univTwinSatFit$DZ.expMeanDZ

#-----------------------------------------------------------------------#

## Univariate ACE Model ##

univACEModel <- mxModel("univACE",
	mxModel("ACE",
        	mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.7, name="a", lbound=-5, ubound=5 ),
        	mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.3, name="c", lbound=-5, ubound=5 ),
       	mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.4, name="e", lbound=-5, ubound=5 ),

		# Matrices A, C, and E compute variance components
		mxAlgebra( expression=a %*% t(a), name="A" ),  # Variance Components
		mxAlgebra( expression=c %*% t(c), name="C" ),  # Variance Components
		mxAlgebra( expression=e %*% t(e), name="E" ),  # Variance Components
		mxAlgebra( expression=A+C+E, name="V" ),	     # Calculates the Phenotypic Variance
		mxAlgebra( expression=cbind(A/V,C/V,E/V),name="stndVCs"), # Standardized Variance Components

     		# 95% Confidence Intervals
		mxCI(c("A","C","E","stndVCs")),	# Calculates the 95% confidence intervals for "stndVCs"

		# Means
		mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=0, label="mean", name="expMean" ),

    		# Algebra for expected variance/covariance matrix in MZ
   		mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
							  cbind(A+C   , A+C+E)), name="expCovMZ" ),

    		# Algebra for expected variance/covariance matrix in DZ 
   		mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
 						        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )),

	mxModel("MZ",
	    		mxData( observed=MZdata, type="raw" ),
	            mxExpectationNormal( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selvars),  
	  		mxFitFunctionML()),

	mxModel("DZ",
	    		mxData( observed=DZdata, type="raw" ),
	            mxExpectationNormal( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selvars),  
	  		mxFitFunctionML()),

    		mxAlgebra( expression=MZ.objective + DZ.objective, name="twin" ),
	      mxFitFunctionAlgebra("twin")	# Calculates the -2LL
	)

	univACEFit <- mxRun(univACEModel, intervals=TRUE) 
	univACESumm <- summary(univACEFit)
	univACESumm

tableFitStatistics(univTwinSatFit , univACEFit )

#-----------------------------------------------------------------------#

# Fitting AE model #
univAEModel <- mxModel("univAE",
	mxModel("ACE",
        	mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.7, name="a", lbound=-1, ubound=1 ),
        	mxMatrix( type="Full", nrow=nv, ncol=nv, free=FALSE, values=0, name="c", lbound=-1, ubound=1 ),
       	mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.4, name="e", lbound=-1, ubound=1 ),

		# Matrices A, C, and E compute variance components
		mxAlgebra( expression=a %*% t(a), name="A" ),  # Variance Components
		mxAlgebra( expression=c %*% t(c), name="C" ),  # Variance Components
		mxAlgebra( expression=e %*% t(e), name="E" ),  # Variance Components
		mxAlgebra( expression=A+C+E, name="V" ),	     # Calculates the Phenotypic Variance
		mxAlgebra( expression=A/V, name="stndAVCs"), # Standardized A Variance Component
		#mxAlgebra( expression=C/V, name="stndCVCs"), # Standardized C Variance Component
		mxAlgebra( expression=E/V, name="stndEVCs"), # Standardized E Variance Component

     		# 95% Confidence Intervals
		mxCI(c("A","C","E","stndAVCs","stndEVCs")),	# Calculates the 95% confidence intervals for "stndVCs"

		# Means
		mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=0, label="mean", name="expMean" ),

    		# Algebra for expected variance/covariance matrix in MZ
   		mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),
							  cbind(A+C   , A+C+E)), name="expCovMZ" ),

    		# Algebra for expected variance/covariance matrix in DZ 
   		mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
 						        cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" )),

	mxModel("MZ",
	    		mxData( observed=MZdata, type="raw" ),
	            mxExpectationNormal( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selvars),  
	  		mxFitFunctionML()),

	mxModel("DZ",
	    		mxData( observed=DZdata, type="raw" ),
	            mxExpectationNormal( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selvars),  
	  		mxFitFunctionML()),

    		mxAlgebra( expression=MZ.objective + DZ.objective, name="twin" ),
	      mxFitFunctionAlgebra("twin")	# Calculates the -2LL
	)

	univAEFit <- mxRun(univAEModel, intervals=TRUE) 
	univAESumm <- summary(univAEFit)
	univAESumm

tableFitStatistics(univTwinSatFit , c(univACEFit,univAEFit))

tableFitStatistics(univACEFit,univAEFit)

#-------------------------------------------#
# Alternative way of testing a nested model #

# Fitting CE model
	univCEModel <- mxRename(univACEModel, "univCE")
	univCEModel$ACE.a <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,name="a")
	univCEFit <- mxRun(univCEModel,intervals=TRUE)
	univCESumm <- summary(univCEFit)
	univCESumm 

tableFitStatistics(univTwinSatFit , c(univACEFit,univAEFit,univCEFit))

tableFitStatistics(univACEFit,c(univAEFit,univCEFit))

#-------------------------------------------#
# Fit E model
	univEModel <- mxRename(univACEModel, "univE")
	univEModel$ACE.a <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,name="a")
	univEModel$ACE.c <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,name="c")
	univEFit <- mxRun(univEModel)
	univESumm <- summary(univEFit)
 	univESumm 

tableFitStatistics(univTwinSatFit , c(univACEFit,univAEFit,univCEFit,univEFit))

tableFitStatistics(univACEFit,c(univAEFit,univCEFit,univEFit))




##-----------------------------## Optional Output Tools ##-----------------------------##

## Summarizing Model Fits ##

# Saturated v ACE Model
tableFitStatistics(univTwinSatFit,univACEFit)

# ACE Model v Nested Models
univACENested <- list(univAEFit, univCEFit, univEFit)
tableFitStatistics(univACEFit,univACENested)

# Saturated v All Models
univACENested2 <- list(univACEFit, univAEFit, univCEFit, univEFit)
tableFitStatistics(univTwinSatFit,univACENested2)

## To see output for any fit, run the summaries we stored ##

univACESumm
univAESumm
univCESumm
univESumm

## See complete number of subjects ##

cpMZ # Complete MZ pairs
sMZ # MZ Singletons
cpDZ # Complete DZ pairs
sDZ # DZ Singletons

