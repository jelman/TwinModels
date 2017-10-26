#-----------------------------------------------------------------------#
# Program: Univariate Saturated and ACE Model with Loop	                #
#  Author: Kelly Spoon                                                  #
#    Date: 07 15 2011                                                   #
#                                                                       #
# Modified by Jeremy Elman 10 07 2015                                 #
#                                                                       #
# Description: Complete univariate analysis (Saturated, ACE model,    	#
#		   and Nested models) for multiple phenotypes		                    #							                 
#												                                                #
#       Input: Raw data file (CSV) with the following columns	 	        #
#  		   IN ORDER: subid, case, twin, zygo09, phenotype 1,2,3....       #
#		   names of columns do not matter					                          #
#		   *twin variable should be defined A or B			                    #
#		   *zygosity should be defined as 1 for MZ and 2 for DZ	            #
#												                                                #
#	   Data: Continuous								                                    #
#												                                                #
#      Output: All univariate statistics are exported to a .csv file	  #
#		   "LOOPoutput.csv" in the working directory                        #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user       #
#												                                                #
#	WARNING: This script allows for a univariate analysis to be run	      #
#		   on multiple phenotyes at once, but does guarantee that 	        #
#		   phenotypes are clean and that the results will be free           #
#		   of errors.  The user needs to properly examine the data 	        #
#		   prior to the twin analyses, and closely inspect the	            #    
#		   resulting output for potential errors.				                    #
#												                                                #
#	   Note: Start values are estimated from the data as part of 	        #
#		   the script.								                                      #
#-----------------------------------------------------------------------#

#*** Set working directory - location of data file and desired output
	setwd("K:/Projects/Cortical_DTI/data")

## Loading Required/Useful Libraries ##

# Load OpenMX
	require(OpenMx)   	#Loads OpenMx
	require(psych)   		#Loads Psych package
	mxOption(NULL, "Default optimizer", "NPSOL")
	
# Load additional GenEpi tools
	source("K:/code/GenEpiHelperFunctions.R")

### Reading in Original Data ###

#*** Read in data file
	twins <- read.csv("V2_MD_CT_AgeSiteAdj.csv",header=T)

#*** Quick check that data file is properly read in
	dim(twins)[1] 		#Should be number of subjects 
	names(twins)[1:4] 	#Should be subid, case, twin, zygosity

#*** Defining number of variables
	tot <- dim(twins)[2]-4 
	tot 				#Should be number of variables to loop through
	names(twins)[5:(4+tot)] #Variables to be used

#*** Summary statistics for the data file
	describe(twins)

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg14"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
	names(newtwins)		#Variable names in the now paired dataset

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg14==1))
	DZdata <- as.data.frame(subset(newtwins,zyg14==2))

# Renaming Data Before Starting Loop
	MZdata1 <- MZdata
	DZdata1 <- DZdata


#-----------------------------------------------------------------------#

### START OF LOOP ###
#*** YOU MUST RUN FROM HERE TO END OF LOOP ***#

#-----------------------------------------------------------------------#


# Initializing Output 
	output <- c()

for (i in 1:tot){
  
# Variable name
varname <- names(twins)[4+i]
  
# Print status
cat("Starting: ", varname, " (", i, "/", tot, ")", sep="")

  ## Defining Variables / Data for use with OpenMX ##

# Making Sparse Data Frame with only phenotype 
	MZdata <- MZdata1[,c(4+i,tot+6+i)]
	DZdata <- DZdata1[,c(4+i,tot+6+i)]

# Getting Number of Twins Used in Analysis
	cpMZ <- sum(complete.cases(MZdata)) # Number of complete MZ pairs
	sMZ <- dim(MZdata)[1]-cpMZ 		# Number of MZ singletons
	cpDZ <- sum(complete.cases(DZdata)) # Number of complete DZ pairs
	sDZ <- dim(DZdata)[1]-cpDZ 		# Number of DZ singletons

# Define Variable to use in OpenMX
	selvars <- c(names(MZdata))
	Vars <- 'var'
	nv <- 1 					# Number of phenotypes
	ntv <- nv*2
	
# Correlation between MZ twins
	MZcor<-cor(MZdata,use="complete")

# Correlation between DZ twins 
	DZcor<-cor(DZdata,use="complete")

# Determining means and standard deviations of the data for use as start values later
	totsd <- sd(na.omit(c(MZdata[,1],MZdata[,2],DZdata[,1],DZdata[,2])))
	totmean <- mean(na.omit(c(MZdata[,1],MZdata[,2],DZdata[,1],DZdata[,2])))

#-----------------------------------------------------------------------#

### Univariate Analysis ###

## Fit Univariate Saturated Model ##

univTwinSatModel <- mxModel("univTwinSat",
    mxModel("MZ",
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=totsd, name="CholMZ" ),
        mxAlgebra( expression=CholMZ %*% t(CholMZ), name="expCovMZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=totmean, name="expMeanMZ" ),
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
        mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=totsd, name="CholDZ" ),
        mxAlgebra( expression=CholDZ %*% t(CholDZ), name="expCovDZ" ),
        mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=totmean, name="expMeanDZ" ),
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

## Univariate ACE Model ##

twinACEModel <- mxModel("ACE",
		# Matrices X, Y, and Z to store a, c, and e path coefficients
		mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=totsd/3, label="a", name="X" ),
		mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=totsd/3, label="c", name="Y" ),
		mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=totsd/3, label="e", name="Z" ),
		# Matrices A, C, and E compute variance components
		mxAlgebra( expression=X %*% t(X), name="A" ),
		mxAlgebra( expression=Y %*% t(Y), name="C" ),
		mxAlgebra( expression=Z %*% t(Z), name="E" ),
		mxAlgebra( expression=A+C+E, name="V" ),
		mxAlgebra( expression=cbind(A/V,C/V,E/V),name="stndVCs"),
     		# Calculate 95% CIs here
		mxCI(c("stndVCs")),
		mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= totmean, label="mean", name="expMean" ),
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
    		mxAlgebra( expression=MZ.objective + DZ.objective, name="twin" ),
		mxAlgebra( expression=MZ.objective + DZ.objective, name="neg2sumll" ),
		mxFitFunctionAlgebra("neg2sumll")
		)

	twinACEFit <- mxRun(twinACEModel,intervals=TRUE) # Additional indicator (interval=T) to calculate 95% CI must be included in mxRun statement
	univACESumm <- summary(twinACEFit)

# Fitting AE model
	twinAEModel <- mxRename(twinACEModel, "twinAE")
	twinAEModel$twinAE.Y <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,label="c")
	twinAEFit <- mxRun(twinAEModel,intervals=TRUE)
	univAESumm <- summary(twinAEFit)

# Fitting CE model
	twinCEModel <- mxRename(twinACEModel, "twinCE")
	twinCEModel$twinCE.X <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,label="a")
	twinCEFit <- mxRun(twinCEModel,intervals=TRUE)
	univCESumm <- summary(twinCEFit)

# Fit E model
	twinEModel <- mxRename(twinACEModel, "twinE")
	twinEModel$twinE.X <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,label="a")
	twinEModel$twinE.Y <- mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=0,label="c")
	twinEFit <- mxRun(twinEModel)
	univESumm <- summary(twinEFit)

#-----------------------------------------------------------------------#
## Generating the values that will go in the output file ##



#Sat Model
sat2ll <- univTwinSatSumm$Minus2LogLikelihood
satdf <- univTwinSatSumm$degreesOfFreedom

#ACE
IFAIL<- twinACEFit@output$status[[1]]
ACEfit<- twinACEFit@output$Minus2LogLikelihood
AICace<- univACESumm$AIC
x2acefit<- ACEfit-sat2ll
pacefit<- pchisq(ACEfit-sat2ll,lower.tail=F,6)

#CE (test of no A)
CEfit<- twinCEFit@output$Minus2LogLikelihood
x2noA<- CEfit-ACEfit
pnoA<- pchisq(CEfit-ACEfit,lower.tail=F,1)
AICce<- univCESumm$AIC

#AE (test of no C)
AEfit<-twinAEFit@output$Minus2LogLikelihood
x2noC<- AEfit-ACEfit
pnoC<- pchisq(AEfit-ACEfit,lower.tail=F,1)
AICae<- univAESumm$AIC

#E Only (test of no A and C)
Efit<-twinEFit@output$Minus2LogLikelihood
x2noAC<- Efit-ACEfit
pnoAC<- pchisq(Efit-ACEfit,lower.tail=F,2)
AICe<- univESumm$AIC

#Expected Twin 1 Means
ExpMean1<-as.numeric(twinACEFit$ACE.expMean@values[1,1])

## Unstandardized ACE Parameter Estimates
#*** To change these, just change twinACEFit to best fit (i.e. twinCEFit)
a1<-twinACEFit@output$algebras$ACE.A
c1<-twinACEFit@output$algebras$ACE.C
e1<-twinACEFit@output$algebras$ACE.E
v1<-twinACEFit@output$algebras$ACE.V

## Standardized ACE Parameter Estimates
a2CI<-cbind(twinACEFit@output$algebras$ACE.stndVCs[1,1],twinACEFit@output$confidenceIntervals[1,1],twinACEFit@output$confidenceIntervals[1,3])
c2CI<-cbind(twinACEFit@output$algebras$ACE.stndVCs[1,2],twinACEFit@output$confidenceIntervals[2,1],twinACEFit@output$confidenceIntervals[2,3])
e2CI<-cbind(twinACEFit@output$algebras$ACE.stndVCs[1,3],twinACEFit@output$confidenceIntervals[3,1],twinACEFit@output$confidenceIntervals[3,3])

##Combine Output
full<-c(varname,MZcor[2,1],DZcor[2,1],IFAIL,ACEfit,x2acefit,pacefit,AICace,CEfit,x2noA,pnoA,AICce,
AEfit,x2noC,pnoC,AICae,Efit,x2noAC,pnoAC,AICe,ExpMean1,a1,c1,e1,v1,a2CI,c2CI,e2CI)

output <- rbind(output,full)

}

### END OF LOOP ###

output <- as.data.frame(output,row.names=c(1:tot))
names(output)<-c("pheno","MZcor","DZcor","IFAIL","ACEfit","x2acefit","pacefit","AICace","CEfit","x2noA","pnoA","AICce",
"AEfit","x2noC","pnoC","AICae","Efit","x2noAC","pnoAC","AICe","ExpMean1","A","C","E","V",
"a2","a2Lower","a2Upper","c2","c2Lower","c2Upper","e2","e2Lower","e2Upper")

#Make an output file for use, will be created in directory specified at beginning
write.csv(output,"../results/LOOPoutput.csv",row.names=F)
