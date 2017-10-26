#-----------------------------------------------------------------------#
# Program: Bivariate Saturated and ACE Models with LOOP 			    #
#	     (One phenotype is fixed, all others cycle through		        # 
#  Author: Kelly Spoon                                                  #
#    Date: 08 01 2011                                                   #
#                                                                       #
# Description: Complete bivariate analysis of Saturated and ACE models, #
#		   including all nested models contraining correlation	        #
#     	   							                                    #
#     	   Loops through all bivariate analyses for one fixed           #
#     	   phenotype (primary phenotypes) with all remaining 	        #
#     	   phenotype variables 				                            #
#     	   							                                    #
#       Input: Raw data file (CSV) with the following columns	 	    #
#  		   IN ORDER: subid, case, twin, zygosity, phenotypes		    #
#		   names of columns do not matter					            #
#		   *primary phenotype should come first, in column 5 		    #
#		   *twin variable should be defined A or B			            #
#		   *zygosity should be defined as 1 for MZ and 2 for DZ	        #
#												                        #
#        Data: Continuous								                #
#												                        #
#      Output: All bivariate statistics                                 #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user     #
#												                        #
#	WARNING: The user needs to closely examine the output in order to   #
#		   verify that all models have run successfully			        #
#-----------------------------------------------------------------------#

#*** Set working directory - location of data file and desired output
	setwd("C:/Documents and Settings/mspanizzon/Desktop/OpenMx/Library")

## Loading Required/Useful Libraries ##
require(OpenMx)		# Loads OpenMX
require(psych)		# Loads Psych package
source("http://www.vipbg.vcu.edu/~vipbg/Tc24/GenEpiHelperFunctions.R") #Loads GenEpi package from VCU

### Reading in Original Data ###
#*** Read in data file
	twins <- read.csv("bivACE_Loop.csv",header=T)

#*** Quick check that data file is properly read in
	dim(twins)[1] 		# should be number of subjects
	names(twins)[1:4] 	# should be subid, case, twin, zygosity
	names(twins)[5] 		# should be primary phenotype, will be included in ALL analyses
	names(twins[6:dim(twins)[2]]) # should be the remaining phenotypes

#*** Defining number of variables
	tot <- dim(twins)[2]-5 
	tot 			# should be number of variables to loop through

#*** Defining number of variables
nv <- 2		#Number of phenotypes (NVAR)
ntv <- nv*2

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg09"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))

# Making data sets of Just MZ & DZ 
	MZdata <- as.data.frame(subset(newtwins,zyg09==1))
	DZdata <- as.data.frame(subset(newtwins,zyg09==2))

# Renaming Data Before Starting Loop
	MZdata1 <- MZdata
	DZdata1 <- DZdata


#-----------------------------------------------------------------------#

### START OF LOOP ###
#*** YOU MUST RUN FROM HERE TO END OF LOOP ***#

#-----------------------------------------------------------------------#


# Initializing Output 
	output <- c()

for (k in 1:tot){

    ## Defining Variables / Data for use with OpenMX ##
    
    # Making Sparse Data Frame with only phenotype 
    	MZdata <- MZdata1[,c(5,5+k,tot+8,tot+8+k)]
    	DZdata <- DZdata1[,c(5,5+k,tot+8,tot+8+k)]
    # Getting Number of Twins Used in Analysis
    	cpMZ <- sum(complete.cases(MZdata)) # Number of complete MZ pairs
    	sMZ <- dim(MZdata)[1]-cpMZ # Number of MZ singletons
    	cpDZ <- sum(complete.cases(DZdata)) # Number of complete DZ pairs
    	sDZ <- dim(DZdata)[1]-cpDZ # Number of DZ singletons
    # Define Variable to use in OpenMX
    	selVars <- c(names(MZdata))
    
    # Getting Number of Twins Used in Analysis
    	cpMZ <- sum(complete.cases(MZdata)) # Number of complete MZ pairs
    	sMZ <- dim(MZdata)[1]-cpMZ # Number of MZ singletons
    	cpDZ <- sum(complete.cases(DZdata)) # Number of complete DZ pairs
    	sDZ <- dim(DZdata)[1]-cpDZ # Number of DZ singletons
    
    	totsd <- sd(na.omit(c(MZdata[,1],MZdata[,2],DZdata[,1],DZdata[,2])))
    	totmean <- mean(na.omit(c(MZdata[,1],MZdata[,2],DZdata[,1],DZdata[,2])))
    
    pheno1 <- as.character((names(twins[5])))
    pheno2 <- as.character((names(twins[5+k])))
    print(pheno2)
    
    # Fit Bivariate Saturated Model
    # -----------------------------------------------------------------------
    meanSV <- c(5,5)			# Start values for the means
    cholSV <- c(.9,0,0,0,.9,0,0,.9,0,.9) # Start values for the variance/covariance matrix
    
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
            mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=cholSV, labels=labM, name="CholMZ" ),
            mxAlgebra( expression=CholMZ %*% t(CholMZ), name="expCovMZ" ),
            mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanSV, labels=meanM, name="expMeanMZ" ),
            mxData( observed=MZdata, type="raw" ),
            mxExpectationNormal( covariance="expCovMZ", means="expMeanMZ", dimnames=selVars),
            mxFitFunctionML(),
        # Algebra's needed for equality constraints    
            mxAlgebra( expression=t(diag2vec(expCovMZ)), name="expVarMZ"),
            mxAlgebra( expression=expVarMZ[1,1:nv], name="expVarMZtA"),
            mxAlgebra( expression=expVarMZ[1,(nv+1):ntv], name="expVarMZtB"),
            mxAlgebra( expression=expMeanMZ[1,1:nv], name="expMeanMZtA"),
            mxAlgebra( expression=expMeanMZ[1,(nv+1):ntv], name="expMeanMZtB")
        ),
        mxModel("DZ",
            mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=TRUE, values=cholSV, labels=labD, name="CholDZ" ),
            mxAlgebra( expression=CholDZ %*% t(CholDZ), name="expCovDZ" ),
            mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=meanSV, labels=meanD, name="expMeanDZ" ),
            mxData( observed=DZdata, type="raw" ),
            mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ", dimnames=selVars),
            mxFitFunctionML(),
        # Algebra's needed for equality constraints    
            mxAlgebra( expression=t(diag2vec(expCovDZ)), name="expVarDZ"),
            mxAlgebra( expression=expVarDZ[1,1:nv], name="expVarDZtA"),
            mxAlgebra( expression=expVarDZ[1,(nv+1):ntv], name="expVarDZtB"),
            mxAlgebra( expression=expMeanDZ[1,1:nv], name="expMeanDZtA"),
            mxAlgebra( expression=expMeanDZ[1,(nv+1):ntv], name="expMeanDZtB")
        ),
        mxAlgebra( MZ.objective + DZ.objective, name="neg2sumll" ),
        mxFitFunctionAlgebra("neg2sumll")
    )
    
    multiTwinSatFit <- mxRun(multiTwinSatModel)
    multiTwinSatSumm <- summary(multiTwinSatFit)
    #multiTwinSatSumm
    
    # Generate Saturated Output
    #parameterSpecifications(multiTwinSatFit)
    #expectedMeansCovariances(multiTwinSatFit)
    #tableFitStatistics(multiTwinSatFit)
    
    # Fit Bivariate ACE Model with Raw Data Input
    # -----------------------------------------------------------------------
    meanSVnv <- rep(5,nv)		# Start values for the means
    cholASVnv <- rep(.3,nv*(nv+1)/2)	# Start values for the genetic and environmental parameter estimates
    cholCSVnv <- rep(.1,nv*(nv+1)/2)
    cholESVnv <- rep(.3,nv*(nv+1)/2)
    
    rows <- c()
    cols <- c()
    for (i in 1:nv){
    row <- c(i:nv)
    col <- rep(i,(nv+1-i))
    rows <- c(rows,row)
    cols <- c(cols,col)
    }
    
    AFac <- paste("A",rows,cols,sep="")
    CFac <- paste("C",rows,cols,sep="")
    EFac <- paste("E",rows,cols,sep="")
    mean <- as.character(c("mean1","mean2"))
    
    multiCholACEModel <- mxModel("ACE",
        # Matrices a, c, and e to store a, c, and e path coefficients
            mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, labels=AFac, name="a" ),
            mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, labels=CFac, name="c" ),
            mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, labels=EFac, name="e" ),
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
            mxExpectationNormal( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selVars ),
            mxFitFunctionML()
        ),
        mxModel("DZ", 
            mxData( observed=DZdata, type="raw" ),
            mxExpectationNormal( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selVars ),
            mxFitFunctionML()
        ),
        mxAlgebra( expression=MZ.objective + DZ.objective, name="neg2sumll" ),
        mxFitFunctionAlgebra("neg2sumll")
    )
    multiCholACEFit <- mxRun(multiCholACEModel,intervals=True)	# Allows for confidence intervals to be calculated
    multiCholACESumm <- summary(multiCholACEFit)
    #multiCholACESumm
    
    # Generate Bivariate Cholesky ACE Output
    #parameterSpecifications(multiCholACEFit)
    #expectedMeansCovariances(multiCholACEFit)
    #tableFitStatistics(multiTwinSatFit, multiCholACEFit)
    
    # Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices
    
    # Set r_C = 0 
    # -----------------------------------------------------------------------
    
    multiCholACEModel_noCcorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_noCcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c" )
    multiCholACEFit_noCcorr <- mxRun(multiCholACEModel_noCcorr)
    multiCholACESumm_noCcorr <- summary(multiCholACEFit_noCcorr)
    #multiCholACESumm_noCcorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_noCcorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_noCcorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_noCcorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    # Set r_A = 0 
    # -----------------------------------------------------------------------
    
    multiCholACEModel_noAcorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_noAcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
    multiCholACEFit_noAcorr <- mxRun(multiCholACEModel_noAcorr)
    multiCholACESumm_noAcorr <- summary(multiCholACEFit_noAcorr)
    #multiCholACESumm_noAcorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_noAcorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_noAcorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_noAcorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    # Set r_E = 0 
    # -----------------------------------------------------------------------
    
    multiCholACEModel_noEcorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_noEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
    multiCholACEFit_noEcorr <- mxRun(multiCholACEModel_noEcorr)
    multiCholACESumm_noEcorr <- summary(multiCholACEFit_noEcorr)
    #multiCholACESumm_noEcorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_noEcorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_noEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_noEcorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    # Set r_A = 0 & r_C = 0
    # -----------------------------------------------------------------------
    
    multiCholACEModel_noACcorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_noACcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
    multiCholACEModel_noACcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c")
    multiCholACEFit_noACcorr <- mxRun(multiCholACEModel_noACcorr)
    multiCholACESumm_noACcorr <- summary(multiCholACEFit_noACcorr)
    #multiCholACESumm_noACcorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_noACcorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_noACcorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_noACcorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    # Set r_A = 0 & r_E = 0
    # -----------------------------------------------------------------------
    
    multiCholACEModel_noAEcorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_noAEcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
    multiCholACEModel_noAEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
    multiCholACEFit_noAEcorr <- mxRun(multiCholACEModel_noACcorr)
    multiCholACESumm_noAEcorr <- summary(multiCholACEFit_noACcorr)
    #multiCholACESumm_noAEcorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_noAEcorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_noAEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_noACcorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    # Set r_E = 0 & r_C = 0
    # -----------------------------------------------------------------------
    
    multiCholACEModel_noCEcorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_noCEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
    multiCholACEModel_noCEcorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c")
    multiCholACEFit_noCEcorr <- mxRun(multiCholACEModel_noCEcorr)
    multiCholACESumm_noCEcorr <- summary(multiCholACEFit_noCEcorr)
    #multiCholACESumm_noCEcorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_noCEcorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_noCEcorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_noCEcorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    # Test overall phenotypic correlation = 0 (i.e, r_A = r_C = r_E = 0)
    # -----------------------------------------------------------------------
    
    multiCholACEModel_nocorr <- mxRename(multiCholACEModel, "ACE")
    multiCholACEModel_nocorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=AFac, name="a")
    multiCholACEModel_nocorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.1,0,.1), labels=CFac, name="c")
    multiCholACEModel_nocorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), labels=EFac, name="e")
    multiCholACEFit_nocorr <- mxRun(multiCholACEModel_nocorr)
    multiCholACESumm_nocorr <- summary(multiCholACEFit_nocorr)
    #multiCholACESumm_nocorr
    #tableFitStatistics(multiCholACEFit,multiCholACEFit_nocorr)
    
    # Cholesky Output (Path Estimates)
    #formatOutputMatrices(multiCholACEFit_nocorr,ACEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    #formatOutputMatrices(multiCholACEFit_nocorr,ACEcovMatrices,ACEcovLabels,Vars,4)
    
    
    # Set C = 0, Bivariate AE Model
    # -----------------------------------------------------------------------
    
    multiCholAEModel <- mxRename(multiCholACEModel, "AE")
    multiCholAEModel$AE.c <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=FALSE, values=0, name="c" )
    multiCholAEModel$AE.Ccorr <- mxAlgebra(0,name="Ccorr")
    multiCholAEFit <- mxRun(multiCholAEModel,intervals=T)
    multiCholAESumm <- summary(multiCholAEFit)
    #multiCholAESumm
    #tableFitStatistics(multiCholACEFit,multiCholAEFit)
    
    # Cholesky Output (Path Estimates)
    AEpathMatrices <- c("AE.a","AE.c","AE.e","AE.iSD","AE.iSD %*% AE.a","AE.iSD %*% AE.c","AE.iSD %*% AE.e")
    #formatOutputMatrices(multiCholAEFit,AEpathMatrices,ACEpathLabels,Vars,4)
    
    # Correlated Factors Output (Cov and Corr Matrices)
    AEcovMatrices <- c("AE.A","AE.C","AE.E","AE.V","AE.A/AE.V","AE.C/AE.V","AE.E/AE.V")
    #formatOutputMatrices(multiCholAEFit,AEcovMatrices,ACEcovLabels,Vars,4)
    
    #-------------------------------------------------------------------#
    ### OUTPUT ###
    #-------------------------------------------------------------------#
    
    sat2ll <- multiTwinSatSumm$Minus2LogLikelihood
    satdf <- multiTwinSatFit$degreesOfFreedom
    
    #ACE
    IFAIL<- multiCholACEFit@output$status[[1]]
    ACEfit<- multiCholACEFit@output$Minus2LogLikelihood
    AICace<- multiCholACESumm$AIC
    x2acefit<- ACEfit-sat2ll
    pacefit<- pchisq(ACEfit-sat2ll,lower.tail=F,17)
    
    
    # Path Estimates, Standardized
    sd <- multiCholACEFit@output$algebra$ACE.iSD 
    a <- multiCholACEFit@output$matrices$ACE.a
    c <- multiCholACEFit@output$matrices$ACE.c
    e <- multiCholACEFit@output$matrices$ACE.e
    Apath <- sd %*% a
    Cpath <- sd %*% c
    Epath <- sd %*% e
    LP1_pheno1_A <- Apath[1,1]
    LP1_pheno2_A <- Apath[2,1]
    LP2_pheno2_A <- Apath[2,2]
    LP1_pheno1_C <- Cpath[1,1]
    LP1_pheno2_C <- Cpath[2,1]
    LP2_pheno2_C <- Cpath[2,2]
    LP1_pheno1_E <- Epath[1,1]
    LP1_pheno2_E <- Epath[2,1]
    LP2_pheno2_E <- Epath[2,2]
    
    # Correlated Factors Output (Cov and Corr Matrices), Standardized
    
    A <- multiCholACEFit@output$algebra$ACE.A 
    C <- multiCholACEFit@output$algebra$ACE.C 
    E <- multiCholACEFit@output$algebra$ACE.E
    V <- multiCholACEFit@output$algebra$ACE.V  
    
    Acov <- A/V
    Ccov <- C/V
    Ecov <- E/V
    A1 <- Acov[1,1]
    A2 <- Acov[2,2]
    A12 <- Acov[1,2]
    C1 <- Ccov[1,1]
    C2 <- Ccov[2,2]
    C12 <- Ccov[1,2]
    E1 <- Ecov[1,1]
    E2 <- Ecov[2,2]
    E12 <- Ecov[1,2]
    
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
    bivariate<-c(pheno1,pheno2,sat2ll,ACEfit,x2acefit,pacefit,AICace,
    LP1_pheno1_A, LP1_pheno2_A, LP2_pheno2_A, LP1_pheno1_C, LP1_pheno2_C, LP2_pheno2_C, LP1_pheno1_E, LP1_pheno2_E, LP2_pheno2_E,
    A1, A12, A2, C1, C12, C2, E1, E12, E2,
    ACEfit_noAcorr,x2noAcorr,pnoAcorr,AICace_noAcorr,
    ACEfit_noCcorr,x2noCcorr,pnoCcorr,AICace_noCcorr,ACEfit_noEcorr,x2noEcorr,pnoEcorr,AICace_noEcorr,
    ACEfit_noACcorr,x2noACcorr,pnoACcorr,AICace_noACcorr,
    ACEfit_noAEcorr,x2noAEcorr,pnoAEcorr,AICace_noAEcorr,
    ACEfit_noCEcorr,x2noCEcorr,pnoCEcorr,AICace_noCEcorr,
    ACEfit_nocorr,x2nocorr,pnocorr,AICace_nocorr)
    
    output <- rbind(output, bivariate)
}

## END OF LOOP ##

output <- as.data.frame(output, row.names = c(1:tot))

names(output)<-c("variable1","variable2","satfit","ACEfit","x2acefit","pacefit","AICace",
"stpathLP1_pheno1_A","stpathLP1_pheno2_A","stpathLP2_pheno2_A",
"stpathLP1_pheno1_C","stpathLP1_pheno2_C","stpathLP2_pheno2_C",
"stpathLP1_pheno1_E","stpathLP1_pheno2_E","stpathLP2_pheno2_E",
"stCovA11", "stCovA12", "stCovA22", "stCovC11", "stCovC12", "stCovC22", "stCovE11", "stCovE12", "stCovE22", 
"ACEfit_noAcorr","x2noAcorr","pnoAcorr","AICace_noAcorr",
"ACEfit_noCcorr","x2noCcorr","pnoCcorr","AICace_noCcorr",
"ACEfit_noEcorr","x2noEcorr","pnoEcorr","AICeace_noEcorr",
"ACEfit_noACcorr","x2noACcorr","pnoACcorr","AICace_noACcorr",
"ACEfit_noAEcorr","x2noAEcorr","pnoAEcorr","AICace_noAEcorr",
"ACEfit_noCEcorr","x2noCEcorr","pnoCEcorr","AICeace_noCEcorr",
"ACEfit_nocorr","x2nocorr","pnocorr","AICace_nocorr")

write.csv(output, "bivariateLOOPout.csv", row.names=F)

