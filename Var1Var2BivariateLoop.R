#-----------------------------------------------------------------------#
# Program: Bivariate AE Model for V1 and V2 MRI Measures                #
#  Author: Kelly Spoon                                                  #
#    Date: 04 02 2013                                                   #
#  Modified by Jeremy Elman 10 12 2015                                  #
#                                                                       #
# Description:							                                            #
#     	   							                                                #
#       Input: Raw data file (CSV) with the following columns	          #
#  		   IN ORDER: subid, case, twin, zygosity, phenotype 1,2,3...      #
#		   names of columns do not matter					                          #
#		   *twin variable should be defined A or B			                    #
#		   *zygosity should be defined as 1 for MZ and 2 for DZ	            #
#												                                                #
#	   Data: Continuous							    	                                #
#												                                                #
#      Output: All bivariate statistics                                 #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user       #
#-----------------------------------------------------------------------#

#*** Set working directory - location of data file and desired output
setwd("K:/Projects/Cortical_DTI/data/Merged/")

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

### Reading in Original Data ###

#*** Read in data file
twins <- read.csv("V2_MD_WMMD_AgeSiteAdj.csv",header=T)

#*** Quick check that data file is properly read in
	dim(twins)[1] 		# should be number of subjects
	dim(twins)[2] 		# should be 6 if formatted correctly
	names(twins)[1:4] 	# should be subid, case, twin, zygosity
	variables <- names(twins)[5:210]
	variables 

  V <- 103

#*** Running Transformations
  Time1vars <- variables[1:V]
  Time2vars <- variables[(V+1):(2*V)]
	variables2 <- cbind(Time1vars,Time2vars) 	
	variables2 # should be the phenotypes of interest, paired!

#*** Defining number of variables
nv <- 2				#Number of phenotypes (NVAR)
ntv <- nv*2

## Creating MZ and DZ data sets ##
# Dividing up twins
	twinA <- twins[twins$twin=="A",]
	twinB <- twins[twins$twin=="B",]

# Remerging by case to create paired data set
	newtwins <- merge(twinA, twinB, by=c("case","zyg14"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))

# Making data sets of Just MZ & DZ 
	mzdata <- as.data.frame(subset(newtwins,zyg14==1))
	dzdata <- as.data.frame(subset(newtwins,zyg14==2))


#initializing storage arrays
	ACESumm <- list()
  biACEout <- c()

#################################################
            ### START OF LOOP ###
#################################################

for(i in 1:V){
  
  # Print status
  cat("Starting: ", variables2[i,1]," ~ ", variables2[i,2], " (", i, "/", V, ")", sep="")
    
  selVars <- paste(variables2[i,],c("_A","_A","_B","_B"),sep="")
  
  # not necessary to make new data frames, but easier if you look at output
  MZdata <- mzdata[,selVars]
  DZdata <- dzdata[,selVars]
  
  # Determining means and standard deviations of the data for use as start values later
  totsd <- sd(c(as.matrix(MZdata),as.matrix(DZdata)), na.rm=TRUE)
  totmean <- mean(c(as.matrix(MZdata),as.matrix(DZdata)), na.rm=TRUE)
  
  #---------------------------------#
  ## Fit Bivariate Saturated Model ##
  #---------------------------------#
  
  meanSV <- c(0,0,0,0)  				# Start values for the means
  
  #-------------------------#
  # Fit Bivariate ACE Model #
  #-------------------------#
  
  multiCholACEModel <- mxModel("ACE",
      # Matrices a, c, and e to store a, c, and e path coefficients
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="a" ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, name="c" ),
          mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.3, name="e" ),
      # Matrices A, C, and E compute variance components
          mxAlgebra( expression=a %*% t(a), name="A" ),
          mxAlgebra( expression=c %*% t(c), name="C" ),
          mxAlgebra( expression=e %*% t(e), name="E" ),
      # Algebra to compute total variances and standard deviations (diagonal only)
          mxAlgebra( expression=A+C+E, name="V" ),
          mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
          mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
      # Confidence intervals portion for covariance matrices
  	    mxAlgebra( expression=(A/V)[1,1],name="A_V1"),
            mxAlgebra( expression=(A/V)[2,2],name="A_V2"),
       	    mxAlgebra( expression=(C/V)[1,1],name="C_V1"),
            mxAlgebra( expression=(C/V)[2,2],name="C_V2"),
            mxAlgebra( expression=(E/V)[1,1],name="E_V1"),
            mxAlgebra( expression=(E/V)[2,2],name="E_V2"),
      # Confidence intervals for Phenotypic correlation pieces
  	    mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V)))[1,2],name="Vcorr"),
  	    mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)))[1,2],name="Acorr"),
        mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C)))[1,2],name="Ccorr"),
  	    mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)))[1,2],name="Ecorr"),
  	    mxCI(c("A_V1","A_V2","C_V1","C_V2","E_V1","E_V2")),	    
  	    mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),
  	    mxCI(c("a","c","e")),
      ## Note that the rest of the mxModel statements do not change for bi/multivariate case
      # Matrix & Algebra for expected means vector
          mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, name="Mean" ),
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
  multiCholACEFit <- mxRun(multiCholACEModel,intervals=TRUE)
  multiCholACESumm <- summary(multiCholACEFit)
  #multiCholACESumm
  
  ACESumm[[i]]<- multiCholACESumm
  
  ACEfit = multiCholACEFit@output$Minus2LogLikelihood
  
  # Set r_A = 0 
  # -----------------------------------------------------------------------
  
  multiCholACEModel_noAcorr <- mxRename(multiCholACEModel, "ACE")
  multiCholACEModel_noAcorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), name="a")
  multiCholACEFit_noAcorr <- mxRun(multiCholACEModel_noAcorr)
  multiCholACESumm_noAcorr <- summary(multiCholACEFit_noAcorr)
  
  #Test of no Acorr
  noAcorrFit <- multiCholACEFit_noAcorr@output$Minus2LogLikelihood
  x2noAcorr<- noAcorrFit-ACEfit
  pnoAcorr<- pchisq(x2noAcorr,lower.tail=F,1)
  AICnoAcorr<- multiCholACESumm_noAcorr$AIC.Mx
  
  
  # Set r_E = 0 
  # -----------------------------------------------------------------------
  
  multiCholACEModel_noEcorr <- mxRename(multiCholACEModel, "ACE")
  multiCholACEModel_noEcorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), name="e")
  multiCholACEFit_noEcorr <- mxRun(multiCholACEModel_noEcorr)
  multiCholACESumm_noEcorr <- summary(multiCholACEFit_noEcorr)
  
  #Test of no Ecorr
  noEcorrFit <- multiCholACEFit_noEcorr@output$Minus2LogLikelihood
  x2noEcorr<- noEcorrFit-ACEfit
  pnoEcorr<- pchisq(x2noEcorr,lower.tail=F,1)
  AICnoEcorr<- multiCholACESumm_noEcorr$AIC.Mx
  
  
  # Test overall phenotypic correlation = 0 (i.e, r_A = r_C = r_E = 0)
  # -----------------------------------------------------------------------
  
  multiCholACEModel_nocorr <- mxRename(multiCholACEModel, "ACE")
  multiCholACEModel_nocorr$ACE.a <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), name="a")
  multiCholACEModel_nocorr$ACE.c <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(FALSE,FALSE,FALSE), values=c(0,0,0), name="c")
  multiCholACEModel_nocorr$ACE.e <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=c(TRUE,FALSE,TRUE), values=c(.3,0,.3), name="e")
  multiCholACEFit_nocorr <- mxRun(multiCholACEModel_nocorr)
  multiCholACESumm_nocorr <- summary(multiCholACEFit_nocorr)
  
  #No Phenotypic Corr
  ACEfit_nocorr<-multiCholACEFit_nocorr@output$Minus2LogLikelihood
  x2nocorr<- ACEfit_nocorr-ACEfit
  pnocorr<- pchisq(x2nocorr,lower.tail=F,2)
  AICace_nocorr<- multiCholACESumm_nocorr$AIC
  
  #-------------------------------------------------------------------#
  ### OUTPUT SUMMARY ###
  #-------------------------------------------------------------------#
  
  #ACE status
  IFAIL<- multiCholACEFit@output$status[[1]]
  
  
  ##Combine Output
  # b<-c(variables2[i,],IFAIL,multiCholACESumm$CI[1,1:3],multiCholACESumm$CI[2,1:3],multiCholACESumm$CI[3,1:3],
  #      multiCholACESumm$CI[4,1:3],multiCholACESumm$CI[5,1:3],multiCholACESumm$CI[6,1:3],multiCholACESumm$CI[7,],
  #      multiCholACESumm$CI[8,],multiCholACESumm$CI[9,],multiCholACESumm$CI[10,],multiCholACESumm$CI[11,],
  #      multiCholACESumm$CI[12,],multiCholACESumm$CI[13,])
  
  # Combine variable names, status, estimates and confidence intervals into single line
  b = multiCholACESumm$CI %>%  # Get CI data.frame
    dplyr::select(1:3) %>%  # Select estimate, lbound, ubound
    mutate(parameter = str_split_fixed(rownames(.),"\\[", n=2)[,1]) %>%  # Add parameter name to data.frame
    slice(1:10) %>%  # Drop path estimates (if retained, do not split parameter names above)
    gather(variable, value, -parameter) %>%  # Put into long format
    unite(temp, parameter, variable) %>% # Join parameter and variable name into one column
    arrange(temp) %>% # Order by parameter
    mutate(Variable1 = variables2[i,1], # Add phenotypic variable1 name 
           Variable2 = variables2[i,2], # Add phenotypic variable1 name 
           Status = IFAIL,              # Add status code
           pnoAcorr = pnoAcorr,         # Add p of no A correlation
           pnoEcorr = pnoEcorr,         # Add p of no E correlation
           pnocorr = pnocorr) %>%       # Add p of no phenotypic correlation
    spread(temp, value) # Transform into wide format (single line in this case)
  
  biACEout <- rbind(biACEout,b)
  

}
  
#################################################
           ### END OF LOOP ###
#################################################
  
write.csv(biACEout,"../../results/BivariateLoopACE_MD_WMMD_adj.csv",row.names=F)

#If you want to look at any particular fit you can use, just replace 1 with the corresponding row in bivariate

# ACESumm[[14]]
# ACESumm[[30]]
# ACESumm[[1]]
# ACESumm_Acorr[[1]]
