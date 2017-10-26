#-----------------------------------------------------------------------#
# Program: Trivariate ACE Model for V1, V2 and V3 MRI Measures          #
#  Author: Jeremy Elman                                                 #
#    Date: 11 04 2015                                                   #
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
twins <- read.csv("./0.8_2mm_Weighted/V2_CT_WMMD_CMD_AgeSiteAdj.csv",header=T)

#*** Quick check that data file is properly read in
dim(twins)[1] 		# should be number of subjects
dim(twins)[2] 		# should be 6 if formatted correctly
names(twins)[1:4] 	# should be subid, case, twin, zygosity
variables <- names(twins)[5:ncol(twins)] # Can set the columns explicitly if files contains extras
variables 


#*** Defining number of variables
nv <- 3				#Number of phenotypes (NVAR)
ntv <- nv*2   #Number of total phenotypes (phenotypes x 2 twins)
V <- length(variables) / nv #Number of measurements of each phenotype. Can set explicitly      

#*** Running Transformations
Time1vars <- variables[1:V]
Time2vars <- variables[(V+1):(2*V)]
Time3vars <- variables[((2*V)+1):(3*V)]
variables2 <- cbind(Time1vars,Time2vars,Time3vars) 	
variables2 # should be the phenotypes of interest, paired!



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
TriACEout <- c()

#################################################
### START OF LOOP ###
#################################################

for(i in 1:V){
  
  # Print status
  cat("Starting: ", variables2[i,1]," ~ ", variables2[i,2], " ~ ", variables2[i,3], " (", i, "/", V, ")", sep="")
  
  selVars <- paste(variables2[i,],c(rep("_A",nv),rep("_B",nv)),sep="")
  
  # not necessary to make new data frames, but easier if you look at output
  MZdata <- mzdata[,selVars]
  DZdata <- dzdata[,selVars]
  
  # Determining means and standard deviations of the data for use as start values later
  totsd <- sd(c(as.matrix(MZdata),as.matrix(DZdata)), na.rm=TRUE)
  totmean <- mean(c(as.matrix(MZdata),as.matrix(DZdata)), na.rm=TRUE)
  
  # Set start values
  meanSVnv <- rep(0, ntv)  				# Start values for the means
  cholASVnv <- diag(.3, nrow=nv, ncol=nv)  #Set value for A variances
  cholASVnv[lower.tri(cholASVnv)] <- .1    #Set value for A covariances
  cholCSVnv <- diag(.3, nrow=nv, ncol=nv)  #Set value for C variances
  cholCSVnv[lower.tri(cholCSVnv)] <- .1    #Set value for C covariances
  cholESVnv <- diag(.3, nrow=nv, ncol=nv)  #Set value for E variances
  cholESVnv[lower.tri(cholESVnv)] <- .1    #Set value for E covariances
  
  # Set labels
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
  
  #--------------------------#
  # Fit Trivariate ACE Model #
  #--------------------------#
  
  multiCholACEModel <- mxModel("ACE",
     # Matrices a, c, and e to store a, c, and e path coefficients
     mxMatrix(type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholASVnv, labels=AFac, name="a", lbound=-3, ubound=3),
     mxMatrix(type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholCSVnv, labels=CFac, name="c", lbound=-3, ubound=3),
     mxMatrix(type="Lower", nrow=nv, ncol=nv, free=TRUE, values=cholESVnv, labels=EFac, name="e", lbound=-3, ubound=3),
     
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
             mxExpectationNormal( covariance="ACE.expCovMZ", means="ACE.expMean", dimnames=selVars),  
             mxFitFunctionML()
     ),
     
     mxModel("DZ", 
             mxData( observed=DZdata, type="raw" ),
             mxExpectationNormal( covariance="ACE.expCovDZ", means="ACE.expMean", dimnames=selVars),  
             mxFitFunctionML()
     ),
     mxAlgebra( expression=MZ.objective + DZ.objective, name="neg2sumll" ),
     mxFitFunctionAlgebra("neg2sumll")
  )
  
  multiCholACEFit <- mxRun(multiCholACEModel,intervals=TRUE)
  multiCholACESumm <- summary(multiCholACEFit)
  #multiCholACESumm
  
  ACESumm[[i]]<- multiCholACESumm
  IFAIL<- multiCholACEFit@output$status[[1]]
  
  
  #-------------------------------------------------------------------#
  ### OUTPUT SUMMARY ###
  #-------------------------------------------------------------------#
  
  ##Combine Output
  # b<-c(variables2[i,],IFAIL,multiCholACESumm$CI[1,1:3],multiCholACESumm$CI[2,1:3],multiCholACESumm$CI[3,1:3],
  #      multiCholACESumm$CI[4,1:3],multiCholACESumm$CI[5,1:3],multiCholACESumm$CI[6,1:3],multiCholACESumm$CI[7,],
  #      multiCholACESumm$CI[8,],multiCholACESumm$CI[9,],multiCholACESumm$CI[10,],multiCholACESumm$CI[11,],
  #      multiCholACESumm$CI[12,],multiCholACESumm$CI[13,])
  
  # Combine variable names, status, estimates and confidence intervals into single line
  b = multiCholACESumm$CI %>%  # Get CI data.frame
    dplyr::select(1:3) %>%  # Select estimate, lbound, ubound
    mutate(parameter = str_split_fixed(rownames(.),"\\[", n=2)[,1]) %>%  # Add parameter name to data.frame
    #slice(1:10) %>%  # Drop path estimates (if retained, do not split parameter names above)
    gather(variable, value, -parameter) %>%  # Put into long format
    unite(temp, parameter, variable) %>% # Join parameter and variable name into one column
    arrange(temp) %>% # Order by parameter
    mutate(Variable1 = variables2[i,1], # Add phenotypic variable1 name 
           Variable2 = variables2[i,2], # Add phenotypic variable2 name 
           Variable3 = variables2[i,3], # Add phenotypic variable3 name 
           Status = IFAIL) %>%          # Add status code
    spread(temp, value) # Transform into wide format (single line in this case)
  
  TriACEout <- rbind(TriACEout,b)
  
  
}

#################################################
### END OF LOOP ###
#################################################

write.csv(TriACEout,"../../results/TrivariateLoopACE_CT_WMMD_CMD_adj.csv",row.names=F)