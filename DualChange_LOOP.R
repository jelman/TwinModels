#-----------------------------------------------------------------------#
# Program: Dual change and ACE Models LOOP, ALL combinations	          #
#		Continuous Data								                                      #	
#  Author: Kelly Spoon                                                  #
#  Edited by: Jeremy Elman                                              #
#    Date: 07 28 2011                                                   #
#                                                                       #
# Description: Complete dual change analyses  ACE models,               #
#		   including all nested models contraining correlation            	#
#     	   							                                                #
#     	   Loops through all possible pairwise combinations of 
#          user-provided phenotype variables. Estimates level and
#          change #
#     	   							                  #
#       Input: Raw data file (CSV) with the following columns	 	#
#  		   IN ORDER: subid, case, twin, zygosity, phenotypes		#
#		   names of columns do not matter					#
#		   *twin variable should be defined A or B			#
#		   *zygosity should be defined as 1 for MZ and 2 for DZ	#
#												#
#      Output: All bivariate statistics                                 #
#                                                                       #
#	  Notes: All comments with *** denote necessary input from user   #
#-----------------------------------------------------------------------#

## Loading Required/Useful Libraries ##
require(OpenMx)		# Loads OpenMX
require(psych)		# Loads Psych package
require(dplyr)
source("http://www.vipbg.vcu.edu/~vipbg/Tc24/GenEpiHelperFunctions.R") 

# Set default Mx optimizer to NPSOL
mxOption(NULL, "Default optimizer", "NPSOL")  

### Reading in Original Data ###

#*** Set working directory - location of data file and desired output
setwd("/home/jelman/netshare/K/Projects/SubcorticalChange")

#*** Read in data file
twins <- read.csv("./data/zV1V2Subcortical.csv",header=T)  

# Define phenotypic regions of interest
phenoNames = c("Cerebellum","Thalamus","Caudate","Putamen","Pallidum","Hippocampus","Amygdala",
                "Accumbens","LatVentricle","InfLatVentricle","Ventricle")
phenoVars = grep(paste(phenoNames, collapse="|"), names(twins), value=TRUE)

#*** Quick check that data file is properly read in
dim(twins)[1] #*** should be number of subjects
names(twins)[1:4] #*** should be subid, case, twin, zygosity
phenoVars #*** should be phenotypes

#*** Select variables to be included in loop and generate table of pairwise combinations
twins = twins %>% dplyr::select(1:4, phenoVars)
pairlist = t(combn(phenoNames, 2))
totPairs = nrow(pairlist) #*** should be number of variables to loop through, less one
totPairs

#*** Defining number of variables for each loop iteration
nv <- 8 # number of phenotypes
ntv  <- 2*nv

## Creating MZ and DZ data sets ##
# Dividing up twins
twinA <- twins[twins$twin=="A",]
twinB <- twins[twins$twin=="B",]
# Remerging by case to create paired data set
newtwins <- merge(twinA, twinB, by=c("case","zyg14"),all.x=TRUE, all.y=TRUE,suffixes=c("_A","_B"))
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

for (i in 1:totPairs){

  ## Defining Variables / Data for use with OpenMX ##    
  pheno1 <- pairlist[i,1]
  pheno2 <- pairlist[i,2]
  variables = c(paste0("adjL",pheno1), paste0("adjR",pheno1), paste0("adjL",pheno1,"_v2"), paste0("adjR",pheno1,"_v2"),
                paste0("adjL",pheno2), paste0("adjR",pheno2), paste0("adjL",pheno2,"_v2"), paste0("adjR",pheno2,"_v2"))
  selVars <- paste(variables,c(rep("_A",nv),rep("_B",nv)),sep="")
  
  # Print status
  cat("Starting: ", pheno1," ~ ", pheno2)

  # Making Sparse Data Frame with only phenotype 
  	MZdata <- MZdata1[,selVars]
  	DZdata <- DZdata1[,selVars]

  # Getting Number of Twins Used in Analysis
  	cpMZ <- sum(complete.cases(MZdata)) # Number of complete MZ pairs
  	sMZ <- dim(MZdata)[1]-cpMZ # Number of MZ singletons
  	cpDZ <- sum(complete.cases(DZdata)) # Number of complete DZ pairs
  	sDZ <- dim(DZdata)[1]-cpDZ # Number of DZ singletons
  

  # ----------------------------------------------------------------------
  # Model 1: ACE Cholesky Model 
  # ----------------------------------------------------------------------
  
  multiCholACEModel <- mxModel("multiCholACE",
       mxModel("ACE",
               # Matrices a, c, and e to store a, c, and e path coefficients
               mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=5, name="a", lbound=-3, ubound=3 ),
               mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.1, name="c", lbound=-3, ubound=3 ),
               mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=5, name="e", lbound=-3, ubound=3 ),
               # Matrices A, C, and E compute variance components
               mxAlgebra( expression=a %*% t(a), name="A" ),
               mxAlgebra( expression=c %*% t(c), name="C" ),
               mxAlgebra( expression=e %*% t(e), name="E" ),
               # Algebra to compute total variances and standard deviations (diagonal only)
               mxAlgebra( expression=A+C+E, name="V" ),
               mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
               mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),
               
               mxAlgebra( expression= V[1,1], name="V11" ),
               mxAlgebra( expression= V[2,2], name="V22" ),
               mxAlgebra( expression= V[3,3], name="V33" ),
               mxAlgebra( expression= V[4,4], name="V44" ),
               
               # Stuff for CIs on diagonals of A/V, C/V, E/V 
               mxAlgebra( expression=diag2vec(A/V),name="stndA"),
               mxAlgebra( expression=diag2vec(C/V),name="stndC"),
               mxAlgebra( expression=diag2vec(E/V),name="stndE"),
               mxCI(c("stndA","stndC","stndE")),
               #mxCI(c("V")),
               
               mxAlgebra((solve(sqrt(I*V)) %*% V %*% solve(sqrt(I*V))),name="Vcorr"),  
               mxAlgebra((solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A))),name="Acorr"),  
               mxAlgebra((solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C))),name="Ccorr"),  
               mxAlgebra((solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E))),name="Ecorr"),  
               mxCI(c("Vcorr","Acorr","Ccorr","Ecorr")),
               
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
  
  multiCholACEModel<- mxOption(multiCholACEModel, "Major iterations", 100000)
  
  multiCholACEFit <- mxRun(multiCholACEModel,intervals=FALSE)
  multiCholACESumm <- summary(multiCholACEFit)

  
  # ---------------------------------------------------------------------------------------------------------------
  # Model 2: Dual Change/Latent Change
  # Based off of McArdle and Nesselroade 1994 Latent Change Model
  # ---------------------------------------------------------------------------------------------------------------
  
  nf <- 4
  
  # Parameter Specifications from Level and Change Factors to Time1 & Time2 Factors
  f1true <- c(FALSE,FALSE,FALSE,FALSE,
              FALSE,FALSE,FALSE,FALSE,
              FALSE,FALSE,FALSE,FALSE,
              FALSE,FALSE,FALSE,FALSE)
  
  # Start Values for parameters from Level and Change Factors to Time1 & Time2 Factors
  f1load <- c(1,1,0,0,
              0,1,0,0,
              0,0,1,1,
              0,0,0,1)
  
  # Parameter Specifications from Time1 & Time2 Factors to Observed Variables
  f2true <- c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
              FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,
              FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,
              FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
              
  # Start Values for parameters from Time1 & Time2 Factors to Observed Variables
  f2load <- c(1,.7,0,0,0,0,0,0,
              0,0,1,.7,0,0,0,0,
              0,0,0,0,1,.7,0,0,
              0,0,0,0,0,0,1,.7)
  
  # Parameter Specifications for Unique Environmental Influences Specific to the Observed Variables 
  Erestrue <- c(TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,
                TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,
                TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,
                TRUE,FALSE,FALSE,FALSE,FALSE,
                TRUE,FALSE,TRUE,FALSE,
                TRUE,FALSE,TRUE,
                TRUE,FALSE,
                TRUE)
  
  # Start Values for Genetic and Common Environmental Influences Specific to the Observed Variables 
  Eresload <- c(.1,0,.1,0,0,0,0,0,
                .1,0,.1,0,0,0,0,
                .1,0,0,0,0,0,
                .1,0,0,0,0,
                .1,0,.1,0,
                .1,0,.1,
                .1,0,
                .1)
  
  DualChangeModel <- mxModel("DualChangeACE",
     mxModel("ACE",
             
                 ####################################################
                 # Matrices a2, c2, and e2 to store unique loadings for the test-specific Level and Change Factors
                 mxMatrix( type="Lower", nrow=nf, ncol=nf, free=TRUE, values=.6, name="a1", lbound=-3, ubound=3),
                 mxMatrix( type="Lower", nrow=nf, ncol=nf, free=TRUE, values=.1, name="c1", lbound=-3, ubound=3),
                 mxMatrix( type="Lower", nrow=nf, ncol=nf, free=TRUE, values=.4, name="e1", lbound=-3, ubound=3),
                 
                 # Matrices compute unique variance components for each test-level phenotype 
                 mxAlgebra( expression=(a1 %*% t(a1)), name="A1" ),
                 mxAlgebra( expression=(c1 %*% t(c1)), name="C1" ),
                 mxAlgebra( expression=(e1 %*% t(e1)), name="E1" ),
                 
                 mxAlgebra( expression= A1 + C1 + E1, name="V1" ),
                 
                 # Algebra to calculate standardized variance components that are unique to the test-specific factors
                 mxAlgebra( expression=diag2vec(A1/V1),name="stndA1"),
                 mxAlgebra( expression=diag2vec(C1/V1),name="stndC1"),
                 mxAlgebra( expression=diag2vec(E1/V1),name="stndE1"),
                 mxAlgebra( expression= V1[1,1], name="V122" ),
                 mxAlgebra( expression= V1[2,2], name="V122" ),
                 #mxCI(c("stndA1","stndE1","V122")),
                 mxCI(c("stndA1","stndC1","stndE1")),
             
                 # Algebra to calculate genetic and environmental correlations among the test-specific factors
                 mxMatrix( type="Iden", nrow=nf, ncol=nf, name="I1"),
                 mxAlgebra((solve(sqrt(I1*A1)) %*% A1 %*% solve(sqrt(I1*A1))),name="A1corr"),
                 mxAlgebra((solve(sqrt(I1*C1)) %*% C1 %*% solve(sqrt(I1*C1))),name="C1corr"),
                 mxAlgebra((solve(sqrt(I1*E1)) %*% E1 %*% solve(sqrt(I1*E1))),name="E1corr"),      
                 mxAlgebra((solve(sqrt(I1*V1)) %*% V1 %*% solve(sqrt(I1*V1))),name="V1corr"),     
                 #mxCI(c("A1corr","E1corr","V1corr")),
                 mxAlgebra( expression= A1corr[2,1], name="A1corr21" ),
                 mxAlgebra( expression= A1corr[4,2], name="A1corr42" ),
                 mxAlgebra( expression= E1corr[2,1], name="E1corr21" ),
                 mxAlgebra( expression= E1corr[4,2], name="E1corr42" ),
                 mxAlgebra( expression= V1corr[2,1], name="V1corr21" ),
                 mxAlgebra( expression= V1corr[4,2], name="V1corr42" ),
                 mxCI(c("A1corr21","A1corr42","E1corr21","E1corr42","V1corr21","V1corr42")),
             
                 
                 ####################################################
                 # Matrix f1 for factor loadings from overall latent phenotype to primary latent phenotypes
                 mxMatrix( type="Full", nrow=nf, ncol=nf, free=f1true, values=f1load, name="f1" ,lbound=-3, ubound=3),
                 
                 # Matrices A2, C2, and E2 compute variance components
                 mxAlgebra( expression=((f1 %*% (a1 %*% t(a1)) %*% t(f1)) ), name="A2" ),
                 mxAlgebra( expression=((f1 %*% (c1 %*% t(c1)) %*% t(f1)) ), name="C2" ),
                 mxAlgebra( expression=((f1 %*% (e1 %*% t(e1)) %*% t(f1)) ), name="E2" ),
                 
                 mxAlgebra( expression= A2 + C2 + E2, name="V2" ),
                 mxAlgebra( expression=diag2vec(A2/V2),name="stndA2"),
                 mxAlgebra( expression=diag2vec(C2/V2),name="stndC2"),
                 mxAlgebra( expression=diag2vec(E2/V2),name="stndE2"),
                 
                 # Matrix and Algebra for constraint on variance of overall latent phenotype
                 mxAlgebra( expression= A2 + C2 + E2, name="V2" ),
                 mxAlgebra( expression= diag2vec(V2), name="CovarFO" ),
                 mxMatrix( type="Full", nrow=nf, ncol=1, values=1, name="Unit"),
                 #mxConstraint(CovarFO == Unit),
                 
                 # Algebra to calculate genetic and environmental correlations among the test-specific factors
                 mxMatrix( type="Iden", nrow=nf, ncol=nf, name="I2"),
                 mxAlgebra((solve(sqrt(I2*A2)) %*% A2 %*% solve(sqrt(I2*A2))),name="A2corr"),
                 #mxAlgebra((solve(sqrt(I2*C2)) %*% C2 %*% solve(sqrt(I2*C2))),name="C2corr"),
                 mxAlgebra((solve(sqrt(I2*E2)) %*% E2 %*% solve(sqrt(I2*E2))),name="E2corr"),      
                 mxAlgebra((solve(sqrt(I2*V2)) %*% V2 %*% solve(sqrt(I2*V2))),name="V2corr"),     
                 #mxCI(c("A2corr","E2corr","V2corr")),
                 
                 # Matrix f2 for factor loadings from primary latent phenotypes to measured phenotypes
                 mxMatrix( type="Full", nrow=nv, ncol=nf, free=f2true, values=f2load, name="f2", lbound=-100, ubound=100),
                 #mxCI(c("f2")),
                 mxAlgebra( expression= f2[2,1], name="f221" ),
                 mxAlgebra( expression= f2[4,2], name="f242" ),
                 mxAlgebra( expression= f2[6,3], name="f263" ),
                 mxAlgebra( expression= f2[8,4], name="f284" ),
                 mxConstraint(f221 == f242),	# Contrain free factor loadings to be equal across time
                 mxConstraint(f263 == f284),	# Contrain free factor loadings to be equal across time
                 
                 # Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, name="as", lbound=-10, ubound=10),
                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=FALSE, values=0, name="cs", lbound=-10, ubound=10),
                 mxMatrix( type="Lower", nrow=nv, ncol=nv, free=Erestrue, values=Eresload, name="es", lbound=-3, ubound=3),
                 
                 # Matrices to calculate the overal variance components for each measured variable
                 mxAlgebra( expression=((f2 %*% f1 %*% (a1%*%t(a1)) %*% t(f1) %*% t(f2)) + (as %*% t(as))), name="A" ),
                 mxAlgebra( expression=((f2 %*% f1 %*% (c1%*%t(c1)) %*% t(f1) %*% t(f2)) + (cs %*% t(cs))), name="C" ),
                 mxAlgebra( expression=((f2 %*% f1 %*% (e1%*%t(e1)) %*% t(f1) %*% t(f2)) + (es %*% t(es))), name="E" ),
                 
                 # Matrices to calculate the SPECIFIC variance components for each measured variable
                 mxAlgebra( expression=(as %*% t(as)), name="As" ),
                 mxAlgebra( expression=(cs %*% t(cs)), name="Cs" ),
                 mxAlgebra( expression=(es %*% t(es)), name="Es" ),
                 
                 # Algebra to compute total variances and standard deviations (diagonal only)
                 mxAlgebra( expression=A+C+E, name="V" ),
                 mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I3"),
                 mxAlgebra( expression=solve(sqrt(I3*V)), name="iSD3"),
                 
                 mxAlgebra( expression=diag2vec(As/V),name="stndVCAs"),
                 mxAlgebra( expression=diag2vec(Cs/V),name="stndVCCs"),
                 mxAlgebra( expression=diag2vec(Es/V),name="stndVCEs"),
                 
                 mxAlgebra( expression=diag2vec(A/V),name="stndVCA"),
                 mxAlgebra( expression=diag2vec(C/V),name="stndVCC"),
                 mxAlgebra( expression=diag2vec(E/V),name="stndVCE"),
                 
                 #mxCI(c("stndVCAs")),
                 
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
     
     ##################################################
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
  
  DualChangeModel <- mxOption(DualChangeModel, "Major iterations", 100000)
  
  DualChangeFit <- mxRun(DualChangeModel , intervals=TRUE)
  DualChangeSumm <- summary(DualChangeFit)
  # DualChangeSumm 
  
  
  # tableFitStatistics(multiCholACEFit,DualChangeFit)
  


  #-------------------------------------------------------------------#
  ### OUTPUT ###
  #-------------------------------------------------------------------#
  
  # #ACE
  # ACEfit<- multiCholACEFit@output$Minus2LogLikelihood
  # AICace<- multiCholACESumm$AIC
  # BICace<- multiCholACESumm$BIC.Mx
  # ACEdof<- multiCholACESumm$degreesOfFreedom
  # 
  # # Dual Change
  # DCfit<- DualChangeFit@output$Minus2LogLikelihood
  # AICdc<- DualChangeSumm$AIC
  # BICdc<- DualChangeSumm$BIC.Mx
  # DCdof<- DualChangeSumm$degreesOfFreedom
  # DCpval<- round(pchisq(DCfit - ACEfit, DCdof - ACEdof,lower.tail=F),3)
  # IFAIL<- DualChangeFit@output$status[[1]]
  
  # # Dual change model output
  # DualChangeFit$ACE.stndA1@result
  # DualChangeFit$ACE.stndC1@result
  # DualChangeFit$ACE.stndE1@result
  # DualChangeFit$ACE.V1@result
  # DualChangeFit$ACE.A1corr@result
  # DualChangeFit$ACE.E1corr@result
  # DualChangeFit$ACE.V1corr@result
  # DualChangeFit$ACE.A1@result
  # DualChangeFit$ACE.C1@result
  # DualChangeFit$ACE.E1@result
  # 
  # DualChangeFit$ACE.f2@values
  # 
  # DualChangeFit$ACE.Es@result

  ##Combine Output
  # bivariate<-c(pheno1,pheno2,sat2ll,ACEfit,x2acefit,pacefit,AICace,
  # LP1_pheno1_A, LP1_pheno2_A, LP2_pheno2_A, LP1_pheno1_C, LP1_pheno2_C, LP2_pheno2_C, LP1_pheno1_E, LP1_pheno2_E, LP2_pheno2_E,
  # A1, A12, A2, C1, C12, C2, E1, E12, E2,
  # ACEfit_noAcorr,x2noAcorr,pnoAcorr,AICace_noAcorr,
  # ACEfit_noCcorr,x2noCcorr,pnoCcorr,AICace_noCcorr,ACEfit_noEcorr,x2noEcorr,pnoEcorr,AICace_noEcorr,
  # ACEfit_noACcorr,x2noACcorr,pnoACcorr,AICace_noACcorr,
  # ACEfit_noAEcorr,x2noAEcorr,pnoAEcorr,AICace_noAEcorr,
  # ACEfit_noCEcorr,x2noCEcorr,pnoCEcorr,AICace_noCEcorr,
  # ACEfit_nocorr,x2nocorr,pnocorr,AICace_nocorr)
  
  dualchangeOut<- DualChangeSumm$CI %>%  
                    dplyr::select(1:3) %>%  # Select estimate, lbound, ubound
                    mutate(parameter = gsub("\\[|\\]|,|ACE\\.", "", rownames(.))) %>%  # Cleanup and add parameter names to data.frame
                    gather(variable, value, -parameter) %>%  # Put into long format
                    unite(tmpcolname, parameter, variable) %>% # Join parameter and variable name into one column
                    arrange(tmpcolname) %>% 
                    mutate(Pheno1 = pheno1,
                           Pheno2 = pheno2, 
                           ACEfit = multiCholACEFit@output$Minus2LogLikelihood, # ACE model output
                           AICace = multiCholACESumm$AIC,
                           BICace = multiCholACESumm$BIC.Mx,
                           ACEdof = multiCholACESumm$degreesOfFreedom,
                           DCfit =DualChangeFit@output$Minus2LogLikelihood, # Dual Change model output
                           AICdc = DualChangeSumm$AIC,
                           BICdc = DualChangeSumm$BIC.Mx,
                           DCdof = DualChangeSumm$degreesOfFreedom,
                           DCpval = round(pchisq(DCfit - ACEfit, DCdof - ACEdof,lower.tail=F),3),
                           IFAIL =DualChangeFit@output$status[[1]]) %>%
                    spread(tmpcolname, value)
  
  output <- rbind(output, dualchangeOut)
  
}

## END OF LOOP ##

output <- as.data.frame(output)

write.csv(output, "./results/Subcort_DualChange_LOOP_out.csv", row.names=F)

