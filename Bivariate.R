# ------------------------------------------------------------------------------
# Program: twinMulAceCon.R  
#  Author: Hermine Maes
#    Date: 03 04 2014 
#
# Multivariate Twin Analysis model to estimate causes of variation
# Matrix style model - Raw data - Continuous data
# -------|---------|---------|---------|---------|---------|---------|


# Load Library
require(OpenMx)
require(psych)
source("myFunctions.R")
source("GenEpiHelperFunctions.R")

# --------------------------------------------------------------------
# PREPARE DATA

# Load Data
nl <- read.table ("DHBQ_bs.dat", header=T, na=-999)
describe(nl, skew=F)

# Recode Data for Analysis - Rescale variables to have variances around 1.0
nl$family1  <- nl$gff1/5
nl$family2  <- nl$gff2/5
nl$happy1   <- nl$hap1/4
nl$happy2   <- nl$hap2/4
nl$life1    <- nl$sat1/5
nl$life2    <- nl$sat2/5
nl$anxdep1  <- nl$AD1/4
nl$anxdep2  <- nl$AD2/4
nl$somatic1 <- nl$SOMA1/2
nl$somatic2 <- nl$SOMA2/2
nl$social1  <- nl$SOC1/1.5
nl$social2  <- nl$SOC2/1.5

# Select Variables for Analysis
Vars      <- c('family','happy') 
nv        <- 2       # number of variables
ntv       <- nv*2    # number of total variables
selVars   <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep="")

# Select random subset to reduce time to fit examples
testData <- head(nl,n=500)

# Select Data for Analysis
mzData    <- subset(testData, zyg2==1, selVars)
dzData    <- subset(testData, zyg2==2, selVars)
describe(mzData, skew=F)
describe(dzData, skew=F)
dim(mzData)
dim(dzData)

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
cor(mzData,use="complete")
cor(dzData,use="complete")


# ------------------------------------------------------------------------------
# PREPARE SATURATED MODEL

# Saturated Model
# Set Starting Values
svMe      <- c(7,5)            # start value for means
svVa      <- valDiag(ntv,1.0)          # start values for variances
lbVa      <- valODiag(ntv,.0001,-10)   # lower bounds for covariances

# -------|---------|---------|---------|---------|---------|---------|
# Algebra for expected Mean Matrices in MZ & DZ twins
meanMZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                       values=svMe, labels=labFull("meMZ",1,ntv), name="expMeanMZ" )
meanDZ    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                       values=svMe, labels=labFull("meDZ",1,ntv), name="expMeanDZ" )

# Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covMZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE,
                       values=svVa, lbound=lbVa, labels=labSymm("vaMZ",ntv), name="expCovMZ" )
covDZ     <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE,
                       values=svVa, lbound=lbVa, labels=labSymm("vaDZ",ntv), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ     <- mxFIMLObjective( covariance="expCovMZ", means="expMeanMZ", dimnames=selVars )
objDZ     <- mxFIMLObjective( covariance="expCovDZ", means="expMeanDZ", dimnames=selVars )

# Combine Groups
modelMZ   <- mxModel( "MZ", meanMZ, covMZ, dataMZ, objMZ )
modelDZ   <- mxModel( "DZ", meanDZ, covDZ, dataDZ, objDZ )
minus2ll  <- mxAlgebra( MZ.objective+ DZ.objective, name="minus2sumloglikelihood" )
obj       <- mxAlgebraObjective( "minus2sumloglikelihood" )
ciCov     <- mxCI( c('MZ.expCovMZ','DZ.expCovDZ' ))
ciMean    <- mxCI( c('MZ.expMeanMZ','DZ.expMeanDZ' ))
twinSatModel   <- mxModel( "twinSat", modelMZ, modelDZ, minus2ll, obj, ciCov, ciMean )

# ------------------------------------------------------------------------------
# RUN SATURATED MODEL

# Run Saturated Model
twinSatFit     <- mxRun( twinSatModel, intervals=F )
twinSatSumm    <- summary( twinSatFit)

# Generate Saturated Model Output
round(twinSatFit$MZ.expMeanMZ@values,2)
round(twinSatFit$DZ.expMeanDZ@values,2)
round(twinSatFit$MZ.expCovMZ@values,2)
round(twinSatFit$DZ.expCovDZ@values,2)

# ------------------------------------------------------------------------------
# RUN SATURATED SUBMODELS

# Constrain expected Means and Variances to be equal across twin order
eqMeVaTwinModel    <- mxModel( twinSatModel, name="eqM&Vtwins" )
eqMeVaTwinModel    <- omxSetParameters( eqMeVaTwinModel, label=labFull("meMZ",1,ntv), 
                       free=TRUE, values=svMe, newlabels=labFull("meMZ",1,nv) )
eqMeVaTwinModel    <- omxSetParameters( eqMeVaTwinModel, label=labFull("meDZ",1,ntv),
                       free=TRUE, values=svMe, newlabels=labFull("meDZ",1,nv) )
eqMeVaTwinModel    <- omxSetParameters( eqMeVaTwinModel, label=labDiag("vaMZ",ntv),
                       free=TRUE, values=diag(svVa), newlabels=labDiag("vaMZ",nv)  )
eqMeVaTwinModel    <- omxSetParameters( eqMeVaTwinModel, label=labDiag("vaDZ",ntv),
                       free=TRUE, values=diag(svVa), newlabels=labDiag("vaDZ",nv) )

eqMeVaTwinFit      <- mxRun( eqMeVaTwinModel, intervals=F )
eqMeVaTwinSumm     <- summary( eqMeVaTwinFit )
mxCompare(twinSatFit, eqMeVaTwinFit)

# Constrain expected Means and Variances to be equal across twin order and zygosity
eqMeVaZygModel     <- mxModel( eqMeVaTwinModel, name="eqM&Vzyg" )
eqMeVaZygModel     <- omxSetParameters( eqMeVaZygModel, labels=c(labFull("meMZ",1,nv),labFull("meDZ",1,nv)),
                       free=TRUE, values=svMe, newlabels=labFull("me",1,nv) )
eqMeVaZygModel     <- omxSetParameters( eqMeVaZygModel, labels=c(labDiag("vaMZ",nv),labDiag("vaDZ",nv)),
                       free=TRUE, values=diag(svVa), newlabels=labDiag("va",nv) )

eqMeVaZygFit       <- mxRun( eqMeVaZygModel, intervals=F )
eqMeVaZygSumm      <- summary( eqMeVaZygFit )
mxCompare(eqMeVaTwinFit, eqMeVaZygFit)
                       
# Print Comparative Fit Statistics
SatNested <- list(eqMeVaTwinFit, eqMeVaZygFit)
mxCompare(twinSatFit, SatNested)

# ------------------------------------------------------------------------------
# PREPARE GENETIC MODEL

# ------------------------------------------------------------------------------
# Cholesky Decomposition ACE Model
# ------------------------------------------------------------------------------
# Set Starting Values
svMe      <- c(7,5)                # start value for means
svPa      <- valDiag(nv,.6)                # start values for parameters on diagonal
lbPa      <- valLUDiag(nv,.0001,-10,NA)    # lower bounds for parameters on diagonal

# Matrices declared to store a, c, and e Path Coefficients
pathA     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE,
                       values=svPa, labels=labLower("a",nv), lbound=lbPa, name="a" )
pathC     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE,
                       values=svPa, labels=labLower("c",nv), lbound=lbPa, name="c" )
pathE     <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE,
                       values=svPa, labels=labLower("e",nv), lbound=lbPa, name="e" )
	
# Matrices generated to hold A, C, and E computed Variance Components
covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances and standard deviations (diagonal only)
covP      <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE,
                       values=svMe, labels=labFull("me",1,nv), name="expMean" )
covMZ     <- mxAlgebra( expression= rbind( cbind(V, A+C), cbind(A+C, V)), name="expCovMZ" )
covDZ     <- mxAlgebra( expression= rbind( cbind(V, 0.5%x%A+C), cbind(0.5%x%A+C, V)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ     <- mxFIMLObjective( covariance="expCovMZ", means="expMean", dimnames=selVars )
objDZ     <- mxFIMLObjective( covariance="expCovDZ", means="expMean", dimnames=selVars )

# Combine Groups
pars      <- list( pathA, pathC, pathE, covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
minus2ll  <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj       <- mxAlgebraObjective( "m2LL" )
CholAceModel  <- mxModel( "CholACE", pars, modelMZ, modelDZ, minus2ll, obj )

# ------------------------------------------------------------------------------
# RUN GENETIC MODEL

# Run Cholesky Decomposition ACE model
CholAceFit    <- mxRun(CholAceModel)

mxCompare(twinSatFit,CholAceFit)
mxCompare(eqMeVaZygFit,CholAceFit)

# Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices
# ACE  Standardized Path Coefficients (pre-multiplied by inverse of standard deviations)
CholACEpathMatrices <- c("iSD %*% a","iSD %*% c","iSD %*% e")
CholACEpathLabels <- c("stPathA","stPathC","stPathE")
formatOutputMatrices(CholAceFit,CholACEpathMatrices,CholACEpathLabels,Vars,4)

# ACE Covariance Matrices & Proportions of Variance Matrices
ACEcovMatrices <- c("A","C","E","V","A/V","C/V","E/V")
ACEcovLabels <- c("covA","covC","covE","Var","stCovA","stCovC","stCovE")
formatOutputMatrices(CholAceFit,ACEcovMatrices,ACEcovLabels,Vars,4)

# ACE Correlation Matrices 
ACEcorMatrices <- c("solve(sqrt(I*A)) %&% A","solve(sqrt(I*C)) %&% C","solve(sqrt(I*E)) %&% E")
ACEcorLabels <- c("corA","corC","corE")
formatOutputMatrices(CholAceFit,ACEcorMatrices,ACEcorLabels,Vars,4)

# ------------------------------------------------------------------------------
# Cholesky Decomposition AE Model
# ------------------------------------------------------------------------------
CholAeModel   <- mxModel( CholAceFit, name="CholAE" )
CholAeModel   <- omxSetParameters( CholAeModel, labels=labLower("c",nv), free=FALSE, values=0 )
CholAeFit     <- mxRun(CholAeModel)
mxCompare(CholAceFit, CholAeFit)


# Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices
# ACE  Standardized Path Coefficients (pre-multiplied by inverse of standard deviations)
CholAEpathMatrices <- c("iSD %*% a","iSD %*% c","iSD %*% e")
CholAEpathLabels <- c("stPathA","stPathC","stPathE")
formatOutputMatrices(CholAeFit,CholAEpathMatrices,CholAEpathLabels,Vars,4)

# ACE Covariance Matrices & Proportions of Variance Matrices
AEcovMatrices <- c("A","C","E","V","A/V","C/V","E/V")
AEcovLabels <- c("covA","covC","covE","Var","stCovA","stCovC","stCovE")
formatOutputMatrices(CholAeFit,AEcovMatrices,AEcovLabels,Vars,4)

# ACE Correlation Matrices 
AEcorMatrices <- c("solve(sqrt(I*A)) %&% A","solve(sqrt(I*E)) %&% E")
AEcorLabels <- c("corA","corE")
formatOutputMatrices(CholAeFit,AEcorMatrices,AEcorLabels,Vars,4)


