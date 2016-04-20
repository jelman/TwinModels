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
Vars      <- c('family','happy', 'life') 
nv        <- 3       # number of variables
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
svMe      <- c(7,5,5)            # start value for means
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
svMe      <- c(7,5, 5)                # start value for means
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


# ------------------------------------------------------------------------------                                           
# FIT GENETIC SUBMODELS

# ------------------------------------------------------------------------------                                           
# Fit Independent Pathway ACE Model
# ------------------------------------------------------------------------------                                           
nf        <- 1       # number of factors

# Matrices ac, cc, and ec to store a, c, and e path coefficients for common factors
pathAc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE,
                       values=.6, labels=labFull("ac",nv,nf), name="ac" )
pathCc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE,
                       values=.6, labels=labFull("cc",nv,nf), name="cc" )
pathEc    <- mxMatrix( type="Full", nrow=nv, ncol=nf, free=TRUE,
                       values=.6, labels=labFull("ec",nv,nf), name="ec" )

# Matrices as, cs, and es to store a, c, and e path coefficients for specific factors
pathAs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE,
                       values=4, labels=labDiag("as",nv), lbound=.00001, name="as" )
pathCs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE,
                       values=4, labels=labDiag("cs",nv), lbound=.00001, name="cs" )
pathEs    <- mxMatrix( type="Diag", nrow=nv, ncol=nv, free=TRUE,
                       values=5, labels=labDiag("es",nv), lbound=.00001, name="es" )

# Matrices A, C, and E compute variance components
covA      <- mxAlgebra( expression=ac %*% t(ac) + as %*% t(as), name="A" )
covC      <- mxAlgebra( expression=cc %*% t(cc) + cs %*% t(cs), name="C" )
covE      <- mxAlgebra( expression=ec %*% t(ec) + es %*% t(es), name="E" )

# Combine Groups
pars      <- list( pathAc, pathCc, pathEc, pathAs, pathCs, pathEs,
                   covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
minus2ll  <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj       <- mxAlgebraObjective( "m2LL" )
IndAceModel  <- mxModel( "IndACE", pars, modelMZ, modelDZ, minus2ll, obj )

# Run IndACE model
IndAceFit    <- mxRun(IndAceModel)

# Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices
IndACEpathMatrices <- c("iSD %*% ac","iSD %*% cc","iSD %*% ec","iSD %*% as","iSD %*% cs","iSD %*% es")
IndACEpathLabels <- c("stPathAc","stPathCc","stPathEc","stPathAs","stPathCs","stPathEs")
formatOutputMatrices(IndAceFit,IndACEpathMatrices,IndACEpathLabels,Vars,4)

# ------------------------------------------------------------------------------                                           
# Fit Independent Pathway AE Model
# ------------------------------------------------------------------------------                                           
IndAeModel   <- mxModel( IndAceFit, name="IndAE" )
IndAeModel   <- omxSetParameters( IndAeModel, labels=c(labFull("cc",nv,nf),labDiag("cs",nv)), free=FALSE, values=0 )
IndAeFit     <- mxRun(IndAeModel)
mxCompare(IndAceFit, IndAeFit)

# -------------------------------------------------------------
# Fit 2A (1C 1E) Factor - Independent Pathway Model
# -------------------------------------------------------------

# Change Dimension of Additive Genetic Factor Matrix Ac
nfA       <- 2
# Create Free and Values for 2 Additive Genetic Factors
#       free         values
#     A1  A2         A1  A2
# P1   T   T     P1  .5  .5
# P2   T   T     P2  .5  .5
# P3   T   T     P3  .5  .5
# P4   T   F     P4  .5   0
# P5   T   F     P5  .5   0
# P6   T   F     P6  .5   0
frAc2     <- c(T,T,T,T,T,T,  T,T,T,F,F,F)
svAc2     <- c(rep(.5,nv),  rep(.5,3),rep(0,3))

pathAc    <- mxMatrix( type="Full", nrow=nv, ncol=nfA, 
                       free=frAc2, values=svAc2, labels=labFull("ac",nv,nfA), name="ac" )

# Rebuild & Run Model
pars      <- list( pathAc, pathCc, pathEc, pathAs, pathCs, pathEs,
                   covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
Ind_2A_1C_1E_Model   <- mxModel( "Ind2ACE", pars, modelMZ, modelDZ, minus2ll, obj )
Ind_2A_1C_1E_Fit     <- mxRun(Ind_2A_1C_1E_Model)

IndNested <- list(IndAceFit, Ind_2A_1C_1E_Fit)
mxCompare(CholAceFit,IndNested)
formatOutputMatrices(Ind_2A_1C_1E_Fit,IndACEpathMatrices,IndACEpathLabels,Vars,4)

# -------------------------------------------------------------
# Fit 3A (1C 1E) Factor - Independent Pathway Model
# -------------------------------------------------------------
nfA <- 3
frAc3     <- c(T,T,T,T,T,T,  T,T,T,F,F,F,  F,F,F,T,T,T)
svAc3     <- c(rep(.5,nv),rep(.5,3),rep(0,3),rep(0,3),rep(.5,3))
pathAc    <- mxMatrix( type="Full", nrow=nv, ncol=nfA, 
                       free=frAc3, values=svAc3, labels=labFull("ac",nv,nfA), name="ac" )
pars      <- list( pathAc, pathCc, pathEc, pathAs, pathCs, pathEs,
                   covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
Ind_3A_1C_1E_Model   <- mxModel( "Ind3ACE", pars, modelMZ, modelDZ, minus2ll, obj )
Ind_3A_1C_1E_Fit     <- mxRun(Ind_3A_1C_1E_Model)

IndNested <- list(IndAceFit, Ind_2A_1C_1E_Fit, Ind_3A_1C_1E_Fit)
mxCompare(CholAceFit,IndNested)
formatOutputMatrices(Ind_3A_1C_1E_Fit,IndACEpathMatrices,IndACEpathLabels,Vars,4)

# -------------------------------------------------------------
# Fit 3A (0C 1E) Factor - Independent Pathway Model
# -------------------------------------------------------------
Ind_3A_0C_1E_Model   <- mxModel( Ind_3A_1C_1E_Fit, name="Ind_3A_0C_1E" )
Ind_3A_0C_1E_Model   <- omxSetParameters( Ind_3A_1C_1E_Model,
                      labels=c(labFull("cc",nv,nf),labDiag("cs",nv)), free=FALSE, values=0 )
Ind_3A_0C_1E_Fit     <- mxRun(Ind_3A_0C_1E_Model)
mxCompare(Ind_3A_1C_1E_Fit, Ind_3A_0C_1E_Fit)

# -------------------------------------------------------------
# Fit 3A (0C 3E) Factor - Independent Pathway Model
# -------------------------------------------------------------
nfE <- 3
frEc3     <- c(T,T,T,T,T,T,  T,T,T,F,F,F,  F,F,F,T,T,T)
svEc3     <- c(rep(.5,nv),rep(.5,3),rep(0,3),rep(0,3),rep(.5,3))
pathEc    <- mxMatrix( type="Full", nrow=nv, ncol=nfE, 
                       free=frEc3, values=svEc3, labels=labFull("ec",nv,nfE), name="ec" )
pars      <- list( pathAc, pathCc, pathEc, pathAs, pathCs, pathEs,
                   covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
Ind_3A_0C_3E_Model   <- mxModel( "Ind3AE", pars, modelMZ, modelDZ, minus2ll, obj )
Ind_3A_0C_3E_Fit     <- mxRun(Ind_3A_0C_3E_Model)

IndNested <- list(IndAceFit, Ind_2A_1C_1E_Fit, Ind_3A_1C_1E_Fit, Ind_3A_0C_3E_Fit)
mxCompare(CholAceFit,IndNested)
formatOutputMatrices(Ind_3A_0C_3E_Fit,IndACEpathMatrices,IndACEpathLabels,Vars,4)


# ------------------------------------------------------------------------------                                           
# Fit Common Pathway ACE Model
# ------------------------------------------------------------------------------                                           
nl        <- 1

# Matrices ac, cc, and ec to store a, c, and e path coefficients for latent phenotype(s)
pathAl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.6,
                       labels=labLower("al",nl), lbound=.00001, name="al" )
pathCl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.6,
                       labels=labLower("cl",nl), lbound=.00001, name="cl" )
pathEl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.6,
                       labels=labLower("el",nl), lbound=.00001, name="el" )

# Matrix and Algebra for constraint on variance of latent phenotype
covarLP   <- mxAlgebra( expression=
                       al %*% t(al) + cl %*% t(cl) + el %*% t(el), name="CovarLP" )
varLP     <- mxAlgebra( expression= diag2vec(CovarLP), name="VarLP" )
unit      <- mxMatrix( type="Unit", nrow=nl, ncol=1, name="Unit")
varLP1    <- mxConstraint( expression=VarLP == Unit, name="varLP1")

# Matrix f for factor loadings on latent phenotype
pathF     <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=TRUE, values=.2,
                       labels=labFull("fl",nv,nl), name="fl" )

# Matrices A, C, and E compute variance components
covA      <- mxAlgebra( expression=fl %&% (al %*% t(al)) + as %*% t(as), name="A" )
covC      <- mxAlgebra( expression=fl %&% (cl %*% t(cl)) + cs %*% t(cs), name="C" )
covE      <- mxAlgebra( expression=fl %&% (el %*% t(el)) + es %*% t(es), name="E" )

# Combine Groups
pars      <- list( pathAl, pathCl, pathEl, covarLP, varLP, unit, pathF,
                   pathAs, pathCs, pathEs, covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
minus2ll  <- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj       <- mxAlgebraObjective( "m2LL" )
ComAceModel  <- mxModel( "ComACE", pars, varLP1, modelMZ, modelDZ, minus2ll, obj )

# Run CholACE model
ComAceFit     <- mxRun(ComAceModel)

# Generate List of Parameter Estimates and Derived Quantities using formatOutputMatrices
ComACEpathMatrices <- c("al","cl","el","iSD %*% fl","iSD %*% as","iSD %*% cs","iSD %*% es")
ComACEpathLabels <- c("stPathAl","stPathCl","stPathEl","stPathFl","stPathAs","stPathCs","stPathEs")
formatOutputMatrices(ComAceFit,ComACEpathMatrices,ComACEpathLabels,Vars,4)

# ------------------------------------------------------------------------------                                           
# Fit Common Pathway ACE Model
# ------------------------------------------------------------------------------                                           
nl      <- 2
pathAl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.6,
                       labels=labLower("al",nl), lbound=.00001, name="al" )
pathCl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.6,
                       labels=labLower("cl",nl), lbound=.00001, name="cl" )
pathEl    <- mxMatrix( type="Lower", nrow=nl, ncol=nl, free=TRUE, values=.6,
                       labels=labLower("el",nl), lbound=.00001, name="el" )
unit      <- mxMatrix( type="Unit", nrow=nl, ncol=1, name="Unit")
varLP1    <- mxConstraint( expression=VarLP == Unit, name="varLP1")
pathF     <- mxMatrix( type="Full", nrow=nv, ncol=nl, free=TRUE, values=.2,
                       labels=labFull("fl",nv,nl), name="fl" )
pars      <- list( pathAl, pathCl, pathEl, covarLP, varLP, unit, pathF,
                   pathAs, pathCs, pathEs, covA, covC, covE, covP, matI, invSD, meanG )
modelMZ   <- mxModel( pars, covMZ, dataMZ, objMZ, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, dataDZ, objDZ, name="DZ" )
Com_2L_AceModel  <- mxModel( "ComACE_2L", pars, varLP1, modelMZ, modelDZ, minus2ll, obj )
Com_2L_AceFit     <- mxRun(Com_2L_AceModel)

ComNested  <- list(ComAceFit,Com_2L_AceFit)
mxCompare(CholAceFit,ComNested)
formatOutputMatrices(ComAceFit,ComACEpathMatrices,ComACEpathLabels,Vars,4)

# ------------------------------------------------------------------------------                                           
                           
# Print Comparative Fit Statistics
CholAceNested <- c( IndNested, ComNested)
mxCompare(CholAceFit,CholAceNested)
