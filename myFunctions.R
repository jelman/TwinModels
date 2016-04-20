fitGofs   <- function(fit) {
            summ <- summary(fit)
            cat(paste("Mx:", fit@name,
                      "  os=", summ$observedStatistics,              # number of observed statistics = non-missing data points
                      "  no=", summ$numObs,                          # number of observations = number of records
                      "  ep=", summ$estimatedParameters,             # number of estimated parameters
                      "  df=", summ$defreesOfFreedom,                # degrees of freedom
                      "  ll=", round(summ$Minus2LogLikelihood,4),    # -2 log likelihood of data
                      "  cpu=", round(summ$cpu,4),                   # cpu time to run model
                      "  ifail=", fit@output$status[[1]],            # status of model fitting
                      "\n",sep=""))
}



fitGofs     <- function(fit) {
            summ <- summary(fit)
            cat(paste("Mx:", fit@name,"  os=", summ$ob,"  ns=", summ$nu,
                      "  ep=", summ$es,"  df=", summ$de, "  ll=", round(summ$Mi,4), 
                      "  cpu=", round(summ$cpu,4), "  ifail=", fit@output$status[[1]], "\n",sep=""))
}


fitGFPs     <- function(fit) {
            summ <- summary(fit)
            output <- as.data.frame(t(c(fit@name, summ$ob, summ$nu, summ$es, summ$de, 
            round(summ$Mi,4), round(summ$AIC,2), round(summ$cpu,4), fit@output$status[[1]] )))
            names(output) <- c("Mx","os","no","ep","df","ll","AIC","cpu","ifail")
            print(output)
            write.table(output,file=paste(filename,".txt",sep=""),sep=",",row.names=FALSE,col.names=FALSE,append=T)
}

fitPars     <- function(fit) {
          print(round(fit@output$estimate, 4))
}


testModel <- function(model, func, nfit) {
#		model  <- mxOption( model, "Function precision", func )
	        fit   <- mxRun(model)
	        fitGofs(fit)
	        fitGFPs(fit)
#	        for (i in 1:nfit) {
#	        fit   <- mxRun(fit)
#	        fitGofs(fit) }
	        return(fit)
}



# Functions to assign labels
labLower   <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="_") }
labSdiag   <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:(nv-1))),rep(1:(nv-1),(nv-1):1),sep="_") }
labFullSq  <- function(lab,nv) { paste(lab,1:nv,rep(1:nv,each=nv),sep="_") }
labDiag    <- function(lab,nv) { paste(lab,1:nv,1:nv,sep="_") } 
labSymm    <- function(lab,nv) { paste(lab,rev(nv+1-sequence(1:nv)),rep(1:nv,nv:1),sep="_") }
labFull    <- function(lab,nr,nc) { paste(lab,1:nr,rep(1:nc,each=nr),sep="_") }
labFull2    <- function(lab,br,nr,bc,nc) { paste(lab,br:nr,rep(bc:nc,each=nr),sep="_") }
labDiag2    <- function(lab,bv,nv) { paste(lab,bv:nv,bv:nv,sep="_") } 

# Functions to assign values
valDiag  <- function(dim, valD) {
valF    <- diag(valD,dim,dim)        # values for diagonal of covariance matrix
valF
}
valODiag  <- function(dim, valD, valOD) {
valF    <- diag(valD,dim,dim)        # values for diagonal of covariance matrix
valF[lower.tri(valF)] <- valOD       # values for below diagonal elements 
valF[upper.tri(valF)] <- valOD       # values for above diagonal elements
valF
}
valLUDiag  <- function(dim, valD, valLD, valUD) {
valF    <- diag(valD,dim,dim)        # values for diagonal of covariance matrix
valF[lower.tri(valF)] <- valLD       # values for below diagonal elements 
valF[upper.tri(valF)] <- valUD       # values for above diagonal elements
valF
}