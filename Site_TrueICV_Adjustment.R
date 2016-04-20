#-----------------------------------------------------------------------#
# Program: Quick Loop to Adjust for age and SITE			                  #
#  Author: Kelly Spoon                                                  #
#    Date: 04 11 2012                                                   #
#                                                                       #
# Description: Simple linear regression storing residuals			          #
#     	   							                                                #
#-----------------------------------------------------------------------#

### Reading in Original Data ###
#*** Set working directory - location of data file and desired output
	setwd("C:/Documents and Settings/mspanizzon/Desktop/Current Projects/V1DTI_Subcortical MD")

#*** Read in data file
	data <- read.csv("V1_ASEG.csv",header=T)

#*** Read variable names from data and store in a vector
	names1 <- names(data)

#*** Grab only variable names with *CTX* as variables of interest
#	index <- which(substr(names1,2,4)=="CTX") #storing location of variables with correct 2nd 3rd and 4th letters
#	variables <- names1[index] #grabbing a list of those variable names

#*** Alternate variable selection, commented out, but here if you want it
	variables <- names1[8:35] # you can select start and end point for the grabbing of variables

#*** Check variables are correct
	variables
	V <- length(variables)	

### Creating Storate Data Frame ###

# Set number of individuals 
	n <- dim(data)[1]
	tot <- dim(data)[2]

# Create Data Frame
	data <- cbind(data,matrix(NA,nrow=n,ncol=V))
	names(data) <- c(names1,paste(variables,"adj",sep=""))

### Running Loop Using lapply ###

# fitting models
models <- lapply(variables, function(x) {
    lm(substitute(i ~ as.factor(SITE) + age + IntracranialVol, list(i = as.name(x))), data = data, na.action=na.exclude)
})

# storing residuals from each model into data frame
for(v in 1:V){
data[,tot+v] <- residuals(models[[v]])
}

#dataR is now your residualized parameters

for(v in 1:V){
data[,tot+v] <- scale(data[tot+v])
}

dataR <- data[,c(1:28,(tot+1):(tot+V))]
write.csv(dataR,"V1DTI_MD_AgeSiteICVAdj.csv")