
devtools::load_all("C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/coxed/cox_epistasis")


num.covariates <- 22
num.inst <- 600
# Lets try user defined beta values (coefficients)
factor <- 25

coefficients <- rep(0.1,num.covariates)
coefficients[3:4] <- coefficients[3:4] * factor



# specify the interactions matrix as a [1:num.inst, 1:num.covariates] matrix
# Initialize the matrix with zeros
inter.mat <- matrix(0, nrow=num.covariates, ncol=num.covariates)

# Set the elements at [1,2] and [2,1] to 1
inter.mat[1,2] <- 1
inter.mat[2,1] <- 1


# Make T bigger and decrease the variance of the X variables
# N = 200: number of observations
# T = 5: maximum time
# xvars = 19: number of covariates
# censor = 0.2: censoring rate (proportion of observations that are censored)
# num.data.frames = 1: number of datasets to simulate

# For interaction effects
#simdata <- coxed::sim.survdata(N=num.inst, T=40, xvars = num.covariates, censor = 0.1, num.data.frames=1, beta=coefficients/2, interactions=TRUE, inter.mat=inter.mat, mu=0, sd=0.5)

# For no interaction effects
#simdata <- coxed::sim.survdata(N=num.inst, T=40, xvars = num.covariates, censor = 0.1, num.data.frames=1, beta=coefficients/2, interactions=FALSE, inter.mat=inter.mat, mu=0, sd=0.5)

#write.csv(simdata$data, "C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData/simdata.csv")


# Create multiple datasets
num.datasets <- 10
for (i in 1:num.datasets) {
    simdata <- coxed::sim.survdata(N=num.inst, T=40, xvars = num.covariates, censor = 0.1, num.data.frames=1, beta=coefficients/2, interactions=FALSE, inter.mat=inter.mat, mu=0, sd=0.5)
    file.name <- paste0("C:/Users/patri/Dropbox/Ph.D/Research/Doctoral_Research/Bioinformatics/sNPDR/data/simulatedData/simdata_", i, ".csv")
    write.csv(simdata$data, file.name)
}





