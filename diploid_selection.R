#Kaylee Rosenberger
#11/23/2022

#Modeling diploid selection 
#This script solves equations for one-locus diploid selection for varying initial proportions
#to determine at which proportions an allele representing resistance can become fixed in a 
#population

#Loading packages
library(ggplot2)

#######################################################################################
#Defining a function to solve the equation: 
diploid_selection = function(generations, p, q, fitness) {
    with(as.list(fitness), {
        for (t in 1:(generations - 1)){
            W = (p[t])^2*W_AA + 2*p[t]*q[t]*W_Aa + (q[t])^2*W_aa
            p[t+1] = (p[t])^2*(W_AA/W) + p[t]*q[t]*(W_Aa/W)
            q[t+1] = 1-p[t]
        }
        list(p)
    })
}

#######################################################################################
#Defining initial conditions, constants, and generations to simulate
#Defining many values to iterate over a lot of scenarios 
generations = 50
p_0 = seq(0, 0.5, 0.01)
q_0 = 1 - p_0
#Relative fitness 
W_AA = c(rep(0.1, 4), 0.2)
W_Aa = c(0.1, rep(0.05, 4))
W_aa = c(0.1, 0.05, 0.025, 0.01, 0.01)

#Container to store solutions in 
solutions = list()

#######################################################################################
#Iterating over scenarios and calculating solutions: 
ctr = 1
for(i in 1:length(p_0)) {
    for(j in 1:length(W_AA)) {
        fitness = list(W_AA = W_AA[j], W_Aa = W_Aa[j], W_aa = W_aa[j])
        solutions[[ctr]] = diploid_selection(generations, p_0[i], q_0[i], fitness)
        ctr = ctr + 1
    }
}

########################################################################################
#Plotting allele frequency dynamics for each scenario
for(i in 1:length(solutions)) {
    #Get data in format for plotting
    df = matrix(nrow = generations, ncol = 3)
    df[,1] = seq(1,generations,1)
    df[,2] = unlist(solutions[[i]])
    df[,3] = 1 - df[,2]
    colnames(df) = c("time", "p", "q")
    df = as.data.frame(df)
    
    #Plotting solution, saving plot as .pdf file 
    png(file=paste("Calculating Biological Quantities/course_project/figures/scenario", i, ".png", sep=""))
    print(ggplot(df, aes(x=time)) +
        geom_point(aes(y=p), color = 'green') +
        geom_point(aes(y=q), color = 'red') +
        ylim(0, 1) +
        theme_bw())
    dev.off()
}
