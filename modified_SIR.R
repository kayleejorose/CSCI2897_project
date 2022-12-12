#Kaylee Rosenberger
#11/23/2022

#Modified SIR model representing the dynamics of DED spread and resistance spread
#This script solves equations based on my modified SIR model 

#######################################################################################
#Loading packages
library(deSolve)
library(ggplot2)

#Defining parameters 
b = 0.15 #birth rate
d = 0.1 #natural death rate
B = 0.5 #infection rate
D = 0.5 #death rate due to DED 
f = 0.5 #scaling factor for death due to DED 
w = 1 #factor for biasing the births of each category to themselves
N = c(s=0.7, is=0.01, to=0.25, it=0.01, r=0.03) #initial population proportions 

#Defining the equations in this function 
str_model = function(t, N, params) {
    with(as.list(N), {
        dS = b*(w*s/((w*s)+to))*s + b*(s/(s+(w*to)+r))*to - d*s - B*s*(is+it)
        dIs = B*s*(is+it) - D*is 
        dT = b*(to/((w*s)+to))*s + b*(w*to/(s+(w*to)))*to + b*(to/(to+(w*r)))*r - d*to - B*to*(is+it)
        dIt = B*to*(is+it) - f*D*it 
        dR = b*(r/((w*to)+r))*to + b*(w*r/(to+(w*r)))*r - d*r
        list(c(dS, dIs, dT, dIt, dR))
    })
}

#defining times to calculate solution
times = seq(1, 150, 0.2)

#calculating solutions 
sol = ode(y=N, times = times, func = str_model)

matplot(x = sol[,1], y = sol[,-1], type = "l", lwd = 2, lty = c("solid"), col = c("red", "orange", "yellow", "green", "blue"), xlab = "Time", ylab = "Proportion of population", ylim = c(0,1))
title("Population disease dynamics")
abline(0,0, lty=2)

#Here, S=red, Is=orange, t=yellow, =green, and r=blue (rainbow order)
#Currently, the proportions of each category don't add to 1
#It's because b > d for R, so over time, it keeps growing. Which isn't necessarily wrong... but it should stop at some point 
#When b=d then R stays constant the entire time which also isn't what we want. 
#could it be because the infected individuals are not contriubting to the populations?
