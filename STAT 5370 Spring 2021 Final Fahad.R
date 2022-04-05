# Read in the data.
alldat = read.table("~/Downloads/STAT_5370__Final_Exam_-_Take_Home_portion/bodytemp.dat", quote="\"", comment.char="")

names(alldat) = c("Temp","Sex","Heartrate")

# To make life easy, give the variables simple names.
x1 = alldat$Heartrate
x2 = alldat$Sex - 1
x3 = x1*x2
y = alldat$Temp

# Some basic variables used throughout.
n = length(y)
p = 3

# Build y.
y = matrix(y,nrow=n,ncol=1)

# Build X.
X = matrix(nrow=n,ncol=4)
X[,1] = rep(1,n)
X[,2] = x1
X[,3] = x2
X[,4] = x3

# Bayesian Prior Parameters

b0 = matrix(ncol=1,nrow=4)
b0[1,1] = 98.6  # beta0 prior mean
b0[2,1] = 0     # beta1 prior mean
b0[3,1] = 0     # beta2 prior mean
b0[4,1] = 0     # beta3 prior mean

B0 = matrix(0,ncol=4,nrow=4)
B0[1,1] = 2   # beta0 prior variance/sigma^2
B0[2,2] = 0.1 # beta1 prior variance/sigma^2
B0[3,3] = 3   # beta2 prior variance/sigma^2
B0[4,4] = 0.1 # beta3 prior variance/sigma^2

n0 = 30    # Prior "sample size"
S0 = 0.3  # Prior "expected value" (approximately) of sigma^2

#legend("topleft", legend=c("Male", "Female", "Male Regression", "Male 95% HPD/PI", "Female Regression", "Female 95% HPD/PI"), lty=c(NA,NA,1,2,1,2), pch=c(1,1,NA,NA,NA,NA), col=c("black","red","black","black","red","red"))


# Find posterior parameters.
B1 = solve(solve(B0)+t(X)%*%X)
b1 = B1 %*% (solve(B0)%*%b0+t(X)%*%y)
n1 = n0 + n
S1 = as.numeric((n0*S0 + (t(b0)%*%solve(B0)%*%b0 + t(y)%*%y - t(b1)%*%solve(B1)%*%b1))/n1)
# print posterior parameters
b1
B1
n1
S1
# Posterior variance/covariance for beta.
bvcovbeta = S1*B1

btcrit = qt(0.975,df=n1)
cat("Posterior mean of beta0 = ", b1[1], " (sd = ", sqrt(bvcovbeta[1,1]), "), 95% HPD = (", b1[1]-btcrit*sqrt(bvcovbeta[1,1]), ", ", b1[1]+btcrit*sqrt(bvcovbeta[1,1]), ")\n", sep="")
cat("Posterior mean of beta1 = ", b1[2], " (sd = ", sqrt(bvcovbeta[2,2]), "), 95% HPD = (", b1[2]-btcrit*sqrt(bvcovbeta[2,2]), ", ", b1[2]+btcrit*sqrt(bvcovbeta[2,2]), ")\n", sep="")
cat("Posterior mean of beta3 = ", b1[3], " (sd = ", sqrt(bvcovbeta[3,3]), "), 95% HPD = (", b1[3]-btcrit*sqrt(bvcovbeta[3,3]), ", ", b1[3]+btcrit*sqrt(bvcovbeta[3,3]), ")\n", sep="")
cat("Posterior mean of beta4 = ", b1[4], " (sd = ", sqrt(bvcovbeta[4,4]), "), 95% HPD = (", b1[4]-btcrit*sqrt(bvcovbeta[4,4]), ", ", b1[4]+btcrit*sqrt(bvcovbeta[4,4]), ")\n", sep="")

cat("Posterior mean of sigma^2 (approx) = ", S1, "\n", sep="")

# subset the data into separate tables/ splitting into male and female
males = alldat[x2==0,]
females = alldat[x2==1,]

plot(males$Heartrate, males$Temp, xlab = "heart rate" ,  ylab = "body temperature", xlim= c(50,90), ylim = c(95,103))
points(females$Heartrate, females$Temp, col="red")

#legend("topleft", legend=c("Male", "Female", "Male Regression", "Male 95% HPD/PI", "Female Regression", "Female 95% HPD/PI"), lty=c(NA,NA,1,2,1,2), pch=c(1,1,NA,NA,NA,NA), col=c("black","red","black","black","red","red"))

# Plot Bayesian stuff now.  
# For males
x.plot_males = seq(min(x1), max(x1), 0.1)
y.plot_males = b1[1] + b1[2]*x.plot_males
lines(x.plot_males, y.plot_males, lty=1, col="black")

lci = vector(mode="numeric", length=length(x.plot_males))
uci = vector(mode="numeric", length=length(x.plot_males))
#lpi = vector(mode="numeric", length=length(x.plot_males))
#upi = vector(mode="numeric", length=length(x.plot_males))
for (i in 1:length(x.plot_males))
{
  x.star = matrix(c(1,x.plot_males[i],0,0*x.plot_males[i]), nrow=1, ncol=4)
  lci[i] = as.numeric(x.star%*%b1 - btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
  uci[i] = as.numeric(x.star%*%b1 + btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
  #lpi[i] = as.numeric(x.star%*%b1 - btcrit*sqrt(S1)*sqrt(1+x.star%*%B1%*%t(x.star)))
  #upi[i] = as.numeric(x.star%*%b1 + btcrit*sqrt(S1)*sqrt(1+x.star%*%B1%*%t(x.star)))
}

lines(x.plot_males, lci, lty=2, col="black")
lines(x.plot_males, uci, lty=2, col="black")
#lines(x.plot, lpi, lty=2, col="red")
#lines(x.plot, upi, lty=2, col="red")

# For Females
x.plot_females = seq(min(x1), max(x1), 0.1)
y.plot_females = b1[1] +b1[3]+ (b1[2]+b1[4])*x.plot_females
lines(x.plot_females, y.plot_females, lty=1, col="red")

lci = vector(mode="numeric", length=length(x.plot_females))
uci = vector(mode="numeric", length=length(x.plot_females))
#lpi = vector(mode="numeric", length=length(x.plot_males))
#upi = vector(mode="numeric", length=length(x.plot_males))
for (i in 1:length(x.plot_females))
{
  x.star = matrix(c(1,x.plot_females[i],1,1*x.plot_females[i]), nrow=1, ncol=4)
  lci[i] = as.numeric(x.star%*%b1 - btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
  uci[i] = as.numeric(x.star%*%b1 + btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
  #lpi[i] = as.numeric(x.star%*%b1 - btcrit*sqrt(S1)*sqrt(1+x.star%*%B1%*%t(x.star)))
  #upi[i] = as.numeric(x.star%*%b1 + btcrit*sqrt(S1)*sqrt(1+x.star%*%B1%*%t(x.star)))
}

lines(x.plot_females, lci, lty=2, col="red")
lines(x.plot_females, uci, lty=2, col="red")
#lines(x.plot, lpi, lty=2, col="red")
#lines(x.plot, upi, lty=2, col="red")

legend("topleft", legend=c("Male", "Female", "Male Regression", "Male 95% HPD/PI", "Female Regression", "Female 95% HPD/PI"), lty=c(NA,NA,1,2,1,2), pch=c(1,1,NA,NA,NA,NA), col=c("black","red","black","black","red","red"), pt.cex=1, cex=.5)

### question-4
# male
x.star = matrix(c(1,75,0,0*75), nrow=1, ncol=4)
#Lower and Upper predictive intervals
lci_male = as.numeric(x.star%*%b1 - btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
uci_male = as.numeric(x.star%*%b1 + btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
# output of the 95% HPD for male 
lci_male
uci_male
# female
x.star = matrix(c(1,75,1,1*75), nrow=1, ncol=4)
#Lower and Upper predictive intervals
lci_female = as.numeric(x.star%*%b1 - btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
uci_female = as.numeric(x.star%*%b1 + btcrit*sqrt(S1)*sqrt(x.star%*%B1%*%t(x.star)))
# output of the 95% HPD for male 
lci_female
uci_female

