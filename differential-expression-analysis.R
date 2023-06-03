#conducting DEA
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)

res

# or to shrink log fold changes association with condition:
resLFC <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
resLFC

#p-values and adjusted p-values
library("BiocParallel")
register(MulticoreParam(4))

resOrdered <- res[order(res$pvalue),]
summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)

#Likelihood-Ratio Test

#function to calculate the likelihood of observing a particular 
#sequence of heads and tails (d) given a probability of flipping 
#heads (p). 

likelihood <- function(d, p){
  l <- 1
  for (i in d){
    if (i == 1){
      l <- l*p
    } else {
      l <- l *(1-p)
    }
  }
  return(l)
}

# creates a vector of potential values for theta, the probability of flipping a heads
potential_thetas <- seq(0, 1, by=.05)

barplot(likelihood(s_quarter, potential_thetas), potential_thetas, names= potential_thetas, width=.9, col = "steelblue", 
        main="Likelihood of Observing Data at Different values of Theta",
  xlab="θ (prob of flipping a heads)", ylab="Likelihood of Data Given θ")

set.seed(7)
s_penny <- sample(c(0,1), 10, replace = TRUE, prob = c(.8,.2))
quarter_then_penny = c(append(s_quarter,s_penny))
cat("quarter flips:", quarter_then_penny[1:10], "\n \n")
cat("penny flips:", quarter_then_penny[11:20], "\n \n")
cat("all flips:", quarter_then_penny)

barplot(likelihood(quarter_then_penny, potential_thetas), potential_thetas, names= potential_thetas, 
        width=.9, col = "steelblue", 
        main="Likelihood of Observing Data at Different Values of θ",
  xlab="θ (probability of flipping a heads)", ylab="Likelihood of Data Given θ")

two_coin_matrix <- (likelihood(s_quarter,potential_thetas)) %o% (likelihood(s_penny, potential_thetas))
# Finds coordinates of maximum value in matrix
xy <- as.vector(which(two_coin_matrix==max(two_coin_matrix), arr.ind=T))


install.packages("plot3D")
library(plot3D)
hist3D(x=(potential_thetas), y=potential_thetas, z=two_coin_matrix, ticktype="detailed", space=0.35, phi = 35, color = rainbow,
      xlab="quarter θ", ylab = "penny θ", zlab = "likelihood of data given q_θ and p_θ", main="Likelihood of Observing Data Given Different Values of q_θ and p_θ")


#Wald Test

library(aod)
wald.test(Sigma = vcov(model), b = coef(model), Terms = 1:2)
