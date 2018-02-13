# HDconfounding

This respository is an R package to implement the methods in "High-dimensional confounding adjustment using continuous spike and slab priors" by Joseph Antonelli, Giovanni Parmigiani, and Francesca Dominici.
The paper can be found at the following link:

https://arxiv.org/pdf/1704.07532.pdf

Currently, the software only allows for fully Bayesian inference with MCMC, but we are in the process of writing software that also estimates posterior modes of treatment effects. This will not provide valid measures of uncertainty, but will estimate quantities with much less computation time and will allow for automatic confounder selection as some confounders will have coefficients hard thresholded to zero. If you have any questions regarding the software don't hesitate to contact Joseph Antonelli at jantonelli111@gmail.com, and please report any bugs that you find!

# How to use the HDconfounding

In this section we will show you how to use the main functions of the HDconfounding R package. We will simulate data under both a binary and continuous outcome and show how to analyze them using the various functions in the package. Throughout we will use a sample size of $$n=200$$ and p=200. In general the software will work for larger sample sizes and covariate dimensions, however, we the code might take substantially longer for covariate dimensions in the thousands. 

### Loading and building the package

To build the package the user must install the devtools R package and then run the following code.

```
library(devtools)
install_github(repo = "jantonelli111/HDconfounding")
library(HDconfounding)
```

# Binary treatments

### Continuous outcome, homogeneous treatment effect

First let's generate some data to practice estimating treatment effects with:

```{r, eval=TRUE}
n = 200
p = 200
x = matrix(rnorm(n*p), n, p)
z = rbinom(n, 1, p=pnorm(0.7*x[,1] + 0.3*x[,2]))
y = rnorm(n, mean=z + 0.3*x[,1] + 0.6*x[,2] + 0.5*x[,3], sd=1)
```

So our sample size is 200, we have 200 covariates, the first 2 covariates are important for confounding adjustment, and the 3rd covariate is a predictor of the outcome. In this case, the treatment effect is homogeneous, i.e does not vary for different values of the covariates. We will now estimate it using both the homogeneous and heterogeneous versions of our approach. Both should yield unbiased answers in this setting, though we expect the heterogeneous approach to be slightly less efficient as it is a more complex model with more parameters.

First, to estimate the homogeneous treatment effect we can run the following

```{r, eval=TRUE}
sslEB = SSL(y=y, z=z, x=x, nScans=3000, burn=1000, thin=2, lambda0="EB")
```

This model estimates the tuning parameter $\lambda_0$ with empirical Bayes. We recommend this as it is difficult to have prior knowledge about what value this parameter should take. If the user wants to specify their own value of $\lambda_0$, then weight must also be provided (w in the main manuscript) and we recommend a value of 0.05.

```{r, eval=TRUE}
ssl = SSL(y=y, z=z, x=x, nScans=3000, burn=1000, thin=2, lambda0=15, weight = 0.05)
```

The empirical Bayes estimate is calculated using MCMC so the running time is a little longer than pre-selecting $\lambda_0$, which is one reason one might want to choose their own value if they have a large data set. It is important to note though that a poorly chosen value of $\lambda_0$ can lead to poor estimation so we highly recommend the EB option. 


### Continuous outcome, heterogeneous treatment effect

Now let's analyze the same data, but allowing for heterogeneous treatment effects. 

```{r, eval=TRUE}
sslEBhetero = SSLhetero(y=y, z=z, x=x, nScans=3000, burn=1000, thin=2, lambda0="EB")
```

We can plot the estimates and intervals to compare the three models we have run. The red line is the truth.
```{r, eval=TRUE, fig.height=6, fig.width=6}
estimates = c(ssl$TreatEffect, sslEB$TreatEffect, 
              sslEBhetero$TreatEffect)
CIlower = c(ssl$TreatEffectCI[1], sslEB$TreatEffectCI[1], 
            sslEBhetero$TreatEffectCI[1])
CIupper = c(ssl$TreatEffectCI[2], sslEB$TreatEffectCI[2],
            sslEBhetero$TreatEffectCI[2])

plot(1:3, estimates, pch=17, cex=2, ylim = range(c(CIlower, CIupper)) + c(-0.05, 0.05),
     xlab="", ylab="Estimates", axes=FALSE)
segments(x0 = 1:3, x1=1:3, y0=CIlower, y1=CIupper, lwd=3)
abline(h = 1, lty=2, col="red")
axis(2)
axis(1, at=1:3, c("homo", "EB-homo", "EB-hetero"))
```

![Alt text](images/Plot1.png)


**References**

Joseph Antonelli, Giovanni Parmigiani, and Francesca Dominici. **High-dimensional confounding adjustment using continuous spike and slab priors**. 2017. arXiv:1704.07532

