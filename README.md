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

<br>

![Alt text](images/weightsPonly.pdf)


**References**

Joseph Antonelli, Giovanni Parmigiani, and Francesca Dominici. **High-dimensional confounding adjustment using continuous spike and slab priors**. 2017. arXiv:1704.07532

