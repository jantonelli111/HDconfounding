# HDconfounding

This respository is an R package to implement the methods in "High-dimensional confounding adjustment using continuous spike and slab priors" by Joseph Antonelli, Giovanni Parmigiani, and Francesca Dominici.
The paper can be found at the following link:

https://arxiv.org/pdf/1704.07532.pdf

This software can be loaded by entering the following into R

```{r, echo=TRUE, message=FALSE}
library(devtools)
install_github(repo = "jantonelli111/HDconfounding")
library(HDconfounding)
```
A brief example of how to run the code is below. For a more detailed description of using the software, use the vignette in the vignettes folder of this repository.

```{r, eval=FALSE}
n = 200
p = 20
x = matrix(rnorm(n*p), n, p)
z = rbinom(n, 1, p=pnorm(0.7*x[,1] + 0.3*x[,2]))
y = rnorm(n, mean=z + 0.3*x[,1] + 0.6*x[,2] + 0.5*x[,3], sd=1)

sslEB = SSL(y=y, z=z, x=x, nScans=3000, burn=1000, thin=2, lambda0="EB")

sslEB$TreatEffect
sslEB$TreatEffectCI
```

If you wish to see the vignette associated with the package that helps to illustrate its usage through examples, use the following lines of code. This will take a few minutes to build.

```{r, echo=TRUE, message=FALSE}
library(devtools)
install_github(repo = "jantonelli111/HDconfounding", build_vignettes = TRUE)
library(HDconfounding)
```
Then, to view the vignette simply type into R

```{r, echo=TRUE, message=FALSE}
vignette("HDconfounding", package="HDconfounding")
```
The PDF of the vignette can also be found in the vignettes folder of this github repository, if the user does not want to build it in their installation.

Currently, the software only allows for fully Bayesian inference with MCMC, but we are in the process of writing software that also estimates posterior modes of treatment effects. This will not provide valid measures of uncertainty, but will estimate quantities with much less computation time and will allow for automatic confounder selection as some confounders will have coefficients hard thresholded to zero. If you have any questions regarding the software don't hesitate to contact Joseph Antonelli at jantonelli111@gmail.com, and please report any bugs that you find!

**References**

Joseph Antonelli, Giovanni Parmigiani, and Francesca Dominici. **High-dimensional confounding adjustment using continuous spike and slab priors**. 2017. arXiv:1704.07532

