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
Currently, the software allows for continuous outcomes only, but we are in the process of adapting it to allow for binary outcomes as well. If you have any questions regarding the software don't hesitate to contact Joseph Antonelli at jantonelli111@gmail.com, and please report any bugs that you find!

