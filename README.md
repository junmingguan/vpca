# vpca
Implementation of the variational PCA algorithm (C. M. Bishop, "Variational principal components," *1999 Ninth International Conference on Artificial Neural Networks ICANN 99* (Conf. Publ. No. 470), Edinburgh, UK, 1999, pp. 509-514 vol.1, doi: [10.1049/cp:19991160](https://doi.org/10.1049/cp:19991160)).

To use the package, do the following steps:

-   Run the command `> git clone https://github.com/junmingguan/vpca.git`;
-   `> cd PATH/TO/vpca`;
-   In R, install the load the package:
```r
    install.packages("PATH/TO/vpca", repos = NULL, type = "source")
    library(vpca)
```

Alternatively,

```r
install.packages("remotes")
remotes::install_github("junmingguan/vpca")
library(vpca)
```

A simple example:
```r
d <- 10
n <- 100
Y <- matrix(data = 0, nrow = d, ncol = n)
alpha <- rep(1, d)
alpha[1:5] <- q:5
for (i in 1:d) {
    Y[i, ] <- rnorm(n, 0, alpha[i])
}
res_vpca <- init_svd(Y, q = 6, include_mu = F, alpha_for_X = F, alpha_MLE = F)
res_vpca <- vpcar(res_vpca, max.iter = 2000, min.iter = 1, epsilon=0.0001)
names(res_vpca)
# contains first and second moments, and ELBO at each iteration
# 'E_W' 'E_X' 'E_mu' 'Sigma_W' 'Sigma_X' 'E_alpha' 'E_tau' 'elbo_list'
```