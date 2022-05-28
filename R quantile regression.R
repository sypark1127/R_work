install.packages("quantreg")
library(quantreg)

install.packages("rqPen")
library(rqPen)

install.packages("mlbench")
library(mlbench)


data("BostonHousing")
head(BostonHousing, 5)

lm(medv~., data = BostonHousing)

X = data.matrix(BostonHousing[, -14])
y = BostonHousing[, 14]

head(y)
head(X)

# rq.lasso.fit
# quantile regression with LASSO penalty
# rq.lasso.fit(x, y, tau =.5, lambda = NULL, weights = NULL,
#              intercept = TRUE, coef.cutoff = .00000001,
#              method = "br" or "fn", penVars = NULL, scalex =TRUE)
rq.lasso.fit(X, y,tau = 0.75, lambda = 0.05, intercept = TRUE,

# QICD
# penalized quantile regression with QICD
# QICD
# QICD produces penalized quantile regression estimates using the QICD algorithm
# this function can handle the LASSO, SCAD, and MCP penalities

# QICD(y, x, tau = .5, lambda, intercept = TRUE, penalty = "SCAD", 
#      initial_beta = NULL,  maxin = 100, maxout = 20, eps = 1e-05,
#      coef.cutoff=1e-08, a = 3.7, scalex = TRUE)

QICD(y, X, tau = 0.75, intercept = TRUE, penalty = "LASSO", lambda = 0.05,
     scalex = TRUE)
