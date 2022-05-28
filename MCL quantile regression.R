library(dplyr)
library(reshape)
library(changepoint)
library(TTR)
library(ggplot2)
library(tidyverse)

# quantile check function
check <- function(x, tau =.5){
  x*(tau - (x<0))
}

# LASSO penalty
lasso <- function(x, lambda = 1){
  lambda*abs(x)
}

# LASSO derivation
lasso_deriv <- function(x, lambda = 1){
  u <- ifelse(x != 0,
              lambda*sign(x),
              0)
  u
}

# SCAD penalty
scad <- function(x, lambda = 1, a = 3.7){
  absx <- abs(x)
  ifelse(absx < lambda,
         lambda*absx,
         ifelse(absx < a*lambda,
                ( (a^2-1)*lambda^2 - (absx - a*lambda)^2) / (2*(a-1) ),
                (a+1)*lambda^2 / 2
                )
         )
}

# SCAD derivation
scad_deriv <- function(x, lambda = 1, a = 3.7){
  absx <- u <- abs(x)
  u[] <- 0
  index <- absx < a*lambda & absx > 0
  u[index] <- ifelse(absx[index] <= lambda,
                     lambda,
                     (a*lambda - absx[index])/ (a-1)
                     )
  u[index] <- u[index]* sign(x[index])
  u[x==0] <- lambda # because we take derivative as x approached 0 from above
  u
}

# MCP penalty
mcp <- function(x, lambda =1, a =3){
  absx <- abs(x)
  ifelse( absx < a*lambda,
          lambda*(absx - absx^2/(2*a*lambda)),
          a*lambda^2/2
          )
}

# MCP derivation
mcp_deriv <- function(x, lambda = 1, a = 3){
  u <- x
  u[] <- 0
  index <- abs(x) < a*lambda
  u[index] <- ifelse(x[index] == 0,
                     lambda,
                     lambda*sign(x[index]) - x[index]/a)
  u
}

# MCL penalty
mcl <- function(x, lambda =1, gamma = 0.5, a = 3){
  absx <- abs(x)
  ifelse( absx < a*(lambda - gamma),
          lambda*absx - absx^2/(2*a),
          gamma*absx + a*(lambda^2 - gamma^2)/2
          )
}

# MCL derivation
mcl_deriv <- function(x, lambda =1, gamma = .5, a = 3){
  u <- x
  u[] <- lambda
  index <- abs(x) > 0
  u[index] <- ifelse(abs(x[index]) < a*(lambda - gamma),
                     lambda*sign(x[index]) - x[index]/a,
                     gamma*sign(x[index])
                     )
  u
}

# break 설정
x <- seq(-4, 4, 0.01)

# penalties 그래프 그리기
par(mfrow = c(2,4))
plot(x, mcl(x), type="l", ylim = c(0, 4), main = "MCL")
plot(x, mcl_deriv(x), type="l", ylim = c(-1, 1), main = "MCL_derivation")
plot(x, mcp(x), type="l", ylim = c(0, 4), main = "MCP")
plot(x, mcp_deriv(x), type="l", ylim = c(-1, 1), main = "MCP_derivation")
plot(x, scad(x), type="l", ylim = c(0, 4), main = "SCAD")
plot(x, scad_deriv(x), type="l", ylim = c(-1, 1), main = "SCAD_derivation")
plot(x, lasso(x), type="l", ylim = c(0, 4), main = "LASSO")
plot(x, lasso_deriv(x), type="l", ylim = c(-1, 1), main = "LASSO_derivation")


# QICD
# Quantile Regression Model by using QICD Algorithm

help(pendrive)
