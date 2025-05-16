# 二值Y 
# Mechanism 2 : R depend on T, X

library(tidyverse)
library(BB)
library(numDeriv)

# 协变量X离散二值时
Generate_data <- function(N, para){
  px <- para[['prop_x']]
  pz <- para[['prop_z']]
  
  X_com <- rbinom(N, size = 1, prob = px)
  
  Z <- rbinom(N, size = 1, prob = pz)
  
  # P(Y = 1 | x, z)
  Y_com <- rbinom(N, size = 1, prob = pnorm(para[['theta_0']] + para[['theta_z']] * Z + para[["theta_x"]] * X_com))
  
  M <- rbinom(N, size = 1, prob = pnorm(para[['alpha_0']] + para[['alpha_z']] * Z + para[['alpha_x']] * X_com))
  
  X_obs <- if_else(M == 1, X_com, NA)
  Y_obs <- if_else(M == 1, Y_com, NA)
  
  tibble(X_com, Y_com, Z, M, X_obs, Y_obs)
}



# proposal distribution for X, density
hx <- function(x, para){
  #px <- 1 / (1 + exp(- para[['prop_x']]))
  
  px <- para[['prop_x']]
  
  
  prop = px ^ x * (1 -px) ^ (1 - x)
  
  prop
}


# proposal distribution for Y, density
hy <- function(x, y, z, para){
  
  pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x) ^ y *
    pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x, lower.tail = FALSE) ^ (1 - y)
  
}


# 填补X
hx_imp <- function(data, para, B, index){
  
  X_imp <- vector(mode = 'list', length(index))
  names(X_imp) <- index
  
  X_imp <- map(X_imp, \(x) rep(c(1, 0), c(B, B))) #
  
  X_imp
}

# 填补Y
hy_imp <- function(data, para, B, index){
  Y_imp <- vector(mode = 'list', length(index))
  names(Y_imp) <- index
  
  Z <- vector(mode = 'list', length(index))
  Z <- map(index, \(x) data$Z[as.numeric(x)] |> rep(2 * B))
  names(Z) <- index
  
  X_imp <- hx_imp(data, para, B, index)
  
  
  Y_imp <- map2(X_imp, Z, \(x, z) rbinom(2 * B, size = 1, prob = pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x))) #
  
  Y_imp
}

Imputation_fun <- function(data, para, B){
  index <- data |> rownames_to_column() |> filter(is.na(X_obs)) |> select(rowname) |> pull()
  Z <- vector(mode = 'list', length(index))
  Z <- map(index, \(x) data$Z[as.numeric(x)] |> rep(2 * B))
  names(Z) <- index
  Z <- Z |> as_tibble() |> pivot_longer(cols = everything(), names_to = 'rowname', values_to = 'Z') |> arrange(as.numeric(rowname))
  
  X_imp <- hx_imp(data, para, B, index)
  Y_imp <- hy_imp(data, para, B, index)
  
  X_imp <- X_imp |> as_tibble() |> pivot_longer(cols = everything(), names_to = 'rowname', values_to = 'X_imp') |> arrange(as.numeric(rowname))
  Y_imp <- Y_imp |> as_tibble() |> pivot_longer(cols = everything(), names_to = 'rowname', values_to = 'Y_imp') |> arrange(as.numeric(rowname))
  
  
  
  data_imp <- bind_cols(Z, X_imp, Y_imp, .name_repair = 'minimal') |> select(1, Z, X_imp, Y_imp) |>  select(rowname, everything())
  
}

# 联合密度函数f(x, y, m)
joint_density <- function(para, data){
  x <- data$X; y <- data$Y; m <- data$M; z = data$Z
  
  px <- para[['prop_x']]
  pz <- para[['prop_z']]
  
  px ^ x * (1 - px) ^ (1 - x) * 
    pz ^ z * (1 - pz) ^ (1 - z) *
    pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x) ^ y * 
    (pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x, lower.tail = FALSE)) ^ (1 - y) *
    (pnorm(para[['alpha_0']] + para[['alpha_z']] * z + para[['alpha_x']] * x)) ^ m *
    (pnorm(para[['alpha_0']] + para[['alpha_z']] * z + para[['alpha_x']] * x, lower.tail = FALSE)) ^ (1 - m)
  
}

# 得分函数
score_fun <- function(para, data){
  x <- data$X; y <- data$Y; m <- data$M; z = data$Z
  px <- para[['prop_x']]
  pz <- para[['prop_z']]
  
  result <- z * log(pz) + (1 - z) * log(1 - pz) +
    x * log(px) + (1 - x) * log(1 - px) +
    y * (pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x, log.p = TRUE)) +
    (1 - y) * (pnorm(para[['theta_0']] + para[['theta_z']] * z + para[["theta_x"]] * x, lower.tail = FALSE, log.p = TRUE)) +
    m * (pnorm(para[['alpha_0']] + para[['alpha_z']] * z + para[['alpha_x']] * x, log.p = TRUE)) +
    (1 - m) * (pnorm(para[['alpha_0']] + para[['alpha_z']] * z + para[['alpha_x']] * x, lower.tail = FALSE, log.p = TRUE))
  
  result
}


# 权重函数更新

# bind_cols(XX_imp, Y_imp, .name_repair = 'minimal')|> select(-1) |> select(rowname, everything())

Data_final <- function(data, para, B){
  data_imp <- Imputation_fun(data, para, B) |> mutate(w = NA, M = 0) |> rename(X = X_imp, Y = Y_imp)
  
  data <- data |> rownames_to_column()
  data_obs <- data |> select(-X_com, -Y_com) |> filter(M == 1) |> mutate(w = NA) |> rename(X = X_obs, Y = Y_obs)
  
  data_final <- bind_rows(data_obs, data_imp) 
  data_final
}


weight_update <- function(para, data_final, h){
  f <- joint_density(para, data_final)
  # h <- hx(para) * hy(x = data_final$X, y = data_final$Y, para)
  
  weight <- data_final |> mutate(w = f / h) |> group_by(rowname) |> mutate(w = w / sum(w)) |> pull(w)
  
  weight
}

score_equation <- function(para, data_final, weight){
  
  #score <- jacobian(score_fun, x = para, data = data_final)
  score <- score_fun(para = para, data = data_final)
  ll <- (weight * score) |> sum()
  
  -ll
  
}


N <- 1000
parameters <- list(prop_x = 0.5, prop_z = 0.5, theta_0 = -1 ,  theta_z = 2 , theta_x = 1,  alpha_0 = 1,  alpha_z = -1, alpha_x = 1.5)


# 填补之后的最终数据，在初始值para下
B <- 50

lo <- rep(c(0.01, -Inf), c(2, 6))
hi <- rep(c(0.99, Inf), c(2, 6))

#h <- 1

n <- 500 # 500次
Estimators_binary_y <- vector(mode = 'list', length = n)

for (i in 1 : n) {
  set.seed(i)
  
  data <- Generate_data(N, para = parameters)
  
  # 原始数据下的极大似然估计
  px <- mean(data$X_obs, na.rm = TRUE)
  pz <- mean(data$Z, na.rm = TRUE)
  model <- glm(Y_obs ~ Z + X_obs, data = data, family = binomial(link = "probit"))
  
  
  
  # 填补之后的最终数据，在初始值para下
  para <- list(prop_x = px, prop_z = pz, theta_0 = model$coefficients[[1]],  theta_z = model$coefficients[[2]], theta_x = model$coefficients[[3]], 
               alpha_0 = runif(1, -5, 5), alpha_z = runif(1, -5, 5), alpha_x = runif(1, -5, 5)) |> unlist()
  lo <- rep(c(0.01, -Inf), c(2, 6))
  hi <- rep(c(0.99, Inf), c(2, 6))
  
  data_final <- Data_final(data, para, B)
  t <- 0
  temp_para <- rep(0, 8)
  h <- hy(x = data_final$X, y = data_final$Y, z = data_final$Z, para) 
  
  while(sum((para[2:5] - temp_para[2:5])^2) >= 1e-06 & t <= 1000){
    num = 0
    t = t + 1
    temp_para <- para
    weight <- weight_update(para, data_final, h = h)
    
    while(sum(any(weight[weight != 1] > 5 / B)) > 1) {
      if(t == 1 | num > 0) {
        B = B + 10
        print(B)
        num = 0
      }
      data_final <- Data_final(data, para, B)
      h <- hy(x = data_final$X, y = data_final$Y, z = data_final$Z, para) 
      
      weight <- weight_update(para, data_final, h = h)
      num = num + 1
    }
    
    #para <- optim(par = para, fn = score_equation, data_final = data_final, weight = weight, lower = lo, upper = hi,  method = "L-BFGS-B")$par
    para <- BBoptim(par = para, fn = score_equation, data_final = data_final, 
                    weight = weight, lower = lo, upper = hi, control = list(trace = FALSE))$par
    cat(sprintf("第%d次优化后的参数为:\n", t))
    print(para)
  } 
  Estimators_binary_y[[i]] <- para
  cat(sprintf("\033[1m\033[31m第 %d 次循环完成\033[0m\n", i))
}


Estimators_binary_y_mle <- vector(mode = 'list', length = n)
for (i in 1 : n) {
  set.seed(i)
  
  data <- Generate_data(N, para = parameters)
  
  # 原始数据下的极大似然估计
  px <- mean(data$X_obs, na.rm = TRUE)
  pz <- mean(data$Z, na.rm = TRUE)
  model <- glm(Y_obs ~ Z + X_obs, data = data, family = binomial(link = "probit"))
  
  # 填补之后的最终数据，在初始值para下
  para <- list(prop_x = px, prop_z = pz, theta_0 = model$coefficients[[1]],  theta_z = model$coefficients[[2]], theta_x = model$coefficients[[3]]) |> unlist()
  
  Estimators_binary_y_mle[[i]] <- para
  cat(sprintf("\033[1m\033[31m第 %d 次循环完成\033[0m\n", i))
}




