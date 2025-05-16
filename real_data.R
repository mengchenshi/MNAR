library(tidyverse)
library(haven)
library(BB)
library(numDeriv)
job_data |> summary()

job_data |> filter(otherwelf_miss == 1) |> nrow()
job_data |> filter(is.na(earny4)) |> nrow()
job_data |> filter(otherwelf_miss == 1, is.na(earny4)) |> nrow()




real_data <-  job_data |> select(earny4, treatment, otherwelf, otherwelf_miss)
real_data <- real_data %>%
  mutate(across(where(is.labelled), as.numeric))

real_data <-  real_data |> mutate(
  Y_obs = if_else(otherwelf_miss == 1, NA, log(earny4 + 1)),
  X_obs = if_else(is.na(Y_obs), NA, otherwelf),
  Z = treatment,
  M = if_else(!is.na(X_obs), 1, 0)
) 


set.seed(1)
real_data <- real_data |> slice_sample(prop = 0.6) 


ggplot(real_data, aes(x = Y_obs)) + geom_density()



# proposal distribution for X, density
hx <- function(x, para){
  #px <- 1 / (1 + exp(- para[['prop_x']]))
  px <- para[['prop_x']]
  
  prop = px ^ x * (1 - px) ^ (1 - x)
  
  prop
}


# proposal distribution for Y, density
hy <- function(x, y, z, para){
  
  (1 / sqrt(2 * pi * para[['sigma_x1']] ^ 2) * exp(- (y - (para[['theta_10']] + para[['theta_1z']] * z)) ^ 2 / (2 * para[['sigma_x1']] ^ 2))) ^ x *
    (1 / sqrt(2 * pi * para[['sigma_x0']] ^ 2) * exp(- (y - (para[['theta_00']] + para[['theta_0z']] * z)) ^ 2 / (2 * para[['sigma_x0']] ^ 2))) ^ (1 - x) 
  
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
  
  
  Y_imp <- map2(X_imp, Z, \(x, z) I(x == 1) * rnorm(2 * B, mean = para[['theta_10']] + para[['theta_1z']] * z, sd = para[['sigma_x1']]) +
                  I(x == 0) * rnorm(2 * B, mean = para[['theta_00']] + para[['theta_0z']] * z, sd = para[['sigma_x0']])) #
  
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
  
  ((px * (1 / sqrt(2 * pi * para[['sigma_x1']] ^ 2) * exp(- (y - (para[['theta_10']] + para[['theta_1z']] * z)) ^ 2 / (2 * para[['sigma_x1']] ^ 2))) * 
      pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_10']] + para[['theta_1z']] * z))) / para[['sigma_x1']])) ^ x *
      ((1 - px) * (1 / sqrt(2 * pi * para[['sigma_x0']] ^ 2) * exp(- (y - (para[['theta_00']] + para[['theta_0z']] * z)) ^ 2 / (2 * para[['sigma_x0']] ^ 2))) * 
         pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_00']] + para[['theta_0z']] * z))) / para[['sigma_x0']])) ^ (1 - x)) ^ m *
    ((px * (1 / sqrt(2 * pi * para[['sigma_x1']] ^ 2) * exp(- (y - (para[['theta_10']] + para[['theta_1z']] * z)) ^ 2 / (2 * para[['sigma_x1']] ^ 2))) * 
        pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_10']] + para[['theta_1z']] * z))) / para[['sigma_x1']], lower.tail = FALSE)) ^ x *
       ((1 - px) * (1 / sqrt(2 * pi * para[['sigma_x0']] ^ 2) * exp(- (y - (para[['theta_00']] + para[['theta_0z']] * z)) ^ 2 / (2 * para[['sigma_x0']] ^ 2))) * 
          pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_00']] + para[['theta_0z']] * z))) / para[['sigma_x0']], lower.tail = FALSE)) ^ (1 - x)) ^ (1 - m) *
    pz ^ z * (1 - pz) ^ (1 - z)
  
  
  
  # 
  # px ^ x * (1 - px) ^ (1 - x) * 
  #   pz ^ z * (1 - pz) ^ (1 - z) *
  #   (1 / sqrt(2 * pi * para[['sigma_x1']] ^ 2) * exp(- (y - (para[['theta_10']] + para[['theta_1z']] * z)) ^ 2 / (2 * para[['sigma_x1']] ^ 2))) ^ x *
  #   (1 / sqrt(2 * pi * para[['sigma_x0']] ^ 2) * exp(- (y - (para[['theta_00']] + para[['theta_0z']] * z)) ^ 2 / (2 * para[['sigma_x0']] ^ 2))) ^ (1 - x) *
  #   (pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_10']] + para[['theta_1z']] * z))) / para[['sigma_x1']]) ^ x *
  #      pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_00']] + para[['theta_0z']] * z))) / para[['sigma_x0']]) ^ (1 - x) ) ^ m *
  #   (pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_10']] + para[['theta_1z']] * z))) / para[['sigma_x1']], lower.tail = FALSE) ^ x *
  #      pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_00']] + para[['theta_0z']] * z))) / para[['sigma_x0']], lower.tail = FALSE) ^ (1 - x) ) ^ (1 - m)
  
}

# 得分函数




score_fun <- function(para, data){
  x <- data$X; y <- data$Y; m <- data$M; z = data$Z
  
  
  result <- m * (x * (dbinom(x, 1, para[['prop_x']], log = TRUE) + (dnorm(y, mean = para[['theta_10']] + para[['theta_1z']] * z, sd = para[['sigma_x1']], log = TRUE)) + (pnorm(((para[['alpha_0']] + z * para[['alpha_z']])  + para[['beta']] * (y - (para[['theta_10']] + para[['theta_1z']] * z))) / para[['sigma_x1']], log.p = TRUE))) + 
                   (1 - x) * (dbinom(x, 1, para[['prop_x']], log = TRUE) + (dnorm(y, mean = para[['theta_00']] + para[['theta_0z']] * z, sd = para[['sigma_x0']], log = TRUE)) + (pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_00']] + para[['theta_0z']] * z))) / para[['sigma_x0']], log.p = TRUE)))) +
    (1 - m) * (x * (dbinom(x, 1, para[['prop_x']], log = TRUE) + (dnorm(y, mean = para[['theta_10']] + para[['theta_1z']] * z, sd = para[['sigma_x1']], log = TRUE)) + (pnorm(((para[['alpha_0']] + z * para[['alpha_z']])  + para[['beta']] * (y- (para[['theta_10']] + para[['theta_1z']] * z))) / para[['sigma_x1']], log.p = TRUE, lower.tail = FALSE))) + 
                 (1 - x) * (dbinom(x, 1, para[['prop_x']], log = TRUE) + (dnorm(y, mean = para[['theta_00']] + para[['theta_0z']] * z, sd = para[['sigma_x0']], log = TRUE)) + (pnorm(((para[['alpha_0']] + z * para[['alpha_z']]) + para[['beta']] * (y - (para[['theta_00']] + para[['theta_0z']] * z))) / para[['sigma_x0']], log.p = TRUE, lower.tail = FALSE)))) + 
    z * log(para[['prop_z']]) + (1 - z) * log(1 - para[['prop_z']])
  
  
  result
}


# 权重函数更新

# bind_cols(XX_imp, Y_imp, .name_repair = 'minimal')|> select(-1) |> select(rowname, everything())

Data_final <- function(data, para, B){
  data_imp <- Imputation_fun(data, para, B) |> mutate(w = NA, M = 0) |> rename(X = X_imp, Y = Y_imp)
  
  data <- data |> rownames_to_column()
  data_obs <- data  |> filter(M == 1) |> mutate(w = NA) |> rename(X = X_obs, Y = Y_obs)
  
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




# 填补之后的最终数据，在初始值para下
B <- 50
lo <- rep(c(0.01, -Inf, 0.01, -Inf), c(2, 4, 2, 3))
hi <- rep(c(0.99, Inf, Inf, Inf), c(2, 4, 2, 3))

n <- 500 # 500次
Estimators_continuous_y <- vector(mode = 'list', length = n)


for (i in 1 : n) {
  set.seed(i)
  Index <- sample(1:nrow(real_data), replace = TRUE)
  data <- real_data[Index, ]
  # 原始数据下的极大似然估计
  px <- mean(data$X_obs, na.rm = TRUE)
  pz <- mean(data$Z, na.rm = TRUE)
  model_x1 <- lm(Y_obs ~ Z, data = filter(data, X_obs == 1))
  model_x0 <- lm(Y_obs ~ Z, data = filter(data, X_obs == 0))
  
  
  para <-list(prop_x = px, prop_z = pz, theta_10 = model_x1$coefficients[[1]], theta_1z = model_x1$coefficients[[2]], 
              theta_00 = model_x0$coefficients[[1]], theta_0z = model_x0$coefficients[[2]], sigma_x1 = sd(model_x1$residuals), sigma_x0 = sd(model_x0$residuals), 
              alpha_0 = runif(1, -2, 2), alpha_z = runif(1, -2, 2), beta = runif(1, -2, 2)) |> unlist()

  data_final <- Data_final(data, para, B)
  t <- 0
  temp_para <- rep(0, 11)
  h <- hy(x = data_final$X, y = data_final$Y, z = data_final$Z, para) 
  
  while(sum((para - temp_para)^2) >= 1e-04 & t <= 1000){
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
    para <- BBoptim(par = para, fn = score_equation, data_final = data_final, weight = weight, lower = lo, upper = hi, control = list(trace = FALSE))$par
    cat(sprintf("第%d次优化后的参数为:\n", t))
    print(para)
  }
  Estimators_continuous_y[[i]] <- para
  cat(sprintf("\033[1m\033[31m第 %d 次循环完成\033[0m\n", i))
}
  
 

Estimators_continuous_y_mle <- vector(mode = 'list', length = n)

for (i in 1 : n) {
  set.seed(i)
  Index <- sample(1:nrow(real_data), replace = TRUE)
  data <- real_data[Index, ]
  # 原始数据下的极大似然估计
  px <- mean(data$X_obs, na.rm = TRUE)
  pz <- mean(data$Z, na.rm = TRUE)
  
  model_x1 <- lm(Y_obs ~ Z , data = filter(data, X_obs == 1))
  model_x0 <- lm(Y_obs ~ Z , data = filter(data, X_obs == 0))
  
  
  para <-list(prop_x = px, prop_z = pz, theta_10 = model_x1$coefficients[[1]], theta_1z = model_x1$coefficients[[2]],
              theta_00 = model_x0$coefficients[[1]],theta_0z = model_x0$coefficients[[2]], sigma_x1 = sd(model_x1$residuals), sigma_x0 = sd(model_x0$residuals)
  ) |> unlist()
  
  Estimators_continuous_y_mle[[i]] <- para
  cat(sprintf("\033[1m\033[31m第 %d 次循环完成\033[0m\n", i))
}


















































