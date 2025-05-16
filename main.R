library(tidyverse)
M1_binary <- Estimators_binary_y |> reduce(rbind) |>  as_tibble()
M1_binary <- M1_binary |> mutate(CRD = pnorm(theta_0 + theta_z + theta_x) * prop_x + pnorm(theta_0 + theta_z) * (1 - prop_x) -
                      pnorm(theta_0 + theta_x) * prop_x - pnorm(theta_0) * (1 - prop_x), 
                    CRD_1 = pnorm(theta_0 + theta_z + theta_x) -
                      pnorm(theta_0 + theta_x),
                    CRD_0 = pnorm(theta_0 + theta_z) -
                      pnorm(theta_0))
pnorm(-1 + 2 + 1) * 0.5 + pnorm(-1 + 2) * (1 - 0.5) -
  pnorm(-1 + 1) * 0.5 - pnorm(-1) * (1 - 0.5)

pnorm(-1 + 2 + 1) -
  pnorm(-1 + 1)

pnorm(-1 + 2) -
  pnorm(-1)



mean_bias <- colMeans(M1_binary) - c(0.5, 0.5, -1, 2, 1, 1.5, -1.5, 0.5799697, 0.4772499, 0.6826895)
sd <- apply(M1_binary, 2, sd) 

M1_binary_cc <- Estimators_binary_y_mle |> reduce(rbind) |>  as_tibble()
M1_binary_cc <- M1_binary_cc |> mutate(CRD = pnorm(theta_0 + theta_z + theta_x) * prop_x + pnorm(theta_0 + theta_z) * (1 - prop_x) -
                                   pnorm(theta_0 + theta_x) * prop_x - pnorm(theta_0) * (1 - prop_x), 
                                 CRD_1 = pnorm(theta_0 + theta_z + theta_x) -
                                   pnorm(theta_0 + theta_x),
                                 CRD_0 = pnorm(theta_0 + theta_z) -
                                   pnorm(theta_0))
mean_bias_cc <- colMeans(M1_binary_cc) - c(0.5, 0.5, -1, 2, 1, 0.5799697, 0.4772499, 0.6826895)
sd_cc <- apply(M1_binary_cc, 2, sd) 


summary_M1_binary <- data.frame(
  Bias_EM = mean_bias,
  Sd_EM = sd,
  Bias_cc = c(mean_bias_cc[1 : 5], NA, NA, mean_bias_cc[6 : 8]),
  Sd_cc = c(sd_cc[1 : 5], NA, NA, sd_cc[6 : 8])
) 
  

M1_continuous <- Estimators_continuous_y |> reduce(rbind) |>  as_tibble()
M1_continuous <- M1_continuous |> mutate(
  CRD = (exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) * prop_x + exp(theta_0 + theta_z + sigma_x ^ 2 / 2) * (1 - prop_x) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2) * prop_x -
    exp(theta_0 + sigma_x ^ 2 / 2) * (1 - prop_x)) ,
  CRD_1 =  exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2),
  CRD_0 =  exp(theta_0 + theta_z + sigma_x ^ 2 / 2) - exp(theta_0 + sigma_x ^ 2 / 2)
)

(exp(-1 + 2 - 1 + 1 ^ 2 / 2) + exp(-1 + 2 + 1 ^ 2 / 2)) / 2 - ( exp(-1 - 1 + 1 ^ 2 / 2) + exp(-1 + 1 ^ 2 / 2)) / 2
exp(-1 + 2 - 1 + 1 ^ 2 / 2) - exp(-1 - 1 + 1 ^ 2 / 2)
exp(-1 + 2 + 1 ^ 2 / 2) - exp(-1 + 1 ^ 2 / 2)

mean_bias <- colMeans(M1_continuous) - c(0.5, 0.5, -1, 2, -1, 1, 1.5, -2, 2.650375, 1.425591, 3.875158)
sd <- apply(M1_continuous, 2, sd) 

M1_continuous_cc <- Estimators_continuous_y_mle |> reduce(rbind) |>  as_tibble()
M1_continuous_cc <- M1_continuous_cc |> mutate(
  CRD = (exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) * prop_x + exp(theta_0 + theta_z + sigma_x ^ 2 / 2) * (1 - prop_x) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2) * prop_x -
           exp(theta_0 + sigma_x ^ 2 / 2) * (1 - prop_x)),
  CRD_1 =  exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2),
  CRD_0 =  exp(theta_0 + theta_z + sigma_x ^ 2 / 2) - exp(theta_0 + sigma_x ^ 2 / 2)
)
mean_bias_cc <- colMeans(M1_continuous_cc) - c(0.5, 0.5, -1, 2, -1, 1, 2.650375, 1.425591, 3.87515)
sd_cc <- apply(M1_continuous_cc, 2, sd) 

summary_M1_continuous <- data.frame(
  Bias_EM = mean_bias,
  Sd_EM = sd,
  Bias_cc = c(mean_bias_cc[1 : 6], NA, NA, mean_bias_cc[7 : 9]),
  Sd_cc = c(sd_cc[1 : 6], NA, NA, sd_cc[7 : 9])
) 








M2_binary <- Estimators_binary_y |> reduce(rbind) |>  as_tibble()
M2_binary <- M2_binary |> mutate(CRD = pnorm(theta_0 + theta_z + theta_x) * prop_x + pnorm(theta_0 + theta_z) * (1 - prop_x) -
                                   pnorm(theta_0 + theta_x) * prop_x - pnorm(theta_0) * (1 - prop_x), 
                                 CRD_1 = pnorm(theta_0 + theta_z + theta_x) -
                                   pnorm(theta_0 + theta_x),
                                 CRD_0 = pnorm(theta_0 + theta_z) -
                                   pnorm(theta_0))

mean_bias <- colMeans(M2_binary) - c(0.5, 0.5, -1, 2, 1, 1.5, -1.5, 1.5, 0.5799697, 0.4772499, 0.6826895)
sd <- apply(M2_binary, 2, sd) 

M2_binary_cc <- Estimators_binary_y_mle |> reduce(rbind) |>  as_tibble()
M2_binary_cc <- M2_binary_cc |> mutate(CRD = pnorm(theta_0 + theta_z + theta_x) * prop_x + pnorm(theta_0 + theta_z) * (1 - prop_x) -
                                         pnorm(theta_0 + theta_x) * prop_x - pnorm(theta_0) * (1 - prop_x), 
                                       CRD_1 = pnorm(theta_0 + theta_z + theta_x) -
                                         pnorm(theta_0 + theta_x),
                                       CRD_0 = pnorm(theta_0 + theta_z) -
                                         pnorm(theta_0))
mean_bias_cc <- colMeans(M2_binary_cc) - c(0.5, 0.5, -1, 2, 1, 0.5799697, 0.4772499, 0.6826895)
sd_cc <- apply(M2_binary_cc, 2, sd) 

summary_M2_binary <- data.frame(
  Bias = mean_bias,
  Sd = sd,
  Bias_cc = c(mean_bias_cc[1 : 5], NA, NA, NA, mean_bias_cc[6 : 8]),
  Sd_cc = c(sd_cc[1 : 5], NA, NA, NA, sd_cc[6 : 8])
) 


M2_continuous <- Estimators_continuous_y |> reduce(rbind) |>  as_tibble()
M2_continuous <- M2_continuous |> mutate(
  CRD = (exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) * prop_x + exp(theta_0 + theta_z + sigma_x ^ 2 / 2) * (1 - prop_x) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2) * prop_x -
           exp(theta_0 + sigma_x ^ 2 / 2) * (1 - prop_x)),
  CRD_1 =  exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2),
  CRD_0 =  exp(theta_0 + theta_z + sigma_x ^ 2 / 2) - exp(theta_0 + sigma_x ^ 2 / 2)
)

mean_bias <- colMeans(M2_continuous) - c(0.5, 0.5, -1, 2, -1, 1, 1.5, -2, 1, 2.650375, 1.425591, 3.87515)
sd <- apply(M2_continuous, 2, sd) 

M2_continuous_cc <- Estimators_continuous_y_mle |> reduce(rbind) |>  as_tibble()
M2_continuous_cc <- M2_continuous_cc |> mutate(
  CRD = (exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) * prop_x + exp(theta_0 + theta_z + sigma_x ^ 2 / 2) * (1 - prop_x) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2) * prop_x -
           exp(theta_0 + sigma_x ^ 2 / 2) * (1 - prop_x)),
  CRD_1 =  exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2),
  CRD_0 =  exp(theta_0 + theta_z + sigma_x ^ 2 / 2) - exp(theta_0 + sigma_x ^ 2 / 2)
)
mean_bias_cc <- colMeans(M2_continuous_cc) - c(0.5, 0.5, -1, 2, -1, 1, 2.650375, 1.425591, 3.87515)
sd_cc <- apply(M2_continuous_cc, 2, sd) 

summary_M2_continuous <- data.frame(
  Bias = mean_bias,
  Sd = sd,
  Bias_cc = c(mean_bias_cc[1 : 6], NA, NA, NA, mean_bias_cc[7 : 9]),
  Sd_cc = c(sd_cc[1 : 6], NA, NA, NA, sd_cc[7 : 9])
) 





M3_binary <- Estimators_binary_y |> reduce(rbind) |>  as_tibble()
M3_binary <- M3_binary |> mutate(CRD = pnorm(theta_0 + theta_z + theta_x) * prop_x + pnorm(theta_0 + theta_z) * (1 - prop_x) -
                                   pnorm(theta_0 + theta_x) * prop_x - pnorm(theta_0) * (1 - prop_x), 
                                 CRD_1 = pnorm(theta_0 + theta_z + theta_x) -
                                   pnorm(theta_0 + theta_x),
                                 CRD_0 = pnorm(theta_0 + theta_z) -
                                   pnorm(theta_0))

mean_bias <- colMeans(M3_binary) - c(0.5, 0.5, -1, 2, 1, 1, -1.5, 1.5, 0.5799697, 0.4772499, 0.6826895)
sd <- apply(M3_binary, 2, sd) 

M3_binary_cc <- Estimators_binary_y_mle |> reduce(rbind) |>  as_tibble()
M3_binary_cc <- M3_binary_cc |> mutate(CRD = pnorm(theta_0 + theta_z + theta_x) * prop_x + pnorm(theta_0 + theta_z) * (1 - prop_x) -
                                         pnorm(theta_0 + theta_x) * prop_x - pnorm(theta_0) * (1 - prop_x), 
                                       CRD_1 = pnorm(theta_0 + theta_z + theta_x) -
                                         pnorm(theta_0 + theta_x),
                                       CRD_0 = pnorm(theta_0 + theta_z) -
                                         pnorm(theta_0))
mean_bias_cc <- colMeans(M3_binary_cc) - c(0.5, 0.5, -1, 2, 1, 0.5799697, 0.4772499, 0.6826895)
sd_cc <- apply(M3_binary_cc, 2, sd) 

summary_M3_binary <- data.frame(
  Bias = mean_bias,
  Sd = sd,
  Bias_cc = c(mean_bias_cc[1 : 5], NA, NA, NA, mean_bias_cc[6 : 8]),
  Sd_cc = c(sd_cc[1 : 5], NA, NA, NA, sd_cc[6 : 8])
) 


M3_continuous <- Estimators_continuous_y |> reduce(rbind) |>  as_tibble()
M3_continuous <- M3_continuous |> mutate(
  CRD = exp(theta_0 + theta_z + theta_x * prop_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x * prop_x + sigma_x ^ 2 / 2),
  CRD_1 =  exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2),
  CRD_0 =  exp(theta_0 + theta_z + sigma_x ^ 2 / 2) - exp(theta_0 + sigma_x ^ 2 / 2)
)

mean_bias <- colMeans(M3_continuous) - c(0.5, 0.5, -1, 2, -1, 1, 2, -2, 1, 2.350402, 1.425591, 3.87515)
sd <- apply(M3_continuous, 2, sd) 

M3_continuous_cc <- Estimators_continuous_y_mle |> reduce(rbind) |>  as_tibble()
M3_continuous_cc <- M3_continuous_cc |> mutate(
  CRD = exp(theta_0 + theta_z + theta_x * prop_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x * prop_x + sigma_x ^ 2 / 2),
  CRD_1 =  exp(theta_0 + theta_z + theta_x + sigma_x ^ 2 / 2) - exp(theta_0 + theta_x + sigma_x ^ 2 / 2),
  CRD_0 =  exp(theta_0 + theta_z + sigma_x ^ 2 / 2) - exp(theta_0 + sigma_x ^ 2 / 2)
)
mean_bias_cc <- colMeans(M3_continuous_cc) - c(0.5, 0.5, -1, 2, -1, 1, 2.350402, 1.425591, 3.87515)
sd_cc <- apply(M3_continuous_cc, 2, sd) 

summary_M3_continuous <- data.frame(
  Bias = mean_bias,
  Sd = sd,
  Bias_cc = c(mean_bias_cc[1 : 6], NA, NA, NA, mean_bias_cc[7 : 9]),
  Sd_cc = c(sd_cc[1 : 6], NA, NA, NA, sd_cc[7 : 9])
) 




M5_continuous <- Estimators_continuous_y |> reduce(rbind) |>  as_tibble()
M5_continuous <- M5_continuous |> mutate(
  CRD = (theta_1z) * prop_x + (theta_0z) * (1 - prop_x),
  CRD_1 =   theta_1z,
  CRD_0 =   theta_0z
)

mean_bias <- colMeans(M5_continuous) - c(0.5, 0.5, 1, 2, -1, 1.5, 1, 2, 1.5, 1.5, 1, 1.75, 2, 1.5)
sd <- apply(M5_continuous, 2, sd) 

M5_continuous_cc <- Estimators_continuous_y_mle |> reduce(rbind) |>  as_tibble()
M5_continuous_cc <- M5_continuous_cc |> mutate(
  CRD = theta_1z * prop_x + theta_0z * (1 - prop_x),
  CRD_1 =  theta_1z,
  CRD_0 =  theta_0z
)
mean_bias_cc <- colMeans(M5_continuous_cc) - c(0.5, 0.5, 1, 2, -1, 1.5, 1, 2, 1.75, 2, 1.5)
sd_cc <- apply(M5_continuous_cc, 2, sd) 

summary_M5_continuous <- data.frame(
  Bias = mean_bias,
  Sd = sd,
  Bias_cc = c(mean_bias_cc[1 : 8], NA, NA, NA, mean_bias_cc[9 : 11]),
  Sd_cc = c(sd_cc[1 : 8], NA, NA, NA, sd_cc[9 : 11])
) 


M4_binary <- Estimators_binary_y |> reduce(rbind) |>  as_tibble()
CEE <- CE |> reduce(rbind) |>  as_tibble()
M4_binary <- cbind(M4_binary, CEE) |> as_tibble()

mean_bias <- colMeans(M4_binary) - c(0.5, 0.5, -1, 2, 1, 1, -2, 1.5, 0.3685667, 1.527505)
sd <- apply(M4_binary, 2, sd) 

M4_binary_cc <- Estimators_binary_y_mle |> reduce(rbind) |>  as_tibble()
M4_binary_cc <- M4_binary_cc |> mutate(CRR1_0 = pnorm(theta_0 + theta_z + theta_x) *  pnorm(theta_0) /
                                         (pnorm(theta_0 + theta_z) *  pnorm(theta_0 + theta_x)),
                                       COR1_0 = CRR1_0 / (pnorm(theta_0 + theta_z + theta_x, lower.tail = FALSE) *  pnorm(theta_0, lower.tail = FALSE) /
                                                   (pnorm(theta_0 + theta_z, lower.tail = FALSE) *  pnorm(theta_0 + theta_x, lower.tail = FALSE))))
mean_bias_cc <- colMeans(M4_binary_cc) - c(0.5, 0.5, -1, 2, 1, 0.3685667, 1.527505)
sd_cc <- apply(M4_binary_cc, 2, sd) 

summary_M4_binary <- data.frame(
  Bias = mean_bias,
  Sd = sd,
  Bias_cc = c(mean_bias_cc[1 : 5], NA, NA, NA, mean_bias_cc[6 : 7]),
  Sd_cc = c(sd_cc[1 : 5], NA, NA, NA, sd_cc[6 : 7])
) 


