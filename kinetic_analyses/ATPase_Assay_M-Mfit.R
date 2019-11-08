# library
library(tidyverse)
library(broom)

# load the data
dfone <- read.delim("~/research/hNS3_project/RNA_atpase_activity_assay/fitting_of_atpase_assay/bl21_expression_in_hepes/wt/rep1_wt_bl21_hepes_atpase_assay.dat", header=FALSE, comment.char="#")

dftwo <- read.delim("~/research/hNS3_project/RNA_atpase_activity_assay/fitting_of_atpase_assay/bl21_expression_in_hepes/wt/rep2_wt_bl21_hepes_atpase_assay.dat", header=FALSE, comment.char="#")

dfthree <- read.delim("~/research/hNS3_project/RNA_atpase_activity_assay/fitting_of_atpase_assay/bl21_expression_in_hepes/wt/rep3_wt_bl21_hepes_atpase_assay.dat", header=FALSE, comment.char="#")

# name the columns
colnames(dfone) <- c("time", 10, 25, 50, 100, 400, 1000)

colnames(dftwo) <- c("time", 10, 25, 50, 100, 400, 1000)

colnames(dfthree) <- c("time", 10, 25, 50, 100, 400, 1000)


# get initial rates for each data set individually for time < 400
kin3 <- dfthree %>%
  gather(conc, obs, 2:7) %>%
  filter(time < 400) %>%
  group_by(conc) %>%
  do(tidy(lm(obs~time, data=.)))

kin2 <- dftwo %>%
  gather(conc, obs, 2:7) %>%
  filter(time < 400) %>%
  group_by(conc) %>%
  do(tidy(lm(obs~time, data=.)))

kin1 <- dfone %>%
  gather(conc, obs, 2:7) %>%
  filter(time < 400) %>%
  group_by(conc) %>%
  do(tidy(lm(obs~time, data=.)))

# get rid of intercept values only keep the rates
rates3 <- kin3 %>%
  filter(estimate < 0.01)

rates2 <- kin2 %>%
  filter(estimate < 0.01)

rates1 <- kin1 %>%
  filter(estimate < 0.01)

# convert conc column to number
rates1$conc <- as.numeric(rates1$conc)

rates2$conc <- as.numeric(rates2$conc)

rates3$conc <- as.numeric(rates3$conc)

# calculate the M-M values for each data set
# fit the M-M equation
mmfit1 <- nls(estimate~(v*conc)/(k+conc*(1+conc/i)), data = rates1, start = list(v = 1e-4, k = 20, i = 500), trace = TRUE)

mmfit2 <- nls(estimate~(v*conc)/(k+conc*(1+conc/i)), data = rates2, start = list(v = 1e-4, k = 20, i = 500), trace = TRUE)

mmfit3 <- nls(estimate~(v*conc)/(k+conc*(1+conc/i)), data = rates3, start = list(v = 1e-4, k = 20, i = 500), trace = TRUE)

mmfittidy1 <- tidy(mmfit1)

mmfittidy2 <- tidy(mmfit2)

mmfittidy3 <- tidy(mmfit3)

# predict the fits
dpre1 <- tibble("conc" = seq(0, 2000, by = 10))

dpre1$estimate <- predict(mmfit1, newdata = dpre1)

# plot data 1
ggplot(rates1, aes(x = conc, y = estimate)) +
  geom_point() +
  geom_line(data = dpre1, aes(x = conc, y = estimate)) +
  geom_errorbar(aes(x = conc, ymin = estimate - std.error, ymax = estimate + std.error)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1e-3)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1300)) +
  theme_bw() +
  labs(x = "ATP concentration", y = "rate", caption = "dataset 1")

write_csv(mmfittidy1, "wt_data1_fit.csv")


# predict data 2
dpre2 <- tibble("conc" = seq(0, 2000, by = 10))

dpre2$estimate <- predict(mmfit2, newdata = dpre2)

# plot data 2
ggplot(rates2, aes(x = conc, y = estimate)) +
  geom_point() +
  geom_line(data = dpre2, aes(x = conc, y = estimate)) +
  geom_errorbar(aes(x = conc, ymin = estimate - std.error, ymax = estimate + std.error)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1e-3)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1300)) +
  theme_bw() +
  labs(x = "ATP concentration", y = "rate", caption = "dataset 2")

write_csv(mmfittidy2, "wt_data2_fit.csv")


# predict data 3
dpre3 <- tibble("conc" = seq(0, 2000, by = 10))

dpre3$estimate <- predict(mmfit3, newdata = dpre3)

# plot data 2
ggplot(rates3, aes(x = conc, y = estimate)) +
  geom_point() +
  geom_line(data = dpre3, aes(x = conc, y = estimate)) +
  geom_errorbar(aes(x = conc, ymin = estimate - std.error, ymax = estimate + std.error)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1e-3)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1300)) +
  theme_bw() +
  labs(x = "ATP concentration", y = "rate", caption = "dataset 3")

write_csv(mmfittidy3, "wt_data3_fit.csv")


# do a global fit of all three data sets together
# combine the date_frames
df_all <- bind_rows(rates1, rates2, rates3)

mmfit.all <- nls(estimate~(v*conc)/(k+conc*(1+conc/i)), data = df_all, start = list(v = 1e-4, k = 20, i = 500), trace = TRUE)

mmfittidyall <- tidy(mmfit.all)

# predict data 3
dpreall <- tibble("conc" = seq(0, 2000, by = 10))

dpreall$estimate <- predict(mmfit.all, newdata = dpreall)

# plot data 2
ggplot(df_all, aes(x = conc, y = estimate)) +
  geom_point() +
  geom_line(data = dpreall, aes(x = conc, y = estimate)) +
  #geom_errorbar(aes(x = conc, ymin = estimate - std.error, ymax = estimate + std.error)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1e-3)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1300)) +
  theme_bw() +
  labs(x = "ATP concentration", y = "rate", caption = "combined data")

write_csv(mmfittidyall, "wt_data_all_fit.csv")

