pacman::p_load(
  rio,
  here,
  dlnm,
  splines,
  survival,
  survminer,
  coxme,
  sjPlot,
  tidyverse
)


data <- import(here("~/MGH/Thesis/Merged_IR+GPS_all 8 countries_Linet_20.03.2025.dta"))

data_subset <- data %>%
  select(
    unique_caseid,
    dhscc,
    anemia,
    pregnancy_date,
    climate_date,
    max_temp_monthly,
    interview_date,
    age_binary,
    v025,
    edu_binary,
    wealth_category,
    anc_binary,
    v481,
    depression_status,
    tobacco_consumption,
    ironTab_binary,
    ironTab_bin_duration,
    case_control#deworming Tab is missing & I assume this is the end-date
  )

km_data <- data_subset %>%
  group_by(unique_caseid) %>%
  # Count number of month since pregnancy_date, which is the number of row per id
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc)
  ) %>%
  ungroup() %>%
  # Keep only id with 2 months after and under 1 year postpartum -> lag = 2
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)


km_fit <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ dhscc,
  data = km_data
)


ggsurvplot(
  km_fit,
  data = km_data,
  fun = "cumhaz",
  pval = TRUE, # show p-value for log-rank test
  conf.int = FALSE, # confidence intervals
  risk.table = TRUE, # show number at risk below plot
  legend.title = "Country",
  xlab = "Months since pregnancy",
  ylab = "Survival probability (no delay)",
  palette = "Set1",# or use any ggplot2 color palette
)

#Step1: DLNM unadjusted and all countries together
library(dplyr)
library(survival)
library(dlnm)
library(sjPlot)

max_lag <- 21
temp_matrix <- data_subset %>%
  semi_join(km_data, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>% # Create lag variable starting at 0
  filter(n() <= max_lag) %>% # Keep only IDs with at least 'lag' rows
  slice(1:max_lag) %>% # Only keep the first 'lag' rows per ID
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()
  #as.matrix()

# Check the dimensions of temp_matrix
dim(temp_matrix)  # Should return 432 rows (unique_caseid count) and 21 columns (lag_0 to lag_20)
str(temp_matrix)  # To confirm number of columns


cb_temp <- crossbasis(
  temp_matrix,
  lag = 20, # max lag in months (example)
  argvar = list(fun = "ns", df = 3), # exposure-response basis
  arglag = list(fun = "ns", df = 2) # lag-response basis
)

# Check cb_temp for number of lags
dim(cb_temp)

model_cox_normal <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp + strata(dhscc), #  can add more other_covariates here
  data = km_data
)

sjPlot::tab_model(model_cox_normal)

# Create predictions from the crossbasis
pred <- crosspred(cb_temp, model_cox_normal, by = 0.1, cumul = TRUE)

colnames(km_data)



# Plot the overall cumulative effect (adjust parameters as needed)
plot(
  pred,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM effect of temperature on maternal anaemia."
)

#Plot 3D Graph of Temperature-Lag-HR
plot(
  pred,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR"
)

#Plot Contour of Temperature-Lag_Hazard Ratio
plot(
  pred,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR"
  ),
  key.title = title("Hazard Ratio")
)



#Step2: DLNM adjusted and all country
km_data_adj <- data_subset %>%
  group_by(unique_caseid) %>%
  # Count number of month since pregnancy_date, which is the number of row per id
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc),
    age_binary =  first(age_binary),
    v025 = first(v025),
    edu_binary = first(edu_binary),  # Ensure edu_binary is retained
    wealth_category = first(wealth_category),
    anc_binary = first(anc_binary),
    v481 = first(v481),
    tobacco_consumption = first(tobacco_consumption),
    depression_status = first(depression_status),
    ironTab_binary = first(ironTab_binary)
  ) %>%
  ungroup() %>%
  # Keep only id with 2 months after and under 1 year postpartum -> lag = 2
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)


km_fit_adj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ dhscc,
  data = km_data_adj
)


ggsurvplot(
  km_fit_adj,
  data = km_data_adj,
  fun = "cumhaz",
  pval = TRUE, # show p-value for log-rank test
  conf.int = FALSE, # confidence intervals
  risk.table = TRUE, # show number at risk below plot
  legend.title = "Country",
  xlab = "Months since pregnancy",
  ylab = "Survival probability (no delay)",
  palette = "Set1",# or use any ggplot2 color palette
)

library(dplyr)
library(survival)
library(dlnm)
library(sjPlot)

max_lag_adj <- 21
temp_matrix_adj <- data_subset %>%
  semi_join(km_data_adj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>% # Create lag variable starting at 0
  filter(n() <= max_lag_adj) %>% # Keep only IDs with at least 'lag' rows
  slice(1:max_lag_adj) %>% # Only keep the first 'lag' rows per ID
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()
#as.matrix()

# Check the dimensions of temp_matrix
dim(temp_matrix_adj)  # Should return 432 rows (unique_caseid count) and 21 columns (lag_0 to lag_20)
str(temp_matrix_adj)  # To confirm number of columns

cb_temp_adj <- crossbasis(
  temp_matrix_adj,
  lag = 20, # max lag in months (example)
  argvar = list(fun = "ns", df = 2), # exposure-response basis
  arglag = list(fun = "ns", df = 2) # lag-response basis
)

# Check cb_temp for number of lags
dim(cb_temp_adj)

colnames(km_data_adj)
#Adjusted cox model
model_cox_normal_adj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_adj + strata(dhscc) + age_binary + v025 + edu_binary + wealth_category + anc_binary + v481 + depression_status + tobacco_consumption + ironTab_binary,
  data = km_data_adj
)

sjPlot::tab_model(model_cox_normal_adj)

# Create predictions from the crossbasis
pred_adj <- crosspred(cb_temp_adj, model_cox_normal_adj, by = 0.1, cumul = TRUE)


# Plot the overall cumulative effect (adjust parameters as needed) adjusted
plot(
  pred_adj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM effect of temperature on maternal anaemia (adjusted)."
)

#Plot 3D Graph of Temperature-Lag-HR adjusted
plot(
  pred_adj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (adjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio adjusted
plot(
  pred_adj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (adjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step3: DLNM for Tanzania unadjusted
data_subset_TZ <- data %>%
  filter(dhscc == "TZ") %>%  # Keep only Tanzania
  select(
    unique_caseid,
    dhscc,
    anemia,
    pregnancy_date,
    climate_date,
    max_temp_monthly,
    interview_date,
    age_binary,
    v025,
    edu_binary,
    wealth_category,
    anc_binary,
    v481,
    depression_status,
    tobacco_consumption,
    ironTab_binary,
    ironTab_bin_duration,
    case_control
  )
km_data_TZ_uadj <- data_subset_TZ %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc)  # Already filtered for Tanzania
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_TZ_uadj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Tanzania)
  data = km_data_TZ_uadj
)


max_lag_TZ <- 21
temp_matrix_TZ_uadj <- data_subset_TZ %>%
  semi_join(km_data_TZ_uadj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_TZ) %>%
  slice(1:max_lag_TZ) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_TZ_uadj <- crossbasis(
  temp_matrix_TZ_uadj,
  lag = 20,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_TZ_uadj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_TZ_uadj,  # dhscc removed, since only Tanzania is included
  data = km_data_TZ_uadj
)

sjPlot::tab_model(model_cox_normal_TZ_uadj)
pred_TZ_uadj <- crosspred(cb_temp_TZ_uadj, model_cox_normal_TZ_uadj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for TZ unadjusted
plot(
  pred_TZ_uadj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Tanzania-unadjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for TZ unadjusted
plot(
  pred_TZ_uadj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Tanzania-unadjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for TZ unadjusted
plot(
  pred_TZ_uadj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Tanzania-unadjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step3.1 DLNM for Tanzania adjusted

km_data_TZ_adj <- data_subset_TZ %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc),
    age_binary =  first(age_binary),
    v025 = first(v025),
    edu_binary = first(edu_binary),  # Ensure edu_binary is retained
    wealth_category = first(wealth_category),
    anc_binary = first(anc_binary),
    v481 = first(v481),
    tobacco_consumption = first(tobacco_consumption),
    ironTab_binary = first(ironTab_binary)
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_TZ_adj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Tanzania)
  data = km_data_TZ_adj
)


max_lag_TZ_adj <- 21
temp_matrix_TZ_adj <- data_subset_TZ %>%
  semi_join(km_data_TZ_adj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_TZ_adj) %>%
  slice(1:max_lag_TZ_adj) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_TZ_adj <- crossbasis(
  temp_matrix_TZ_adj,
  lag = 20,
  argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_TZ_adj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_TZ_adj + age_binary + v025 + edu_binary + wealth_category + anc_binary + v481 + tobacco_consumption + ironTab_binary, # dhscc removed, since only Tanzania is included
  data = km_data_TZ_adj
)

sjPlot::tab_model(model_cox_normal_TZ_adj)
pred_TZ_adj <- crosspred(cb_temp_TZ_adj, model_cox_normal_TZ_adj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for TZ adjusted
plot(
  pred_TZ_adj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Tanzania adjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for TZ adjusted
plot(
  pred_TZ_adj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Tanzania-adjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for TZ adjusted
plot(
  pred_TZ_adj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Tanzania-adjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step4:DLNM for Nepal unadjusted
data_subset_NP <- data %>%
  filter(dhscc == "NP") %>%  # Keep only Tanzania
  select(
    unique_caseid,
    dhscc,
    anemia,
    pregnancy_date,
    climate_date,
    max_temp_monthly,
    interview_date,
    age_binary,
    v025,
    edu_binary,
    wealth_category,
    anc_binary,
    v481,
    depression_status,
    tobacco_consumption,
    ironTab_binary,
    ironTab_bin_duration,
    case_control
  )
km_data_NP_uadj <- data_subset_NP %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc)  # Already filtered for Tanzania
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_NP_uadj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Nepal)
  data = km_data_NP_uadj
)


max_lag_NP <- 21
temp_matrix_NP_uadj <- data_subset_NP %>%
  semi_join(km_data_NP_uadj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_NP) %>%
  slice(1:max_lag_NP) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_NP_uadj <- crossbasis(
  temp_matrix_NP_uadj,
  lag = 20,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_NP_uadj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_NP_uadj,  # dhscc removed, since only Tanzania is included
  data = km_data_NP_uadj
)

sjPlot::tab_model(model_cox_normal_NP_uadj)
pred_NP_uadj <- crosspred(cb_temp_NP_uadj, model_cox_normal_NP_uadj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for NP unadjusted
plot(
  pred_NP_uadj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Nepal-unadjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for NP unadjusted
plot(
  pred_NP_uadj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Nepal-unadjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for NP unadjusted
plot(
  pred_NP_uadj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Nepal-unadjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step4.1 DLNM for Nepal adjusted

km_data_NP_adj <- data_subset_NP %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc),
    age_binary =  first(age_binary),
    v025 = first(v025),
    edu_binary = first(edu_binary),  # Ensure edu_binary is retained
    wealth_category = first(wealth_category),
    anc_binary = first(anc_binary),
    v481 = first(v481),
    tobacco_consumption = first(tobacco_consumption),
    ironTab_binary = first(ironTab_binary)
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_NP_adj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Nepal)
  data = km_data_NP_adj
)


max_lag_NP_adj <- 21
temp_matrix_NP_adj <- data_subset_NP %>%
  semi_join(km_data_NP_adj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_NP_adj) %>%
  slice(1:max_lag_NP_adj) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_NP_adj <- crossbasis(
  temp_matrix_NP_adj,
  lag = 20,
  argvar = list(fun = "ns", df = 1),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_NP_adj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_NP_adj + age_binary + v025 + edu_binary + wealth_category + anc_binary + v481 + tobacco_consumption + ironTab_binary, # dhscc removed, since only Tanzania is included
  data = km_data_NP_adj
)

sjPlot::tab_model(model_cox_normal_NP_adj)
pred_NP_adj <- crosspred(cb_temp_NP_adj, model_cox_normal_NP_adj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for NP adjusted
plot(
  pred_NP_adj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Nepal adjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for NP adjusted
plot(
  pred_NP_adj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Nepal-adjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for NP adjusted
plot(
  pred_NP_adj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Nepal-adjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step5: DLNM for Mozambique unadjusted
data_subset_MZ <- data %>%
  filter(dhscc == "MZ") %>%  # Keep only Mozambique
  select(
    unique_caseid,
    dhscc,
    anemia,
    pregnancy_date,
    climate_date,
    max_temp_monthly,
    interview_date,
    age_binary,
    v025,
    edu_binary,
    wealth_category,
    anc_binary,
    v481,
    depression_status,
    tobacco_consumption,
    ironTab_binary,
    ironTab_bin_duration,
    case_control
  )
km_data_MZ_uadj <- data_subset_MZ %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc)  # Already filtered for MZ
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_MZ_uadj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Mozambique)
  data = km_data_MZ_uadj
)


max_lag_MZ <- 21
temp_matrix_MZ_uadj <- data_subset_MZ %>%
  semi_join(km_data_MZ_uadj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_MZ) %>%
  slice(1:max_lag_MZ) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_MZ_uadj <- crossbasis(
  temp_matrix_MZ_uadj,
  lag = 20,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_MZ_uadj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_MZ_uadj,  # dhscc removed, since only Mozambique is included
  data = km_data_MZ_uadj
)

sjPlot::tab_model(model_cox_normal_MZ_uadj)
pred_MZ_uadj <- crosspred(cb_temp_MZ_uadj, model_cox_normal_MZ_uadj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for MZ unadjusted
plot(
  pred_MZ_uadj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Mozambique-unadjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for MZ unadjusted
plot(
  pred_MZ_uadj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Mozambique-unadjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for MZ unadjusted
plot(
  pred_MZ_uadj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Mozambique-unadjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step5.1 DLNM for Mozambique adjusted

km_data_MZ_adj <- data_subset_MZ %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc),
    age_binary =  first(age_binary),
    v025 = first(v025),
    edu_binary = first(edu_binary),  # Ensure edu_binary is retained
    wealth_category = first(wealth_category),
    anc_binary = first(anc_binary),
    v481 = first(v481),
    depression_status = first(depression_status),
    tobacco_consumption = first(tobacco_consumption),
    ironTab_binary = first(ironTab_binary)
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_MZ_adj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Mozambique)
  data = km_data_MZ_adj
)


max_lag_MZ_adj <- 21
temp_matrix_MZ_adj <- data_subset_MZ %>%
  semi_join(km_data_MZ_adj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_MZ_adj) %>%
  slice(1:max_lag_MZ_adj) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_MZ_adj <- crossbasis(
  temp_matrix_MZ_adj,
  lag = 20,
  argvar = list(fun = "ns", df = 2),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_MZ_adj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_MZ_adj + age_binary + v025 + edu_binary + wealth_category + anc_binary + v481 + depression_status + tobacco_consumption + ironTab_binary, # dhscc removed, since only Tanzania is included
  data = km_data_MZ_adj
)

sjPlot::tab_model(model_cox_normal_MZ_adj)
pred_MZ_adj <- crosspred(cb_temp_MZ_adj, model_cox_normal_MZ_adj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for MZ adjusted
plot(
  pred_MZ_adj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Mozambique adjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for MZ adjusted
plot(
  pred_MZ_adj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Mozambique-adjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for MZ adjusted
plot(
  pred_MZ_adj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Mozambique-adjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step6: DLNM for Lesotho unadjusted
data_subset_LS <- data %>%
  filter(dhscc == "LS") %>%  # Keep only Lesotho
  select(
    unique_caseid,
    dhscc,
    anemia,
    pregnancy_date,
    climate_date,
    max_temp_monthly,
    interview_date,
    age_binary,
    v025,
    edu_binary,
    wealth_category,
    anc_binary,
    v481,
    depression_status,
    tobacco_consumption,
    ironTab_binary,
    ironTab_bin_duration,
    case_control
  )
km_data_LS_uadj <- data_subset_LS %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc)  # Already filtered for LS
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_LS_uadj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Lesotho)
  data = km_data_LS_uadj
)


max_lag_LS <- 21
temp_matrix_LS_uadj <- data_subset_LS %>%
  semi_join(km_data_LS_uadj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_LS) %>%
  slice(1:max_lag_LS) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_LS_uadj <- crossbasis(
  temp_matrix_LS_uadj,
  lag = 20,
  argvar = list(fun = "ns", df = 1),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_LS_uadj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_LS_uadj,  # dhscc removed, since only Lesotho is included
  data = km_data_LS_uadj
)

sjPlot::tab_model(model_cox_normal_LS_uadj)
pred_LS_uadj <- crosspred(cb_temp_LS_uadj, model_cox_normal_LS_uadj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for LS unadjusted
plot(
  pred_LS_uadj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Lesotho-unadjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for LS unadjusted
plot(
  pred_LS_uadj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Lesotho-unadjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for LS unadjusted
plot(
  pred_LS_uadj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Lesotho-unadjusted)"
  ),
  key.title = title("Hazard Ratio")
)

#Step6.1 DLNM for Lesotho adjusted

km_data_LS_adj <- data_subset_LS %>%
  group_by(unique_caseid) %>%
  summarise(
    eventmonth_since_pregnancy = n(),
    anemia = max(anemia),
    dhscc = first(dhscc),
    age_binary =  first(age_binary),
    v025 = first(v025),
    edu_binary = first(edu_binary),  # Ensure edu_binary is retained
    wealth_category = first(wealth_category),
    anc_binary = first(anc_binary),
    v481 = first(v481),
    depression_status = first(depression_status),
    tobacco_consumption = first(tobacco_consumption),
    ironTab_binary = first(ironTab_binary)
  ) %>%
  ungroup() %>%
  filter(eventmonth_since_pregnancy >= 2 & eventmonth_since_pregnancy <= 21)

km_fit_LS_adj <- survfit(
  Surv(eventmonth_since_pregnancy, anemia) ~ 1,  # No need for dhscc (only Lesotho)
  data = km_data_LS_adj
)


max_lag_LS_adj <- 21
temp_matrix_LS_adj <- data_subset_LS %>%
  semi_join(km_data_LS_adj, by = "unique_caseid") %>%
  group_by(unique_caseid) %>%
  mutate(lag = row_number() - 1) %>%
  filter(n() <= max_lag_LS_adj) %>%
  slice(1:max_lag_LS_adj) %>%
  ungroup() %>%
  select(unique_caseid, lag, max_temp_monthly) %>%
  pivot_wider(
    names_from = lag,
    values_from = max_temp_monthly,
    names_prefix = "lag"
  ) %>%
  select(-unique_caseid) %>%
  as.data.frame()

cb_temp_LS_adj <- crossbasis(
  temp_matrix_LS_adj,
  lag = 20,
  argvar = list(fun = "ns", df = 1),
  arglag = list(fun = "ns", df = 2)
)

model_cox_normal_LS_adj <- coxph(
  Surv(eventmonth_since_pregnancy, anemia) ~ cb_temp_LS_adj + age_binary + v025 + edu_binary + wealth_category + anc_binary + v481 + depression_status + tobacco_consumption + ironTab_binary, # dhscc removed, since only Tanzania is included
  data = km_data_LS_adj
)

sjPlot::tab_model(model_cox_normal_LS_adj)
pred_LS_adj <- crosspred(cb_temp_LS_adj, model_cox_normal_LS_adj, by = 0.1, cumul = TRUE)

# Plot the overall cumulative effect for LS adjusted
plot(
  pred_LS_adj,
  exp = TRUE,
  "overall",
  xlab = "Temperature (°C)",
  ylab = "Hazard Ratio",
  main = "Cumulative DLNM Effect of Temperature on Maternal Anemia (Lesotho adjusted)"
)


#Plot 3D Graph of Temperature-Lag-HR for LS adjusted
plot(
  pred_LS_adj,
  exp = TRUE,
  xlab = "Temperature (°C)",
  zlab = "Hazard Ratio",
  theta = 220, # rotation angle
  phi = 35, # elevation angle
  ltheta = -120, # lighting angle
  main = "3D Graph of Temperature-Lag-HR (Lesotho-adjusted)"
)

#Plot Contour of Temperature-Lag_Hazard Ratio for LS adjusted
plot(
  pred_LS_adj,
  exp = TRUE,
  "contour",
  plot.title = title(
    xlab = "Temperature (°C)",
    ylab = "Lag (months)",
    main = "Contour Plot of Temperature-Lag-HR (Lesotho-adjusted)"
  ),
  key.title = title("Hazard Ratio")
)
