#09_overall_counts_for_paper

library(tidyverse)

cc_data <- readRDS("Data/basefile.Rds")
cc_data <- cc_data %>% 
  filter(!(hcw_hhd == "hcw" & teacher_flag ==1))


teachers_control <- sum(!duplicated(cc_data$anon_id[cc_data$hcw_teacher == "Teacher" & cc_data$is_case ==0]))
hcws_control <- sum(!duplicated(cc_data$anon_id[cc_data$hcw_teacher == "Healthcare Worker" & cc_data$is_case ==0]))
tots_control <- sum(!duplicated(cc_data$anon_id[cc_data$is_case ==0 & cc_data$hcw_teacher == "Neither"]))
# c(teachers_control, hcws_control, tots_control)

teachers <- sum(!duplicated(cc_data$anon_id[cc_data$hcw_teacher == "Teacher"]))
hcws <- sum(!duplicated(cc_data$anon_id[cc_data$hcw_teacher == "Healthcare Worker, patient facing"]))
npfs <- sum(!duplicated(cc_data$anon_id[cc_data$hcw_teacher == "Healthcare worker, non-patient facing"]))

tots <- sum(!duplicated(cc_data$anon_id)) - teachers - hcws
hcw_hhd <- sum(!duplicated(cc_data$anon_id[cc_data$hhd ==1]))
tch_hhd <- sum(!duplicated(cc_data$anon_id[cc_data$hhd_tch ==1]))
hcw_tch_hhd <- sum(!duplicated(cc_data$anon_id[cc_data$hhd ==1 & cc_data$hhd_tch ==1]))
hcw_pop <- readRDS("output/number_hcws_hhds.Rds")


rm(cc_data)
save.image("output/overall_counts.Rdata")