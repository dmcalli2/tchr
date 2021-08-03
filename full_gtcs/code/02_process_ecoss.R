#02_process__ecoss

## Drops all tests subsequent to first positive test. Also drop results which are mis-dated

library(dplyr)
library(janitor)

### obtain testing summaries ----
ecoss <-  readRDS(
  "TOPLEVEL/Case_control/data/2021-07-28/CC_ecoss_tests__2021-07-28.rds") %>%
  clean_names() %>% 
  # filter(test_type == "PCR") %>% 
  select(anon_id, test_date = date_ecoss_specimen, ncov_result = test_result)  %>% 
  mutate(test_date = as.Date(test_date))

## drop ecoss before 2020
ecoss <- ecoss %>% 
  mutate(test_date = as.Date(test_date)) %>% 
  filter(test_date >= as.Date("2020-01-01")) %>% 
  mutate(ncov_result = case_when(
    ncov_result == "NEGATIVE" ~ "Negative",
    ncov_result == "POSITIVE" ~ "Positive",
    TRUE ~ ncov_result
  ))

# Drop all tests after first positive test
# First create dataset with first positive test date
firstpos <- ecoss %>% 
  filter(ncov_result == "Positive") %>% 
  arrange(anon_id, test_date) %>% 
  select(anon_id, test_date) %>% 
  distinct(anon_id, .keep_all = TRUE) %>% 
  rename(first_pos = test_date) 
ecoss <- ecoss %>% 
  left_join(firstpos) %>% 
  mutate(first_pos = if_else(is.na(first_pos), as.Date("9999-01-01"), first_pos)) 
ecoss %>% 
  filter(test_date > first_pos) %>% 
  summarise(tests = length(anon_id),
            people = sum(!duplicated(anon_id)))
## drop any tests after the first positive or on the same day as the first positive which is not itself positive
ecoss <- ecoss %>% 
  filter(test_date < first_pos  |
           (test_date == first_pos & ncov_result == "Positive")) %>% 
  distinct(anon_id, test_date, ncov_result)

rm(firstpos)

