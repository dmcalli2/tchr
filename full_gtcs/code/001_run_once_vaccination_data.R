# 0001_vaccination data.frame.R

# Read in vaccination data
library(data.table)

vax <- readRDS("TOPLEVEL/Case_control/data/2021-07-28/CC_VACC_2021-07-28.rds")
vax <- as.data.table(vax)
vax <- vax[ , .(anon_id, vacc_occurence_time, vacc_product_name)]
vax[ , vax_date := as.Date(vacc_occurence_time)]

## drop duplicates
vax <- unique(vax)

## Drop vaccines before 8th December 2020 and non-applicable which are thought most likely to where had appointment but individual was were not vaccinated
vax <- vax[ vax_date >= as.Date("2020-12-08") &
          !vacc_product_name == "Not Applicable"&
          !anon_id == "ANONNANA",]

## two different vaccines on same day
setkey(vax, anon_id)

## Collapse multiple vaccine names into one
vax <- vax[, .(vacc_product_name = paste(vacc_product_name, collapse = "|")),
           by = .(anon_id, vax_date)]

## Order vaccines and count
vax <- setorder(vax, anon_id, vax_date)

## Explore vax where 3 or more doses 
vax3 <- (vax[, vax_seq := seq_len(.N),
             by = .(anon_id)]) 
vax3 <- vax3[vax_seq >=3,]
vax3 <- vax[anon_id %in% vax3$anon_id,]
vax3 <- vax3[ , .SD[c(1, .N)]  , anon_id]


## add back in to main vaccine data
vax <- vax [!anon_id %in% vax3$anon_id,]
vax <- rbind(vax, vax3)

# Assign as sequence one or two for first or second vaccine
vax <- setorder(vax, anon_id, vax_date)
vax <- vax[, vax_seq := seq_len(.N),
           by = .(anon_id)]
## Arrange vaccination data to wide format in other table
saveRDS(vax, "Data/vaccination_cleaned.Rds")
