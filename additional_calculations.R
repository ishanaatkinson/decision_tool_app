

##### READ IN DATA #####

country_code <- "CMR"

full_data <- readRDS(paste0(getwd(), "/", country_code, ".rds"))

# convert all "-" and spaces into "_" to ease in analysis later on and remove any accents from admin-1 area names
full_data$sites$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$sites$name_1), id = "Latin-ASCII")
full_data$prevalence$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$prevalence$name_1), id = "Latin-ASCII")
full_data$interventions$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$interventions$name_1), id = "Latin-ASCII")
full_data$population$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$population$name_1), id = "Latin-ASCII")
full_data$vectors$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$vectors$name_1), id = "Latin-ASCII")
full_data$pyrethroid_resistance$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$pyrethroid_resistance$name_1), id = "Latin-ASCII")
full_data$seasonality$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$seasonality$name_1), id = "Latin-ASCII")
full_data$eir$name_2 <- stri_trans_general(str=gsub("-", "_", full_data$eir$name_1), id = "Latin-ASCII")

full_data$sites$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$sites$name_2), id = "Latin-ASCII")
full_data$prevalence$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$prevalence$name_2), id = "Latin-ASCII")
full_data$interventions$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$interventions$name_2), id = "Latin-ASCII")
full_data$population$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$population$name_2), id = "Latin-ASCII")
full_data$vectors$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$vectors$name_2), id = "Latin-ASCII")
full_data$pyrethroid_resistance$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$pyrethroid_resistance$name_2), id = "Latin-ASCII")
full_data$seasonality$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$seasonality$name_2), id = "Latin-ASCII")
full_data$eir$name_2 <- stri_trans_general(str=gsub(" ", "_", full_data$eir$name_2), id = "Latin-ASCII")

area_names <- unique(full_data$sites$name_2)



##### COMMON PARAMETERS/VARIABLES #####

human_population <- 1e6 # population size in model
min_age_to_model <- 0 # minimum age to model (in days)
max_age_to_model <- 2.5*365 # enter in multiples of 0.5, maximum age to model (in days)

# set the time span over which to simulate
# NOTE: currently 23 years is the max as we only have ITN/IRS data for this length
year <- 365; years <- 1; sim_length_for_data <- year * years

years_proj_forward <- 1 # number of years to project forward following known data
years_of_simulation <- years + years_proj_forward # simulation length (years)
sim_length <- (years + years_proj_forward) * year # simulation length (days)

# interval between modeled age groups (1 = days, 7 = weeks etc)
step_length <- 1

# vector of min and max ages to be used in each age bracket
age_min <- seq(min_age_to_model, max_age_to_model, step_length)
age_max <- seq(min_age_to_model, max_age_to_model, step_length) + step_length

# six month ages (in days)
sixmonth_intervals <- c(0, 183, 365, 548, 730, 913, 1095, 1278, 1460, 1643,
                        1825, 2008, 2190, 2373, 2555, 2738, 2920, 3285, 3468,
                        3650, 4015)

# time steps used in model (number of rows in the simulations output data frame)
timesteps<-seq(0,round(sim_length),1)

# vector for the midpoint age for each age bracket (useful for graphing)
age_in_days <- seq(age_min[1], age_max[length(age_max)], length.out=length(age_min) + 1)
age_in_days_midpoint <- age_in_days[-length(age_in_days)] + diff(age_in_days)/2

# number of 6 month interval age groups to model
no_sixmonth_intervals <- round(max(age_max)/182.5)

# empty vectors for the column names of interest
age_group_names <- c() # age groups
age_group_names_sixmonth <- c() # 6 month age groups
sixmonth_intervals_midpoint <- c() # 6 month age group midpoints
clin_inc_cols <- c() # clinical incidence column names
sev_inc_cols <- c() # severe incidence column names
tot_inc_cols <- c() # total incidence column names
asym_inc_cols <- c() # asymptomatic incidence column names

# fill empty vectors
for (i in 1:(length(age_min))) {
  age_group_names <- append(age_group_names, paste0("n_age_", as.character(age_min[i]), "_", as.character(age_max[i])))
}

for (i in 1:(length(sixmonth_intervals) - 1)) {
  
  if (i == (length(sixmonth_intervals) - 1)) {
    age_group_names_sixmonth <- append(age_group_names_sixmonth, paste0("n_age_", as.character(sixmonth_intervals[i]), "_", as.character(sixmonth_intervals[i+1])))
    sixmonth_intervals_midpoint <- append(sixmonth_intervals_midpoint, (sixmonth_intervals[i]+sixmonth_intervals[i+1])/2)
  }
  
  if (i != (length(sixmonth_intervals) - 1)) {
    age_group_names_sixmonth <- append(age_group_names_sixmonth, paste0("n_age_", as.character(sixmonth_intervals[i]), "_", as.character(sixmonth_intervals[i+1] - 1)))
    sixmonth_intervals_midpoint <- append(sixmonth_intervals_midpoint, (sixmonth_intervals[i]+(sixmonth_intervals[i+1] - 1))/2)
  }
  
}

# write column names
for (i in 1:length(age_min)) {
  clin_inc_cols <- append(clin_inc_cols, paste0("n_inc_clinical_", as.character(age_min[i]), "_", as.character(age_max[i])))
  sev_inc_cols <- append(sev_inc_cols, paste0("n_inc_severe_", as.character(age_min[i]), "_", as.character(age_max[i])))
  tot_inc_cols <- append(tot_inc_cols, paste0("n_inc_", as.character(age_min[i]), "_", as.character(age_max[i])))
  asym_inc_cols <- append(asym_inc_cols, paste0("n_inc_asym_", as.character(age_min[i]), "_", as.character(age_max[i])))
}

##### READ IN IMPERIAL MODEL OUTPUT #####

incidence_ppy_df <- read.csv(paste0(getwd(), "/simulation_results/incidence_ppy_df_merged.csv"))
population_df <- read.csv(paste0(getwd(), "/simulation_results/population_df_merged.csv"))

##### READ IN PROPORTIONS OF HAPLOTYPES #####

haplotype_data <- read_xlsx(paste0(getwd(), "/haplotypes_combined_forIshana2.xlsx"))

# converts admin1 names into clean format that R can understand
haplotype_data$NAME_2 <- stri_trans_general(str=gsub("-", "_", haplotype_data$NAME_1), id = "Latin-ASCII")
haplotype_data$NAME_2 <- stri_trans_general(str=gsub(" ", "_", haplotype_data$NAME_2), id = "Latin-ASCII")

haplotype_proportions <- data.frame()

I_AKA_ <- I_GKA_ <- I_GEA_ <- I_GEG_ <- V_GKA_ <- V_GKG_ <- c()

# loops over each country and admin1 area and calculates proportions for areas
# which have 6 codon data and appends to data frame

country_haplotype_data <- haplotype_data %>% filter(iso_code == country_code)
admin1 <- unique(country_haplotype_data$NAME_2)

# start using admin1 over area_names so we only calculate impact for where
# we have haplotype data

for (j in 1:length(admin1)) {
  
  # area haplotype data
  hap <- country_haplotype_data %>% filter(NAME_2 == admin1[j])
  
  denom <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS,
                            hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS,
                            hap$ISGEAA, hap$IAGEAA, hap$ISGEAS,
                            hap$ISGEGA, hap$IAGEGA,
                            hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS,
                            hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA,
                            hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA
  )), na.rm=TRUE)
  
  I_AKA_ <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  I_GKA_ <- sum(as.integer(c(hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  I_GEA_ <- sum(as.integer(c(hap$ISGEAA, hap$IAGEAA, hap$ISGEAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  I_GEG_ <- sum(as.integer(c(hap$ISGEGA, hap$IAGEGA)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  V_GKA_ <- sum(as.integer(c(hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  V_GKG_ <- sum(as.integer(c(hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  other <- sum(as.integer(c(hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA)), na.rm=TRUE) / denom#sum(as.integer(hap$denom_haplo))
  
  props_5sf <- round(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other), digits=5)
  sum_of_prop <- sum(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other))
  
  # check proportions add to 1
  if (is.nan(sum_of_prop) == TRUE) {
    print(paste0("Summation of proportions for ", admin1[j], " in ", country_code, " cannot be calculated due to lack of 6-codon data"))
  } else {
    print(paste0("Summation of proportions for ", admin1[j], " in ", country_code, " is: ", sum_of_prop))
    
    haplotype_proportions <- rbind(haplotype_proportions, c(country_code, admin1[j], mean(hap$longitude, na.rm = TRUE), mean(hap$latitude, na.rm = TRUE),  props_5sf, sum_of_prop))
  }
  
  
  
}


# change column names
colnames(haplotype_proportions) <- c("iso_code", "NAME_2", "longitude", "latitude", "I_AKA_", "I_GKA_", "I_GEA_", "I_GEG_", "V_GKA_", "V_GKG_", "other", "Sum of proportions" )



##### CALCULATE INCIDENCE FOR WHOLE COUNTRY #####

# NOTE: population weights are calculated out of the admin-1 units that have
# haplotype data so admin-1 units that have no data take the average

admin1_with_haplotype_data <- c()

for (i in 1:length(admin1)){
  proportions <- haplotype_proportions %>% filter(iso_code == country_code, NAME_2 == admin1[i])
  
  if (dim(proportions)[1] != 0) {
    
    admin1_with_haplotype_data <- c(admin1_with_haplotype_data, admin1[i])
    
  }
  
}

admin1 <- admin1_with_haplotype_data

incidence_ppy_df_whole_country_weighted <- (incidence_ppy_df %>% filter(area %in% admin1))
population_weights <- c()

pop_size_across_admin1 <- sum((full_data$population %>% filter(year==2023, name_2 %in% admin1))$pop)

for (i in 1:length(admin1)){
  
  pop <- sum((full_data$population %>% filter(year == 2023, name_2 == admin1[i]))$pop)
  
  pop_weight <- pop / pop_size_across_admin1
  population_weights <- c(population_weights, pop_weight)
}

population_weights <- rep(population_weights, each=4*length(age_in_days_midpoint))
incidence_ppy_df_whole_country_weighted$value <- incidence_ppy_df_whole_country_weighted$value * population_weights

incidence_ppy_df_whole_country <- data.frame(age_in_days_midpoint)
incidence_ppy_df_whole_country_age_sum_clin <- incidence_ppy_df_whole_country_age_sum_sev <- incidence_ppy_df_whole_country_age_sum_tot <- incidence_ppy_df_whole_country_age_sum_asym <- c()

SP_protection_additional_doses_whole_country <- c()
SP_protection_whole_country <- c()

for (i in 1:length(age_in_days_midpoint)) {
  
  # incidence 
  incidence_ppy_df_whole_country_age_sum_clin <- c(incidence_ppy_df_whole_country_age_sum_clin, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "clinical", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  incidence_ppy_df_whole_country_age_sum_sev <- c(incidence_ppy_df_whole_country_age_sum_sev, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "severe", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  incidence_ppy_df_whole_country_age_sum_tot <- c(incidence_ppy_df_whole_country_age_sum_tot, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "total", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  incidence_ppy_df_whole_country_age_sum_asym <- c(incidence_ppy_df_whole_country_age_sum_asym, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "asymptomatic", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  
}

incidence_ppy_df_whole_country$clinical <- incidence_ppy_df_whole_country_age_sum_clin
incidence_ppy_df_whole_country$severe <- incidence_ppy_df_whole_country_age_sum_sev
incidence_ppy_df_whole_country$total <- incidence_ppy_df_whole_country_age_sum_tot
incidence_ppy_df_whole_country$asymptomatic <- incidence_ppy_df_whole_country_age_sum_asym

##### WRITE CSVs #####

write.csv(incidence_ppy_df_whole_country, paste0(getwd(), "/simulation_results/incidence_ppy_df_whole_country.csv"), row.names=FALSE)

