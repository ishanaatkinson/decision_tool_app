
#schedule <- c(10*7, 14*7, 6*30, 9*30, 15*30)
schedule_additional_doses <- c(10*7, 14*7, 6*30, 9*30, 12*30, 15*30, 18*30, 24*30)
#cov <- c(0.82, 0.73, 0.26, 0.52, 0.24)
cov_additional_doses <- c(0.82, 0.73, 0.26, 0.52, 0.15 , 0.24 ,0.09, 0.03)

schedule <- 7*schedule # convert weeks into days


##### SET UP VECTORS FOR PMC SCHEDULE AND DOSES #####


# Get the order of indices that would sort schedule_additional_doses
sorted_indices <- order(schedule)
sorted_indices_additional_doses <- order(schedule_additional_doses)

# Use sorted_indices to reorder both vectors
sorted_schedule <- schedule[sorted_indices]
sorted_cov <- cov[sorted_indices]

sorted_schedule_additional_doses <- schedule_additional_doses[sorted_indices_additional_doses]
sorted_cov_additional_doses <- cov_additional_doses[sorted_indices_additional_doses]

schedule <- sorted_schedule
schedule_additional_doses <- sorted_schedule_additional_doses
cov <- sorted_cov
cov_additional_doses <- sorted_cov_additional_doses


##### READ IN DATA #####
full_data <- readRDS(paste0(getwd(), "/", country, "/", country, ".rds"))

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


list_of_possible_country_codes <- c("CMR", "MOZ")
area_names <- unique(full_data$sites$name_2)

##### PRINT ERROR MESSAGES #####
if (length(schedule) != length(cov)) {
  print("Number of PMC doses does not equal number of coverage values inputted.")
}

if (length(schedule_additional_doses) != length(cov_additional_doses)) {
  print("Number of PMC doses for the additional schedule does not equal number of coverage values inputted.")
}

if (!(country %in% list_of_possible_country_codes)) {
  print("Country code is unknown.")
}

# if (!(chosen_area %in% area_names)) {
#   paste0("Chosen area '", chosen_area, "' is not an offical admin-1 unit in ", country, ".")
# }

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

incidence_ppy_df <- read.csv(paste0(getwd(), "/", country, "/simulation_results/incidence_ppy_df_merged.csv"))
population_df <- read.csv(paste0(getwd(), "/", country, "/simulation_results/population_df_merged.csv"))
incidence_ppy_df_whole_country <- read.csv(paste0(getwd(), "/", country, "/simulation_results/incidence_ppy_df_whole_country.csv"))


##### READ IN PROPORTIONS OF HAPLOTYPES #####

source(paste0(getwd(), "/", country, "/HAPLOTYPE_PROPORTIONS.R"), local=TRUE)

country_haplotype_data <- haplotype_data %>% filter(iso_code == country)
admin1 <- unique(country_haplotype_data$NAME_2)


##### EFFICACY CURVES WEIBULL PARAMETERS #####

# Weibull scale parameters for each haplotype
lambda_trip<-59.57659
lambda_quadr<-33.05391
lambda_quint<-18.55328
lambda_sext<-12.81186
lambda_other<-23
lambda_VAGKAA<-18.55328 # assumed to be same as the QUINTUPLE
lambda_VAGKGS<-12.81186 # assumed to be same as the SEXTUPLE

# Weibull shape parameters for each haplotype
w_trip<- 8.435971
w_quadr<-4.862126
w_quint<-2.4840752
w_sext<-3.691953
w_other<-4.5
w_VAGKAA<-1.7
w_VAGKGS<- 4.1


##### SET UP EFFICACY CURVED (NORMAL PMC SCHEDULE) #####

# set up df with the timesteps
if (sim_length > length(age_in_days_midpoint)) {
  df<- data.frame(time=1:sim_length)
}

if (sim_length < length(age_in_days_midpoint)) {
  df<- data.frame(time=1:length(age_in_days_midpoint))
}

# initalise empty columns
df$prot_trip<-NA
df$prot_quadr<-NA
df$prot_quint<-NA
df$prot_sext<-NA
df$prot_other<-NA
df$prot_VAGKAA<-NA
df$prot_VAGKGS<-NA


# construct weibull curves to model the efficacy of SP through time for
# each haplotype

if (length(schedule) > 0) {
  for (t in 1: schedule[1]) {  # day 0 to dose 1 on day 70
    df$prot_trip[t]<-0
    df$prot_quadr[t]<-0
    df$prot_quint[t]<-0
    df$prot_sext[t]<- 0
    df$prot_other[t]<-0
    df$prot_VAGKAA[t]<- 0
    df$prot_VAGKGS[t]<- 0
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 1) {
  for (t in (schedule[1]+1) : schedule[2])  {   # day 71 to dose 2 on day 98
    df$prot_trip[t]<- exp(-(df$time[t-schedule[1]]/lambda_trip)^w_trip) * cov[1]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[1]]/lambda_quadr)^w_quadr)* cov[1]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[1]]/lambda_quint)^w_quint)* cov[1]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[1]]/lambda_sext)^w_sext)* cov[1]
    df$prot_other[t]<- exp(-(df$time[t-schedule[1]]/lambda_other)^w_other)* cov[1]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[1]]/lambda_VAGKAA)^w_VAGKAA)* cov[1]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[1]]/lambda_VAGKGS)^w_VAGKGS)* cov[1]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 2) {
  for (t in (schedule[2]+1) : schedule[3])  {  #  day 99 to dose 3 on day 180
    df$prot_trip[t]<- exp(-(df$time[t-schedule[2]]/lambda_trip)^w_trip)* cov[2]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[2]]/lambda_quadr)^w_quadr)* cov[2]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[2]]/lambda_quint)^w_quint)* cov[2]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[2]]/lambda_sext)^w_sext)* cov[2]
    df$prot_other[t]<- exp(-(df$time[t-schedule[2]]/lambda_other)^w_other)* cov[2]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[2]]/lambda_VAGKAA)^w_VAGKAA)* cov[2]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[2]]/lambda_VAGKGS)^w_VAGKGS)* cov[2]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 3) {
  for (t in (schedule[3]+1) : schedule[4])  {  # day 181 to dose 4 on day 270
    df$prot_trip[t]<- exp(-(df$time[t-schedule[3]]/lambda_trip)^w_trip)* cov[3]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[3]]/lambda_quadr)^w_quadr)* cov[3]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[3]]/lambda_quint)^w_quint)* cov[3]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[3]]/lambda_sext)^w_sext)* cov[3]
    df$prot_other[t]<- exp(-(df$time[t-schedule[3]]/lambda_other)^w_other)* cov[3]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[3]]/lambda_VAGKAA)^w_VAGKAA)* cov[3]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[3]]/lambda_VAGKGS)^w_VAGKGS)* cov[3]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 4) {
  for (t in (schedule[4]+1) : schedule[5])  {  # day 271 to dose 5 on day 360
    df$prot_trip[t]<- exp(-(df$time[t-schedule[4]]/lambda_trip)^w_trip)* cov[4]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[4]]/lambda_quadr)^w_quadr)* cov[4]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[4]]/lambda_quint)^w_quint)* cov[4]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[4]]/lambda_sext)^w_sext)* cov[4]
    df$prot_other[t]<- exp(-(df$time[t-schedule[4]]/lambda_other)^w_other)* cov[4]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[4]]/lambda_VAGKAA)^w_VAGKAA)* cov[4]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[4]]/lambda_VAGKGS)^w_VAGKGS)* cov[4]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 5) {
  for (t in (schedule[5]+1) : schedule[6])  {  # day 361 to dose 6 on day 450
    df$prot_trip[t]<- exp(-(df$time[t-schedule[5]]/lambda_trip)^w_trip)* cov[5]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[5]]/lambda_quadr)^w_quadr)* cov[5]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[5]]/lambda_quint)^w_quint)* cov[5]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[5]]/lambda_sext)^w_sext)* cov[5]
    df$prot_other[t]<- exp(-(df$time[t-schedule[5]]/lambda_other)^w_other)* cov[5]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[5]]/lambda_VAGKAA)^w_VAGKAA)* cov[5]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[5]]/lambda_VAGKGS)^w_VAGKGS)* cov[5]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 6) {
  for (t in (schedule[6]+1) : schedule[7])  {  # day 451 to dose 7 on day 540
    df$prot_trip[t]<- exp(-(df$time[t-schedule[6]]/lambda_trip)^w_trip)* cov[6]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[6]]/lambda_quadr)^w_quadr)* cov[6]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[6]]/lambda_quint)^w_quint)* cov[6]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[6]]/lambda_sext)^w_sext)* cov[6]
    df$prot_other[t]<- exp(-(df$time[t-schedule[6]]/lambda_other)^w_other)* cov[6]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[6]]/lambda_VAGKAA)^w_VAGKAA)* cov[6]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[6]]/lambda_VAGKGS)^w_VAGKGS)* cov[6]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 7) {
  for (t in (schedule[7]+1) : schedule[8])  {  # day 541 to dose 8 on day 720
    df$prot_trip[t]<- exp(-(df$time[t-schedule[7]]/lambda_trip)^w_trip)* cov[7]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[7]]/lambda_quadr)^w_quadr)* cov[7]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[7]]/lambda_quint)^w_quint)* cov[7]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[7]]/lambda_sext)^w_sext)* cov[7]
    df$prot_other[t]<- exp(-(df$time[t-schedule[7]]/lambda_other)^w_other)* cov[7]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[7]]/lambda_VAGKAA)^w_VAGKAA)* cov[7]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[7]]/lambda_VAGKGS)^w_VAGKGS)* cov[7]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}

if (length(schedule) > 8) {
  for (t in (schedule[8]+1) : nrow(df))  {  # day 721 to end of simulation
    df$prot_trip[t]<- exp(-(df$time[t-schedule[8]]/lambda_trip)^w_trip)* cov[8]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[8]]/lambda_quadr)^w_quadr)* cov[8]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[8]]/lambda_quint)^w_quint)* cov[8]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[8]]/lambda_sext)^w_sext)* cov[8]
    df$prot_other[t]<- exp(-(df$time[t-schedule[8]]/lambda_other)^w_other)* cov[8]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[8]]/lambda_VAGKAA)^w_VAGKAA)* cov[8]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[8]]/lambda_VAGKGS)^w_VAGKGS)* cov[8]
  }
} else {
  for (t in (schedule[length(schedule)]+1) : nrow(df))  {
    df$prot_trip[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_trip)^w_trip)* cov[length(cov)]
    df$prot_quadr[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quadr)^w_quadr)* cov[length(cov)]
    df$prot_quint[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_quint)^w_quint)* cov[length(cov)]
    df$prot_sext[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_sext)^w_sext)* cov[length(cov)]
    df$prot_other[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_other)^w_other)* cov[length(cov)]
    df$prot_VAGKAA[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKAA)^w_VAGKAA)* cov[length(cov)]
    df$prot_VAGKGS[t]<- exp(-(df$time[t-schedule[length(schedule)]]/lambda_VAGKGS)^w_VAGKGS)* cov[length(cov)]
  }
}




##### SET UP EFFICACY CURVED (ADDITIONAL PMC SCHEDULE) #####

# set up df with the timesteps
if (sim_length > length(age_in_days_midpoint)) {
  df_additional_doses<- data.frame(time=1:sim_length)
}

if (sim_length < length(age_in_days_midpoint)) {
  df_additional_doses<- data.frame(time=1:length(age_in_days_midpoint))
}

# initalise empty columns
df_additional_doses$prot_trip<-NA
df_additional_doses$prot_quadr<-NA
df_additional_doses$prot_quint<-NA
df_additional_doses$prot_sext<-NA
df_additional_doses$prot_other<-NA
df_additional_doses$prot_VAGKAA<-NA
df_additional_doses$prot_VAGKGS<-NA


# construct weibull curves to model the efficacy of SP through time for
# each haplotype

if (length(schedule_additional_doses) > 0) {
  for (t in 1: schedule_additional_doses[1]) {  # day 0 to dose 1 on day 70
    df_additional_doses$prot_trip[t]<-0
    df_additional_doses$prot_quadr[t]<-0
    df_additional_doses$prot_quint[t]<-0
    df_additional_doses$prot_sext[t]<- 0
    df_additional_doses$prot_other[t]<-0
    df_additional_doses$prot_VAGKAA[t]<- 0
    df_additional_doses$prot_VAGKGS[t]<- 0
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 1) {
  for (t in (schedule_additional_doses[1]+1) : schedule_additional_doses[2])  {   # day 71 to dose 2 on day 98
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_trip)^w_trip) * cov_additional_doses[1]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_quadr)^w_quadr)* cov_additional_doses[1]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_quint)^w_quint)* cov_additional_doses[1]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_sext)^w_sext)* cov_additional_doses[1]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_other)^w_other)* cov_additional_doses[1]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[1]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[1]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[1]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 2) {
  for (t in (schedule_additional_doses[2]+1) : schedule_additional_doses[3])  {  #  day 99 to dose 3 on day 180
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_trip)^w_trip)* cov_additional_doses[2]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_quadr)^w_quadr)* cov_additional_doses[2]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_quint)^w_quint)* cov_additional_doses[2]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_sext)^w_sext)* cov_additional_doses[2]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_other)^w_other)* cov_additional_doses[2]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[2]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[2]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[2]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 3) {
  for (t in (schedule_additional_doses[3]+1) : schedule_additional_doses[4])  {  # day 181 to dose 4 on day 270
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_trip)^w_trip)* cov_additional_doses[3]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_quadr)^w_quadr)* cov_additional_doses[3]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_quint)^w_quint)* cov_additional_doses[3]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_sext)^w_sext)* cov_additional_doses[3]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_other)^w_other)* cov_additional_doses[3]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[3]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[3]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[3]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 4) {
  for (t in (schedule_additional_doses[4]+1) : schedule_additional_doses[5])  {  # day 271 to dose 5 on day 360
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_trip)^w_trip)* cov_additional_doses[4]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_quadr)^w_quadr)* cov_additional_doses[4]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_quint)^w_quint)* cov_additional_doses[4]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_sext)^w_sext)* cov_additional_doses[4]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_other)^w_other)* cov_additional_doses[4]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[4]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[4]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[4]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 5) {
  for (t in (schedule_additional_doses[5]+1) : schedule_additional_doses[6])  {  # day 361 to dose 6 on day 450
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_trip)^w_trip)* cov_additional_doses[5]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_quadr)^w_quadr)* cov_additional_doses[5]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_quint)^w_quint)* cov_additional_doses[5]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_sext)^w_sext)* cov_additional_doses[5]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_other)^w_other)* cov_additional_doses[5]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[5]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[5]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[5]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 6) {
  for (t in (schedule_additional_doses[6]+1) : schedule_additional_doses[7])  {  # day 451 to dose 7 on day 540
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_trip)^w_trip)* cov_additional_doses[6]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_quadr)^w_quadr)* cov_additional_doses[6]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_quint)^w_quint)* cov_additional_doses[6]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_sext)^w_sext)* cov_additional_doses[6]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_other)^w_other)* cov_additional_doses[6]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[6]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[6]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[6]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 7) {
  for (t in (schedule_additional_doses[7]+1) : schedule_additional_doses[8])  {  # day 541 to dose 8 on day 720
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_trip)^w_trip)* cov_additional_doses[7]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_quadr)^w_quadr)* cov_additional_doses[7]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_quint)^w_quint)* cov_additional_doses[7]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_sext)^w_sext)* cov_additional_doses[7]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_other)^w_other)* cov_additional_doses[7]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[7]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[7]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[7]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

if (length(schedule_additional_doses) > 8) {
  for (t in (schedule_additional_doses[8]+1) : nrow(df_additional_doses))  {  # day 721 to end of simulation
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_trip)^w_trip)* cov_additional_doses[8]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_quadr)^w_quadr)* cov_additional_doses[8]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_quint)^w_quint)* cov_additional_doses[8]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_sext)^w_sext)* cov_additional_doses[8]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_other)^w_other)* cov_additional_doses[8]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[8]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[8]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[8]
  }
} else {
  for (t in (schedule_additional_doses[length(schedule_additional_doses)]+1) : nrow(df_additional_doses))  {
    df_additional_doses$prot_trip[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_trip)^w_trip)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quadr[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quadr)^w_quadr)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_quint[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_quint)^w_quint)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_sext[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_sext)^w_sext)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_other[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_other)^w_other)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKAA[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKAA)^w_VAGKAA)* cov_additional_doses[length(cov_additional_doses)]
    df_additional_doses$prot_VAGKGS[t]<- exp(-(df_additional_doses$time[t-schedule_additional_doses[length(schedule_additional_doses)]]/lambda_VAGKGS)^w_VAGKGS)* cov_additional_doses[length(cov_additional_doses)]
  }
}

##### EFFICACY CURVE GRAPHS #####

colors <- c("triple" = "#00C094", "quadruple" = "#00B6EB", "quintuple" = "#FFA500",
            "VAGKAA"= "gold1", "sextuple" = "#F8766D" , "VAGKGS"= "#FB61D7", "other"= "#A58AFF")

# plot graphs for the probability of drug protection for each genotype through
# time during the SP dosing schedule given varying levels of coverage. The
# graph shapes will be the same but since coverage varies, the height differs

# graph1 <- ggplot(data=df)+ theme_bw() +
#   geom_line(aes(x=time, y=prot_trip,color="I_AKA_")) +
#   geom_line(aes(x=time, y=prot_quadr,color="I_GKA_")) +
#   geom_line(aes(x=time, y=prot_quint,color="I_GEA_")) +
#   geom_line(aes(x=time, y=prot_VAGKAA,color= "V_GKA_")) +
#   geom_line(aes(x=time, y=prot_sext, color="I_GEG_")) +
#   geom_line(aes(x=time, y=prot_VAGKGS, color="V_GKG_")) +
#   geom_line(aes(x=time, y=prot_other, color="other")) +
#   geom_vline(xintercept = schedule, color="blue", linetype="dashed") +
#   labs(color="Resistance genotype") +
#   ylab("Probability of drug protection")+xlab("Age(days)")
# 
# graph2 <- ggplot(data=df_additional_doses)+ theme_bw() +
#   geom_line(aes(x=time, y=prot_trip,color="I_AKA_")) +
#   geom_line(aes(x=time, y=prot_quadr,color="I_GKA_")) +
#   geom_line(aes(x=time, y=prot_quint,color="I_GEA_")) +
#   geom_line(aes(x=time, y=prot_VAGKAA,color= "V_GKA_")) +
#   geom_line(aes(x=time, y=prot_sext, color="I_GEG_")) +
#   geom_line(aes(x=time, y=prot_VAGKGS, color="V_GKG_")) +
#   geom_line(aes(x=time, y=prot_other, color="other")) +
#   geom_vline(xintercept = schedule_additional_doses, color="blue", linetype="dashed") +
#   labs(color="Resistance genotype") +
#   ylab("Probability of drug protection")+xlab("Age(days)")

##### CALCULATE OVERALL EFFICACY FOR EACH ADMIN-1 AREA #####

# loop over each admin-1 area and calculate weighted mean for overall efficacy
# based on proportion of each haplotype present

admin1_with_haplotype_data <- c()

for (i in 1:length(admin1)){
  proportions <- haplotype_proportions %>% filter(iso_code == country, NAME_2 == admin1[i])
  
  if (dim(proportions)[1] != 0) {
    
    admin1_with_haplotype_data <- c(admin1_with_haplotype_data, admin1[i])
    
    df[paste0("prot_overall_", admin1[i])] <- as.double(proportions$I_AKA_)*df$prot_trip +
      as.double(proportions$I_GKA_)*df$prot_quadr +
      as.double(proportions$I_GEA_)*df$prot_quint +
      as.double(proportions$I_GEG_)*df$prot_sext +
      as.double(proportions$V_GKA_)*df$prot_VAGKAA +
      as.double(proportions$V_GKG_)*df$prot_VAGKGS +
      as.double(proportions$other)*df$prot_other
    
    df_additional_doses[paste0("prot_overall_", admin1[i])]<- as.double(proportions$I_AKA_)*df_additional_doses$prot_trip +
      as.double(proportions$I_GKA_)*df_additional_doses$prot_quadr +
      as.double(proportions$I_GEA_)*df_additional_doses$prot_quint +
      as.double(proportions$I_GEG_)*df_additional_doses$prot_sext +
      as.double(proportions$V_GKA_)*df_additional_doses$prot_VAGKAA +
      as.double(proportions$V_GKG_)*df_additional_doses$prot_VAGKGS +
      as.double(proportions$other)*df_additional_doses$prot_other
  }
  
}

admin1 <- admin1_with_haplotype_data

##### CALCULATE INCIDENCE FOR WHOLE COUNTRY #####

# NOTE: population weights are calculated out of the admin-1 units that have
# haplotype data so admin-1 units that have no data take the average
# 
# incidence_ppy_df_whole_country_weighted <- (incidence_ppy_df %>% filter(area %in% admin1))
population_weights <- c()

pop_size_across_admin1 <- sum((full_data$population %>% filter(year==2023, name_2 %in% admin1))$pop)

for (i in 1:length(admin1)){

  pop <- sum((full_data$population %>% filter(year == 2023, name_2 == admin1[i]))$pop)

  pop_weight <- pop / pop_size_across_admin1
  population_weights <- c(population_weights, pop_weight)
}

population_weights <- rep(population_weights, each=4*length(age_in_days_midpoint))
# incidence_ppy_df_whole_country_weighted$value <- incidence_ppy_df_whole_country_weighted$value * population_weights
# 
# incidence_ppy_df_whole_country <- data.frame(age_in_days_midpoint)
# incidence_ppy_df_whole_country_age_sum_clin <- incidence_ppy_df_whole_country_age_sum_sev <- incidence_ppy_df_whole_country_age_sum_tot <- incidence_ppy_df_whole_country_age_sum_asym <- c()

SP_protection_additional_doses_whole_country <- c()
SP_protection_whole_country <- c()

for (i in 1:length(age_in_days_midpoint)) {
  
  # incidence 
  # incidence_ppy_df_whole_country_age_sum_clin <- c(incidence_ppy_df_whole_country_age_sum_clin, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "clinical", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  # incidence_ppy_df_whole_country_age_sum_sev <- c(incidence_ppy_df_whole_country_age_sum_sev, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "severe", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  # incidence_ppy_df_whole_country_age_sum_tot <- c(incidence_ppy_df_whole_country_age_sum_tot, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "total", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  # incidence_ppy_df_whole_country_age_sum_asym <- c(incidence_ppy_df_whole_country_age_sum_asym, sum((incidence_ppy_df_whole_country_weighted %>% filter(infection_class == "asymptomatic", age_in_days_midpoint == age_in_days_midpoint[i]))$value))
  # 
  # SP protection (sum of weighted mean (protection in each admin-1 area with population size weights) at each time point)
  
  SP_protection_whole_country <- c(SP_protection_whole_country, sum(unique(population_weights) * df[i,(dim(df)[2] - length(admin1) + 1):dim(df)[2]])) # only admin-1 area protections
  SP_protection_additional_doses_whole_country <- c(SP_protection_additional_doses_whole_country, sum(unique(population_weights) * df_additional_doses[i,(dim(df_additional_doses)[2] - length(admin1) + 1):dim(df_additional_doses)[2]])) # only admin-1 area protections
  
}

# incidence_ppy_df_whole_country$clinical <- incidence_ppy_df_whole_country_age_sum_clin
# incidence_ppy_df_whole_country$severe <- incidence_ppy_df_whole_country_age_sum_sev
# incidence_ppy_df_whole_country$total <- incidence_ppy_df_whole_country_age_sum_tot
# incidence_ppy_df_whole_country$asymptomatic <- incidence_ppy_df_whole_country_age_sum_asym


##### CALCULATE PMC IMPACT FOR EACH ADMIN-1 AREA (both PMC schedules) #####

PMC_impact_ppy_additional_doses <- data.frame()
PMC_impact_ppy <- data.frame()

for (i in 1:length(admin1)){
  
  # PMC impact on incidence (ppy)
  new_PMC_impact_ppy_additional_doses_df <- (incidence_ppy_df %>% filter(area == admin1[i]))
  new_PMC_impact_ppy_df <- (incidence_ppy_df %>% filter(area == admin1[i]))
  
  PMC_impact_ppy_additional_doses_cases <- as.double(new_PMC_impact_ppy_additional_doses_df$value) * rep((1-df_additional_doses[paste0("prot_overall_", admin1[i])][min(age_min):max(age_max),]),4)
  PMC_impact_ppy_cases <- as.double(new_PMC_impact_ppy_df$value) * rep((1-df[paste0("prot_overall_", admin1[i])][min(age_min):max(age_max),]),4)
  
  new_PMC_impact_ppy_additional_doses_df$value <- PMC_impact_ppy_additional_doses_cases
  new_PMC_impact_ppy_df$value <- PMC_impact_ppy_cases
  
  PMC_impact_ppy_additional_doses <- rbind(PMC_impact_ppy_additional_doses, new_PMC_impact_ppy_additional_doses_df)
  PMC_impact_ppy <- rbind(PMC_impact_ppy, new_PMC_impact_ppy_df)

  }

##### PMC IMPACT FOR WHOLE COUNTRY #####


PMC_impact_ppy_additional_doses_whole_country <- data.frame(age_in_days_midpoint)
PMC_impact_ppy_whole_country <- data.frame(age_in_days_midpoint)

PMC_impact_ppy_additional_doses_whole_country$clinical <- incidence_ppy_df_whole_country$clinical * (1-SP_protection_additional_doses_whole_country)
PMC_impact_ppy_additional_doses_whole_country$severe <- incidence_ppy_df_whole_country$severe * (1-SP_protection_additional_doses_whole_country)
PMC_impact_ppy_additional_doses_whole_country$total <- incidence_ppy_df_whole_country$total * (1-SP_protection_additional_doses_whole_country)
PMC_impact_ppy_additional_doses_whole_country$asymptomatic <- incidence_ppy_df_whole_country$asymptomatic * (1-SP_protection_additional_doses_whole_country)

PMC_impact_ppy_whole_country$clinical <- incidence_ppy_df_whole_country$clinical * (1-SP_protection_whole_country)
PMC_impact_ppy_whole_country$severe <- incidence_ppy_df_whole_country$severe * (1-SP_protection_whole_country)
PMC_impact_ppy_whole_country$total <- incidence_ppy_df_whole_country$total * (1-SP_protection_whole_country)
PMC_impact_ppy_whole_country$asymptomatic <- incidence_ppy_df_whole_country$asymptomatic * (1-SP_protection_whole_country)

colours <- c( "#E16E21","#8FC7C7")

# PMC <- 1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value
# no_PMC <- 1000*(incidence_ppy_df %>% filter(area == chosen_area, infection_class == "clinical"))$value
# 
# graph1<- ggplot() +
#   geom_line(aes(x=1:913, y=1000*(incidence_ppy_df %>% filter(area == chosen_area, infection_class == "clinical"))$value, colour = "No PMC")) +
#   geom_line(aes(x=1:913, y=1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value, colour = "PMC"))+
#   scale_color_manual("Incidence", breaks=c("No PMC", "PMC"),
#                      values = colours) + 
#   geom_ribbon(aes(x=1:913, ymin=1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value,
#                   ymax=1000*(incidence_ppy_df %>% filter(area == chosen_area, infection_class == "clinical"))$value, fill="Cases averted"), alpha=0.3) +
#   
#   geom_segment(aes(x = schedule[schedule < max(age_max)], y=max(1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value + 500),
#                    xend = schedule[schedule < max(age_max)], yend=max(1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value) + 100)
#                , arrow = arrow(length= unit(0.1, "inches"))) +
#   geom_text(aes(x=schedule[schedule < max(age_max)], label=paste0((cov[1:length(schedule[schedule < max(age_max)])] * 100), "%"), y=max(1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value) + 700),
#             colour="black", angle=90, size = 2.25, nudge_y=150) +
#   labs(x="Age (days)", y="New clinical infections per 1000 children",
#        colour = "Incidence", fill="Calculation") +
#   ylim(0, max(1000*(PMC_impact_ppy %>% filter(area == chosen_area, infection_class == "clinical"))$value) + 1200) +
#   scale_x_continuous(breaks = seq(0, 913, by = 90)) + 
#   scale_fill_manual("Calculation", values = colours[2]) +
#   theme(#panel.grid.minor = element_line(colour = "gray", size = 0.2),
#     panel.grid.major = element_line(colour = "gray", size = 0.2),
#     panel.background = element_rect(fill='transparent'), #transparent panel bg
#     plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#     legend.background = element_rect(fill='transparent'), #transparent legend bg
#   )
# 
# 
# graph2<-ggplot() +
#   geom_line(aes(x=1:913, y=1000*(incidence_ppy_df %>% filter(area == chosen_area, infection_class == "clinical"))$value, colour = "No PMC")) +
#   geom_line(aes(x=1:913, y=1000*(PMC_impact_ppy_additional_doses %>% filter(area == chosen_area, infection_class == "clinical"))$value, colour = "PMC"))+
#   scale_color_manual("Incidence", breaks=c("No PMC", "PMC"),
#                      values = colours) + 
#   geom_ribbon(aes(x=1:913, ymin=1000*(PMC_impact_ppy_additional_doses %>% filter(area == chosen_area, infection_class == "clinical"))$value,
#                   ymax=1000*(incidence_ppy_df %>% filter(area == chosen_area, infection_class == "clinical"))$value, fill="Cases averted"), alpha=0.3) +
#   
#   geom_segment(aes(x = schedule_additional_doses[schedule_additional_doses < max(age_max)], y=max(1000*(PMC_impact_ppy_additional_doses %>% filter(area == chosen_area, infection_class == "clinical"))$value + 500),
#                    xend = schedule_additional_doses[schedule_additional_doses < max(age_max)], yend=max(1000*(PMC_impact_ppy_additional_doses %>% filter(area == chosen_area, infection_class == "clinical"))$value) + 100)
#                , arrow = arrow(length= unit(0.1, "inches"))) +
#   geom_text(aes(x=schedule_additional_doses[schedule_additional_doses < max(age_max)], label=paste0((cov_additional_doses[1:length(schedule_additional_doses[schedule_additional_doses < max(age_max)])] * 100), "%"), y=max(1000*(PMC_impact_ppy_additional_doses %>% filter(area == chosen_area, infection_class == "clinical"))$value) + 700),
#             colour="black", angle=90, size = 2.25, nudge_y=150) +
#   labs(x="Age (days)", y="New clinical infections per 1000 children",
#        colour = "Incidence", fill="Calculation") +
#   ylim(0, max(1000*(PMC_impact_ppy_additional_doses %>% filter(area == chosen_area, infection_class == "clinical"))$value) + 1200) +
#   scale_x_continuous(breaks = seq(0, 913, by = 90)) + 
#   scale_fill_manual("Calculation", values = colours[2]) +
#   theme(#panel.grid.minor = element_line(colour = "gray", size = 0.2),
#     panel.grid.major = element_line(colour = "gray", linewidth = 0.2),
#     panel.background = element_rect(fill='transparent'), #transparent panel bg
#     plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#     legend.background = element_rect(fill='transparent'), #transparent legend bg
#   )
# 
# merged_graphs <- ggarrange(graph1,graph2)
# merged_graphs<- annotate_figure(merged_graphs, top = text_grob(paste0("Average clinical incidence in ", chosen_area, ", ", country, " during ",
#                                                                      length(schedule), " and ", length(schedule_additional_doses), " dose PMC schedules")))
#  
# print(merged_graphs) 

# graph3 <- ggplot() +
#   geom_line(aes(x=unique(incidence_ppy_df_whole_country$age_in_days_midpoint), y=1000*(incidence_ppy_df_whole_country$clinical), colour = "No PMC")) + 
#   geom_line(aes(x=unique(PMC_impact_ppy_whole_country$age_in_days_midpoint), y=1000*as.double(PMC_impact_ppy_whole_country$clinical), colour = "PMC"))+
#   scale_color_manual("Incidence", breaks=c("No PMC", "PMC"),
#                                            values = colours) +
#   geom_ribbon(aes(x=unique(incidence_ppy_df_whole_country$age_in_days_midpoint), ymin=1000*as.double(PMC_impact_ppy_whole_country$clinical),
#                   ymax=1000*(incidence_ppy_df_whole_country$clinical), fill="Cases averted"), alpha=0.3) + 
#   
#   geom_segment(aes(x = schedule[schedule < max(age_max)], y=max(1000*(PMC_impact_ppy_whole_country$clinical)) + 500,
#                    xend = schedule[schedule < max(age_max)], yend=max(1000*(PMC_impact_ppy_whole_country$clinical)) + 100),
#                arrow = arrow(length= unit(0.1, "inches"))) + 
#   geom_text(aes(x=schedule[schedule < max(age_max)], label=paste0((cov[1:length(schedule[schedule < max(age_max)])] * 100), "%"), y=max(1000*(PMC_impact_ppy_whole_country$clinical) + 700)),
#             colour="black", angle=90, size = 2.25, nudge_y=150) +
#   labs(x="Age (days)", y="New clinical infections per 1000 children",
#        colour = "Incidence", fill="Calculation") +
#   ylim(0, max(1000*(PMC_impact_ppy_whole_country$clinical) + 1200)) + 
#   scale_x_continuous(breaks = seq(0, 913, by = 90)) + 
#   scale_fill_manual("Calculation", values = colours[2]) +
#   theme(#panel.grid.minor = element_line(colour = "gray", size = 0.2),
#     panel.grid.major = element_line(colour = "gray", linewidth = 0.2),
#     panel.background = element_rect(fill='transparent'), #transparent panel bg
#     plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#     legend.background = element_rect(fill='transparent'), #transparent legend bg
#   )
# 
# 
# graph4 <- ggplot() +
#   geom_line(aes(x=unique(incidence_ppy_df_whole_country$age_in_days_midpoint), y=1000*(incidence_ppy_df_whole_country$clinical), colour = "No PMC")) + 
#   geom_line(aes(x=unique(PMC_impact_ppy_additional_doses_whole_country$age_in_days_midpoint), y=1000*as.double(PMC_impact_ppy_additional_doses_whole_country$clinical), colour = "PMC"))+
#   
#   geom_ribbon(aes(x=unique(incidence_ppy_df_whole_country$age_in_days_midpoint), ymin=1000*as.double(PMC_impact_ppy_additional_doses_whole_country$clinical),
#                   ymax=1000*(incidence_ppy_df_whole_country$clinical), fill="Cases averted"), alpha=0.3) + 
#   
#   geom_segment(aes(x = schedule_additional_doses[schedule_additional_doses < max(age_max)], y=max(1000*(PMC_impact_ppy_additional_doses_whole_country$clinical)) + 500,
#                    xend = schedule_additional_doses[schedule_additional_doses < max(age_max)], yend=max(1000*(PMC_impact_ppy_additional_doses_whole_country$clinical)) + 100)
#                , arrow = arrow(length= unit(0.1, "inches"))) + 
#   geom_text(aes(x=schedule_additional_doses[schedule_additional_doses < max(age_max)], label=paste0((cov_additional_doses[1:length(schedule_additional_doses[schedule_additional_doses < max(age_max)])] * 100), "%"), y=max(1000*(PMC_impact_ppy_additional_doses_whole_country$clinical) + 700)),
#             colour="black", angle=90, size = 2.25, nudge_y=150) +
#   labs(x="Age (days)", y="New clinical infections per 1000 children",
#        colour = "Incidence", fill="Calculation") +
#   ylim(0, max(1000*(PMC_impact_ppy_additional_doses_whole_country$clinical) + 1200)) + 
#   scale_x_continuous(breaks = seq(0, 913, by = 90)) + 
#   scale_fill_manual("Calculation", values = colours[2]) +
#   theme(#panel.grid.minor = element_line(colour = "gray", size = 0.2),
#     panel.grid.major = element_line(colour = "gray", linewidth = 0.2),
#     panel.background = element_rect(fill='transparent'), #transparent panel bg
#     plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#     legend.background = element_rect(fill='transparent'), #transparent legend bg
#   )

# merged_graphs <- ggarrange(graph3,graph4)
# merged_graphs<- annotate_figure(merged_graphs, top = text_grob(paste0("Average clinical incidence in",  country, " during ",
#                                                                      length(schedule), " and ", length(schedule_additional_doses), " dose PMC schedules")))
#  
# print(merged_graphs) 