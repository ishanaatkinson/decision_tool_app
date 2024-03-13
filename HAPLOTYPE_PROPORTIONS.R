# library(pacman)
# p_load(malariasimulation, foreSIGHT, malariaEquilibrium, tidyr, dplyr, ggplot2,
#        reshape2, ggpubr, gridExtra, readxl, stringi, scene, XML, maps, readr,
#        here, sf)

# set your working directory to where the haplotype data is found
#setwd(paste0("C:/Users/IshanaAtkinson/OneDrive - London School of Hygiene and Tropical Medicine/Documents/malariasimulation work"))

# can add more countries when I clean excel data further
#iso_codes <- c("Benin", "Mozambique", "Cameroon", "Cote d'Ivoire", "Zambia", "Ghana")
iso_codes <- c("BEN", "MOZ", "CMR", "CIV", "ZMB", "GHA")

# set this to the name of your haplotype data file
haplotype_data <- read_xlsx(paste0(getwd(), "/", country, "/haplotypes_combined_forIshana2.xlsx"))

# converts admin1 names into clean format that R can understand 
haplotype_data$NAME_2 <- stri_trans_general(str=gsub("-", "_", haplotype_data$NAME_1), id = "Latin-ASCII")
haplotype_data$NAME_2 <- stri_trans_general(str=gsub(" ", "_", haplotype_data$NAME_2), id = "Latin-ASCII")

haplotype_proportions <- data.frame()

I_AKA_ <- I_GKA_ <- I_GEA_ <- I_GEG_ <- V_GKA_ <- V_GKG_ <- c()

# loops over each country and admin1 area and calculates proportions for areas 
# which have 6 codon data and appends to data frame

for (i in 1:length(iso_codes)) {
  
  country_haplotype_data <- haplotype_data %>% filter(iso_code == iso_codes[i])
  admin1 <- unique(country_haplotype_data$NAME_2)
  
  # calculate 6 codon data normally but include 5 codon data in these countries 
  if ((iso_codes[i] %in% c("MOZ", "ZMB", "DRC"))) {
    for (j in 1:length(admin1)) {
      
      # area haplotype data
      hap <- country_haplotype_data %>% filter(NAME_2 == admin1[j])
      
      # include 5 codon genotypes for these countries
      denom <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS, hap$SAKAA, hap$AAKAA, hap$FAKAA, hap$SAKAS, hap$AAKAS, hap$FAKAS, hap$YAKAS, hap$SAKAA_AAKAA_FAKAA, hap$AAKAS_FAKAS_YAKAS_SAKAS,
                                hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS, hap$SGKAA, hap$AGKAA, hap$FGKAA, hap$SGKAS, hap$AGKAS, hap$FGKAS, hap$CGKAA_SGKAA_AGKAA, hap$AGKAS_SGKAS,
                                hap$ISGEAA, hap$IAGEAA, hap$ISGEAS, hap$SGEAA, hap$AGEAA, hap$SGEAA_AGEAA,
                                hap$ISGEGA, hap$IAGEGA, hap$AGEGA, hap$SGEGS, hap$SGEGA,
                                hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS,
                                hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA,
                                hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA, hap$SGKGS, hap$AGKGS, hap$AAEGA, hap$SAEGA, hap$AGKGA, hap$SAEAA, hap$SAKGA, hap$SAKGS, hap$AAKGS, hap$SGKGA, hap$AAEAS, hap$FGKGA, hap$FGKGS
      )), na.rm=TRUE)
      
      I_AKA_ <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS,hap$SAKAA, hap$AAKAA, hap$FAKAA, hap$SAKAS, hap$AAKAS, hap$FAKAS, hap$YAKAS, hap$SAKAA_AAKAA_FAKAA, hap$AAKAS_FAKAS_YAKAS_SAKAS)), na.rm=TRUE) / denom
      I_GKA_ <- sum(as.integer(c(hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS,hap$SGKAA, hap$AGKAA, hap$FGKAA, hap$SGKAS, hap$AGKAS, hap$FGKAS, hap$CGKAA_SGKAA_AGKAA, hap$AGKAS_SGKAS)), na.rm=TRUE) / denom
      I_GEA_ <- sum(as.integer(c(hap$ISGEAA, hap$IAGEAA, hap$ISGEAS, hap$SGEAA, hap$AGEAA, hap$SGEAA_AGEAA)), na.rm=TRUE) / denom
      I_GEG_ <- sum(as.integer(c(hap$ISGEGA, hap$IAGEGA, hap$AGEGA, hap$SGEGS, hap$SGEGA)), na.rm=TRUE) / denom
      V_GKA_ <- sum(as.integer(c(hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS)), na.rm=TRUE) / denom
      V_GKG_ <- sum(as.integer(c(hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA)), na.rm=TRUE) / denom
      other <- sum(as.integer(c(hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA, hap$SGKGS, hap$AGKGS, hap$AAEGA, hap$SAEGA, hap$AGKGA, hap$SAEAA, hap$SAKGA, hap$SAKGS, hap$AAKGS, hap$SGKGA, hap$AAEAS, hap$FGKGA, hap$FGKGS)), na.rm=TRUE) / denom
      
      props_5sf <- round(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other), digits=5)
      sum_of_prop <- sum(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other))
      
      # check proportions add to 1 
      if (is.nan(sum_of_prop) == TRUE) {
        print(paste0("Summation of proportions for ", admin1[j], " in ", iso_codes[i], " cannot be calculated due to lack of 6-codon data"))
      } else {
        print(paste0("Summation of proportions for ", admin1[j], " in ", iso_codes[i], " is: ", sum_of_prop))
        
        haplotype_proportions <- rbind(haplotype_proportions, c(iso_codes[i], admin1[j], mean(hap$longitude, na.rm = TRUE), mean(hap$latitude, na.rm = TRUE),  props_5sf, sum_of_prop))
      }
      
    }
  } else {
  
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
    
    I_AKA_ <- sum(as.integer(c(hap$IAAKAA, hap$IFAKAA, hap$ISAKAA, hap$CAKAA, hap$HAKAA, hap$IYAKAS, hap$IFAKAS, hap$ISAKAS, hap$IAAKAS, hap$SAKAA_AAKAA_FAKAA, hap$AAKAS_FAKAS_YAKAS_SAKAS)), na.rm=TRUE) / denom
    I_GKA_ <- sum(as.integer(c(hap$ISGKAA, hap$IAGKAA, hap$IAGKAS, hap$ISGKAS, hap$IFGKAA, hap$IFGKAS)), na.rm=TRUE) / denom
    I_GEA_ <- sum(as.integer(c(hap$ISGEAA, hap$IAGEAA, hap$ISGEAS)), na.rm=TRUE) / denom
    I_GEG_ <- sum(as.integer(c(hap$ISGEGA, hap$IAGEGA)), na.rm=TRUE) / denom
    V_GKA_ <- sum(as.integer(c(hap$VSGKAA, hap$VAGKAA, hap$VAGKAS, hap$VSGKAS)), na.rm=TRUE) / denom
    V_GKG_ <- sum(as.integer(c(hap$VSGKGA, hap$VAGKGS, hap$VSGKGS, hap$VAGKGA)), na.rm=TRUE) / denom
    other <- sum(as.integer(c(hap$VAAKAA,  hap$VSAKAA, hap$VSGEAA, hap$ISGKGA, hap$VSGEGA, hap$VAGEAA, hap$VAGEGA, hap$IAGKGS, hap$ISGKGS, hap$IAAKGS, hap$other, hap$IAGKGA, hap$IAAKGA, hap$ISAEAA)), na.rm=TRUE) / denom
    
    props_5sf <- round(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other), digits=5)
    sum_of_prop <- sum(c(I_AKA_, I_GKA_, I_GEA_, I_GEG_, V_GKA_, V_GKG_, other))
    
    # check proportions add to 1 
    if (is.nan(sum_of_prop) == TRUE) {
      print(paste0("Summation of proportions for ", admin1[j], " in ", iso_codes[i], " cannot be calculated due to lack of 6-codon data"))
    } else {
      print(paste0("Summation of proportions for ", admin1[j], " in ", iso_codes[i], " is: ", sum_of_prop))
      
      haplotype_proportions <- rbind(haplotype_proportions, c(iso_codes[i], admin1[j], mean(hap$longitude, na.rm = TRUE), mean(hap$latitude, na.rm = TRUE),  props_5sf, sum_of_prop))
    }
    
  }
  }
}

# change column names 
colnames(haplotype_proportions) <- c("iso_code", "NAME_2", "longitude", "latitude", "I_AKA_", "I_GKA_", "I_GEA_", "I_GEG_", "V_GKA_", "V_GKG_", "other", "Sum of proportions" )


#haplotype_proportions



