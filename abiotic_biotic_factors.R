##################### Adding Abiotic and Biotic Factors ########################################
library(dplyr)
library(rgbif)
library(dplyr)
library(countrycode)
library(geosphere)

plant_features = read.csv("~/plant_features_database_no_Sphagneticola.csv")
# create a dataframe for collecting environmental and ecological factors for each location and date
locations = unique(plant_features$Location)
countries = NULL
for(loc in locations) {
  country = first(plant_features$Country[which(plant_features$Location == loc)]) #extract first value of a variable for which a location is TRUE
  countries = append(countries,country)
}
articles = NULL
for(loc in locations) {
  article = first(plant_features$source_article[which(plant_features$Location == loc)])
  articles = append(articles,article)
}
species = NULL
for(loc in locations) {
  speci = first(plant_features$Species[which(plant_features$Location == loc)])
  species = append(species,speci)
}

wild.not_wild = NULL
for(loc in locations) {
  w.nw = first(plant_features$wild.not_wild[which(plant_features$Location == loc)])
  wild.not_wild = append(wild.not_wild,w.nw)
}

Habitat = NULL
for(loc in locations) {
  habi = first(plant_features$Habitat[which(plant_features$Location == loc)])
  Habitat = append(Habitat,habi)
}

Sampling_date = NULL
for(loc in locations) {
  sample = first(plant_features$sampling_date[which(plant_features$Location == loc)])
  Sampling_date = append(Sampling_date,sample)
}

factors_df = as.data.frame(cbind(articles,species,countries,locations,wild.not_wild,Habitat,Sampling_date))
write.csv(factors_df,"C:\\Users\\Jamila\\Documents\\factors_df.csv",row.names=FALSE)

###################################################################### Biotic Factors ###########################


########################### test 
india_data = occ_data(countr="IN",year=2006) #finds occurance data for india in 2006
india_data_2 = india_data$data
india_Latitude = india_data_2 %>% 
  dplyr::select(decimalLatitude) #extracts decimal latitude data from the tibble
india_Latitude_vec = pull(india_Latitude) #converts tibble to vector or dataframe

mad_data = occ_data(country="MG",year=2006,hasCoordinate = TRUE) #hasCoordinate means only data with coordinates are obtained

#mad_data_2 =  mad_data$data
#mad_Latitude = mad_data_2 %>% 
#  dplyr::select(decimalLatitude)
#mad_Latitude_vec = pull(mad_Latitude)
#mad_data_3 = as.data.frame(mad_data_2)
#unique(mad_data_3)

mad_data = occ_data(country="MG",year=c(2000:2010),hasCoordinate = TRUE) #hasCoordinate means only data with coordinates are obtained
mad_data_3 = NULL
for(i in 1:length(mad_data)) {
  mad_data_2 = mad_data[i]
  year = names(mad_data)[[i]]
  mad_data_2 = mad_data_2[[year]]$data
  mad_data_2 = as.data.frame(mad_data_2)
  mad_subset = subset(mad_data_2, select=c(scientificName,decimalLatitude,decimalLongitude,datasetKey,taxonRank,species,year))
  mad_data_3 = rbind(mad_data_3,mad_subset)
  
}

years = c(2003:2005)
mad_data_y = mad_data_3[which(mad_data_3$year %in% years),]

# test for calculating lat/long distances
distHaversine(c(-18,48),c(-19,48.65))
# gives rough distance between sets of lat/long points in m
 
###############################

########################### code run ##########################################################################

######### code for plant species richness 
plantae_key = name_backbone("Plantae")
plantae_key = plantae_key %>%
  dplyr::select(kingdomKey)
plantae_key = plantae_key[[1]]

taxon = plantae_key
factors_df = read.csv("~/factors_df.csv")
dataframe = factors_df[2:nrow(factors_df),]
overall_sample = c(1970:2020)
  
  
  richness_table = as.data.frame(matrix(NA, ncol = 8, nrow = nrow(dataframe)))
  data_p = 1

  while(data_p < (nrow(dataframe)+1)) {
    richness_table[data_p,1] = dataframe$articles[[data_p]]
    country = dataframe$countries[[data_p]]
    richness_table[data_p,2] = country
    if(country == "Hawaii") {
      country = "United States of America"
    }
    con_code = countrycode(country,origin ="country.name",destination = "iso2c") #convert country name into iso2c code - accepted code by gbif for identifying countries
    richness_table[data_p,3] = dataframe$locations[[data_p]]
    lat = as.numeric(dataframe$Latitude[[data_p]])
    long = as.numeric(dataframe$Longitude[[data_p]])
    richness_table[data_p,4] = lat
    richness_table[data_p,5] = long
    
    # select occurance data +/- 10 years from the sampling data, or from 1980-2020 if there is no sampling date given
    sample_y = as.numeric(dataframe$Sample_year[[data_p]])
    if(is.na(sample_y) == TRUE) {
      sampling = c(1980:2020)
    } else {
        sampling_low = sample_y-10
        sampling_high = sample_y+10
        if(sampling_high > 2020) {
          sampling_highu = 2020
        } else {
          sampling_highu = sampling_high
        }
        sampling = c(sampling_low,sampling_highu)
    }
    
    if(is.na(lat) == FALSE & is.na(long) == FALSE) {
      
      if(con_code %in% richness_table[,8] == FALSE) {  # if a dataset for a country has not been obtained
        
        country_data = occ_data(country=con_code,year=overall_sample,hasCoordinate = TRUE,kingdomKey=taxon) # obtain dataset of a given kingdom and country, only data with coordinates
        format_country = NULL
        
        for(i in 1:length(country_data)) {
          country_data2 = country_data[i]
          year = names(country_data)[[i]]
          country_data2 = country_data2[[year]]$data
          if(is.null(country_data2) == FALSE & is.null(country_data2$species) == FALSE) {
            country_data2 = as.data.frame(country_data2) # creates a dataframe of downloaded data
            country_subset = subset(country_data2, select=c(scientificName,decimalLatitude,decimalLongitude,datasetKey,taxonRank,species,year))
            format_country = rbind(format_country,country_subset)
          }
          
        }
        country_only_species = format_country[is.na(format_country$species) == FALSE,] #removes non-species records
        assign(con_code,country_only_species) # creates a dataset named a country's iso2c of all its cleaned gbif data
      } else {
        country_only_species = eval(as.name(con_code)) # if a dataset for a country has already been downloaded, this code pulls up this dataset
      }
      
      country_year_range = country_only_species[which(country_only_species$year %in% sampling),]
      
      
      # this function contains 2 methods of data collection 
      # 1) select records closest to the emission profile collection site, data from multiple databases
      # 2) select databases with the average distance of samples taken in said database is closest to the emission profile collection site.
      
      distances_df = as.data.frame(matrix(NA,ncol=6,nrow=nrow(country_year_range)))
      distances_df[,1] = rep(lat,length=nrow(distances_df)) # creates repeats of the emission latitude & longitude data data - not much function
      distances_df[,2] = rep(long,length=nrow(distances_df))
      names(distances_df) = c("emission_lat","emission_long","gbif_lat","gbif_long","distance_from_emission_to_gbif_site","datasetKey")
      
      for(dat in 1:nrow(country_year_range)) {
        
        lat_gbif = country_year_range$decimalLatitude[[dat]]
        long_gbif = country_year_range$decimalLongitude[[dat]]
        distances_df[dat,3] = lat_gbif
        distances_df[dat,4] = long_gbif
        
        distance_m = distHaversine(c(long,lat),c(long_gbif,lat_gbif)) # calculates the distance between lat/long sets of coordinates using the haversine equation
        distance_km = distance_m/1000 # distance_m is in m, code converts it to km
        distances_df[dat,5] = distance_km
        
        database = country_year_range$datasetKey[[dat]] # adds the database a coordinate comes from
        distances_df[dat,6] = database
        
      }
      
      # function 1) calculate species richness from data from all datasets <=75km from the emission profile collection site.
      data_to_use_1 = which(distances_df[,5] <= 75)
      species_1 = country_year_range$species[data_to_use_1]
      unique_species1 = unique(species_1) # identifies unique species <=75km from the emission profile collection site
      species_richness_1 = length(unique_species1) # finds no. of unique species aka species richness
      richness_table[data_p,6] = species_richness_1
      
      
      #function 2) species richness from datasets with a mean sampling distance <=75km from the emission profile collection site.
      uniq_datasets = unique(distances_df[,6])
      means_dataset = NULL
      for(d in uniq_datasets) {
        distance_dataset = distances_df[which(distances_df[,6] == d),5] # extracts data from a single dataset
        mean_dist = mean(distance_dataset) # calculates the mean distance of that dataset's data from the emission profile collection site
        means_dataset = append(means_dataset,mean_dist) 
      }
      data_to_use_2 = uniq_datasets[which(means_dataset <= 75)]
      species_2 = country_year_range$species[which(country_year_range$datasetKey %in% data_to_use_2)] # creates a vector of species from datasets where the mean distance from the emission profile site <=75km
      unique_species2 = unique(species_2)
      species_richness_2 = length(unique_species2) 
      richness_table[data_p,7] = species_richness_2
      
      
    }
    richness_table[data_p,8] = con_code # adds the country iso2c so if a country's data has already been downloaded, it can be identified in following iterations
    data_p = data_p + 1
  }


  
# do species richness by country   

richness_table = richness_table[,1:8]
richness_table = cbind(richness_table,rep(NA,length=nrow(richness_table)))
  
con = 1
while(con < nrow(richness_table)+1) {
  
  con_code2 = richness_table[con,8]
  
  
    sample_y = as.numeric(dataframe$Sample_year[[con]])
    if(is.na(sample_y) == TRUE) {
      sampling = c(1980:2020)
    } else {
      sampling_low = sample_y-10
      sampling_high = sample_y+10
      if(sampling_high > 2020) {
        sampling_highu = 2020
      } else {
        sampling_highu = sampling_high
      }
      sampling = c(sampling_low,sampling_highu)
    }
    if(is.na(con_code2) == FALSE) {
      
      if(is.na(richness_table[con,4]) == TRUE) {
        country_data_2 = occ_data(country=con_code2,hasCoordinate = TRUE,kingdomKey=taxon,year=overall_sample)
        format_country2 = NULL
        
        for(i in 1:length(country_data_2)) {
          country_data2 = country_data_2[i]
          year = names(country_data_2)[[i]]
          country_data2 = country_data2[[year]]$data
          if(is.null(country_data2) == FALSE & is.null(country_data2$species) == FALSE) {
            country_data2 = as.data.frame(country_data2) # creates a dataframe of downloaded data
            country_subset = subset(country_data2, select=c(scientificName,decimalLatitude,decimalLongitude,datasetKey,taxonRank,species,year))
            format_country2 = rbind(format_country2,country_subset)
          }
          
        }
        country_only_species2 = format_country2[is.na(format_country2$species) == FALSE,] #removes non-species records
        assign(con_code2,country_only_species2)
      }   else {
        country_only_species2 = eval(as.name(con_code2))
      }
      country_only_species3 = country_only_species2[which(country_only_species2$year %in% sampling),]
      species_richness_3 = length(unique(country_only_species3$species)) # obtains no. of unique species from an entire country's data on a particular kingdom
      richness_table[con,9] = species_richness_3
  }
  con = con + 1
}
    
plant_biotic_features_bigrange = richness_table
names(plant_biotic_features_bigrange) = c("articles","country","location","latitude","longitude","plantae_SR_20yrecords_under75km","plantae_SR_20ydatasets_meanunder75km","iso2c","plantae_20ySR_country")
plant_biotic_features_bigrange = rbind(rep(NA,length=ncol(plant_biotic_features_bigrange)),plant_biotic_features_bigrange)
write.csv(plant_biotic_features_bigrange,"C:\\Users\\Jamila\\Documents\\better_plant_biotic_data.csv")

######## data for animals speices richness 
# re-run the above code but with 
animalia_key = name_backbone("Animalia")
animalia_key = animalia_key %>%
  dplyr::select(kingdomKey)
animalia_key = animalia_key[[1]]
taxon = animalia_key

# Some of the above code was corrected manually and reverted back to the general state
an_biotic_features_bigrange = richness_table
names(an_biotic_features_bigrange) = c("articles","country","location","latitude","longitude","an_SR_20yrecords_under75km","an_SR_20ydatasets_meanunder75km","iso2c","an_20ySR_country")
an_biotic_features_bigrange = rbind(rep(NA,length=ncol(an_biotic_features_bigrange)),an_biotic_features_bigrange)
 
# Save data  
factors_df = cbind(factors_df,subset(plant_biotic_features_bigrange,select=c(plantae_SR_20yrecords_under75km,plantae_SR_20ydatasets_meanunder75km,plantae_20ySR_country)),
                   subset(an_biotic_features_bigrange,select=c(an_SR_20yrecords_under75km,an_SR_20ydatasets_meanunder75km,an_20ySR_country)),plant_biotic_features_bigrange$iso2c)
write.csv(factors_df,"C:\\Users\\Jamila\\Documents\\factors_df.csv",row.names=FALSE)
# This data was added to the complete dataframe with compounds manually, outside of R

################### first instance of invasive in a country - from the end of final_normalisation.R ########################
final_df_read = readRDS("final_plant_features.rds")

final_df_read$Species[which(final_df_read$Species %in% "Melaleuca quinquenervia ")] = "Melaleuca quinquenervia"
final_df_read = cbind(final_df_read,rep(NA,length=nrow(final_df_read)))
species = unique(final_df_read$Species)
#species = species[-1]

for(spec in species) {
  #spec = species[1]
  
  taxon_key = name_backbone(spec)
  taxon_key = taxon_key %>%
    dplyr::select(speciesKey)
  taxon_key = taxon_key[[1]]

  countries = unique(na.omit(final_df_read$iso2c[which(final_df_read$Species %in% spec)]))
  spec_data = occ_data(taxonKey = taxon_key,hasCoordinate = TRUE, country=c(countries))
  
  
  first_year_by_country = NULL
  iso2cs = names(spec_data)
  for(iso2c in iso2cs) {
    country_data = spec_data[[iso2c]]
    years = unique(country_data$data$year)
    early_year = min(as.numeric(na.omit(years)))
    
    if(early_year == Inf) {
      first_year_by_country = rbind(first_year_by_country,c(iso2c,"not_recorded"))
    } else {
      first_year_by_country = rbind(first_year_by_country,c(iso2c,early_year))
    }
    
  }
  
  for(rec in 1:nrow(first_year_by_country)) {
    iso = first_year_by_country[rec,1]
    first_rec = first_year_by_country[rec,2]
    ind = which(final_df_read$Species == spec & final_df_read$iso2c == iso)
    final_df_read$`rep(NA, length = nrow(final_df_read))`[ind] = first_rec
  }
  assign(spec,first_year_by_country)
}

colnames(final_df_read)[35] = "1st_record_of_invasive"
saveRDS(final_df_read,file="final_plant_features.rds")

