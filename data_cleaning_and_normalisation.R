
library(tidyverse)
library(tidyverse)
library(webchem)
library(taxize)
library(plyr)
library(dplyr)
library(rJava)
library(rcdk)

features_df = read.csv("~/Plant_Features_Database.csv") # created outside of R
species_nonorm = unique(as.character(features_df[,1]))

################################################ Data Tidying ###############################################

df_no_spha = features_df[which(features_df$Species != species_nonorm[12]),] #removes Sphagneticola trilobata as I have not finshed formatting it yet
df_no_spha = df_no_spha[which(df_no_spha$Species != species_nonorm[13]),] #removes empty values included in the csv for some reason

# Removing datasets where location was not found
df_no_spha$Location[df_no_spha$Location == ""] = NA # add NAs in case some rows were blank
rowskeep = which(is.na(df_no_spha$Location)==FALSE)
df_no_spha = df_no_spha[rowskeep,]

######################################################## Data Normalisation ###################################################


# part = vector of characters, all the non-normalised name of one feature value
# replacer = character string of the normalised value to replace the non-normalised value
# df = dataframe being normalised
# header = column name/name of plant feature as a character string
part_norm = function(part,replacer,df,header) { 
  part_index = c()
  for(p in part) {
    part_index = append(part_index,which(df[[header]] == p)) # obtain row indexes that contain the data to be normalised. Creates a vectof
  }
  for(p2 in part_index) {
    df[[header]][[p2]] = replacer # replaces the data point of a given feature with the normalised data value
  }
  return(df) # returns the normalised dataframe
}

# Normalising Plant_Parts
parts_nonorm = unique(df_no_spha$Plant_Part)
df_no_spha = part_norm(c("aerial part","aerial parts","Aerial parts"),"aerial",df_no_spha,"Plant_Part")
df_no_spha = part_norm(c("leaves","leaf","Leaves"),"leaf",df_no_spha,"Plant_Part")
df_no_spha = part_norm(c("flower","flowers"),"flower",df_no_spha,"Plant_Part")
df_no_spha = part_norm(c("fruit","Fruit"),"fruit",df_no_spha,"Plant_Part")
df_no_spha = part_norm(c("root and leaves","leaves and flowers","flower and leaves"),"whole plant",df_no_spha,"Plant_Part")
df_no_spha = part_norm(c("manufactured leaf oil","extracted oil"),"manufactured oil",df_no_spha,"Plant_Part")

# Normalising invasivity, wild status and habitat
unique(df_no_spha$invasivity)
df_no_spha = part_norm(c("native","native "),"native",df_no_spha,"invasivity")
df_no_spha = part_norm(c("not_introduced","cultivated, not invasive","cutlivated, not invasive"),"cultivated",df_no_spha,"invasivity")
df_no_spha = part_norm(c("cultivated or invasive"),"invasive",df_no_spha,"invasivity")

wild_val = unique(df_no_spha$wild.not_wild)
df_no_spha[which(df_no_spha$wild.not_wild == wild_val[3]),] # val = cultivated?
df_no_spha = part_norm("cultivated?","not_wild",df_no_spha,"wild.not_wild")
df_no_spha[which(df_no_spha$wild.not_wild == wild_val[6]),] # val = ?
df_no_spha = part_norm("?","not_wild",df_no_spha,"wild.not_wild")
df_no_spha = part_norm("city - maybe garden","urban",df_no_spha,"Habitat") # where "?" was filled in for the wild value
df_no_spha[which(df_no_spha$wild.not_wild == wild_val[4]),] # val = NA # apparently there are none?
df_no_spha[which(df_no_spha$wild.not_wild == wild_val[5]),] 
df_no_spha = part_norm("","NA",df_no_spha,"wild.not_wild")
df_no_spha = part_norm("NA",NA,df_no_spha,"wild.not_wild")


# Normalising Plant Names

plant_names = unique(df_no_spha$Species)

lan_cam_names = plant_names[1:4]
tsn_lancam = get_tsn(lan_cam_names, accepted = FALSE) # none were found, do a manual search. I think these were part of the formatting in EssoilDB rather than names
df_no_spha = part_norm(plant_names[1:4],"Lantana camara",df_no_spha,"Species")

mel_qui_names = plant_names[5:8]
tsn_melqui = get_tsn(mel_qui_names, accepted = FALSE) # most weren't found. Will do a manual research check
lapply(tsn_melqui, itis_acceptname) # error - not enough accepted names?
# manual search suggests other names are synonamous 
df_no_spha = part_norm(mel_qui_names,"Melaleuca quinquenervia",df_no_spha,"Species")

psi_cat_names = c("Psidium cattleianum","Psidium cattleianum sabine")
tsn_psicat = get_tsn(psi_cat_names,accepted=FALSE) # Psidium cattleianum sabine not found, indicates it is not an accepted name
df_no_spha = part_norm(plant_names[9:11],"Psidium cattleianum",df_no_spha,"Species")

############################################### More Data Tidying ##############################################################
# Remove df annotations
df_no_spha = df_no_spha[,1:12]

# Add NAs to sampling_date
df_no_spha$sampling_date[df_no_spha$sampling_date == ""] = NA

#Save our data
write.csv(df_no_spha,"C:\\Users\\Jamila\\Documents\\plant_features_database_no_Sphagneticola.csv",row.names=FALSE)

################################## Some more data normalisation with the updated dataframe 
df_no_spha = read.csv("~/plant_features_database_no_Sphagneticola.csv")

# Check new data is normalised
unique(df_no_spha$Plant_Part)
df_no_spha = part_norm(c("aerial part"),"aerial",df_no_spha,"Plant_Part")
unique(df_no_spha$Habitat)
df_no_spha = part_norm(c(""),NA,df_no_spha,"Habitat")
df_no_spha = part_norm(c("urban","town"),"urban",df_no_spha,"Habitat")

write.csv(df_no_spha,"C:\\Users\\Jamila\\Documents\\plant_features_database_no_Sphagneticola.csv",row.names=FALSE)
df_no_spha = read.csv("~/plant_features_database_no_Sphagneticola.csv")
#Factorise the appropriate data - when it was all factorised a lot of bad things happened
to_factor = c("Species","Plant_Part","Country","invasivity","wild.not_wild","Habitat")
for(var in to_factor) {
  df_no_spha[[var]] = as.factor(df_no_spha[[var]])
}

write.csv(df_no_spha,"C:\\Users\\Jamila\\Documents\\plant_features_database_no_Sphagneticola.csv",row.names=FALSE)


######################################### Checks on up-to-date df ########################################
plant_features = read.csv("~/plant_features_database_no_Sphagneticola.csv")
# check to see if up to data df needs normalising 
for(col in 1:ncol(plant_features)) {
  print(unique(plant_features[,col]))
}

# invasivity needs normalising
# a spelling mistake in some of the plant names
plant_features$Species[which(plant_features$Species == "Psidium catteianum")] = "Psidium cattleianum"  
plant_features$invasivity[which(plant_features$invasivity== "cultivated")] = "introduced, not invasive" 
write.csv(plant_features, "C:\\Users\\Jamila\\Documents\\plant_features_database_no_Sphagneticola.csv")


################################## Compound Normalisation ####################################################

plant_features = read.csv("~/plant_features_database_no_Sphagneticola.csv")


# normalise data to SMILES
# cir_query() converts identifier = common name to representation = Isomeric SMILES
#  note - some compound common names have multiple SMILES
# returns NA if no smiles is found, so can fill out a table

compounds = plant_features$Chemical

SMILES = data.frame(matrix(NA,ncol=30,nrow = length(compounds))) #creates a too large NA-filled dataframe so different length rows can be added to the datafrane
chem_no = 0 #start here
while(chem_no < (length(compounds))+1) { # create a while loop cos of dodgy internet. If the internet cuts out the dataframe will not start from the beginning
  chem = compounds[chem_no]
  chem_list = NULL
  chem_list = chem
  can_smi = cir_query(identifier = chem, representation = "smiles") #variable takes a named list form, name = common name (identifier), value = smiles (representation)
  can_smi = can_smi[[1]] # assigns the variable with the smiles (or NA) value only
  chem_list = append(chem_list, can_smi) # creates list containing the common name and all canonical smiles (1 compound/common name can have more than 1)
  SMILES[chem_no,1:length(chem_list)] = chem_list #adds list to the dataframe from the dataframe's 1st col
  chem_no = chem_no + 1
}

write.csv(SMILES,"C:\\Users\\Jamila\\Documents\\compounds_norm.csv",row.names=FALSE)
SMILES = read.csv("~/compounds_norm.csv")

# clean data - remove excess columns 
SMILES_2 = SMILES #to make sure not all data is lost
# loop to remove empty columns (containing on NAs) from the dataframe
for(colu in 1:ncol(SMILES_2)) {
  na_check = length(which(is.na(SMILES_2[,colu]))) # identifies the no. of NA values in a column
  if(na_check == length(SMILES_2[,colu])) { # if all values in the column = NA
    SMILES_2 = SMILES_2[,1:(colu-1)] # removes column and all the columns after
    break #end the loop here
  }
}

# Chemicals with multiple SMILES appear to be repeats of the same SMILES
# Check code to see if this is true
# If true, additionally SMILES can be removed
multi_SMILE_ind = which(is.na(SMILES_2[,3]) == FALSE)
no_rep_BOOl = 0

for(index in multi_SMILE_ind) {
  if(SMILES_2[index,2] == SMILES_2[index,3]) {
    SMILES_2[index,3] = NA # sets redundant SMILES to NA in the 3rd column, removing repeats
  }
  else if(SMILES_2[index,2] != SMILES_2[index,3]) {
    no_rep_BOOl = no_rep_BOOl + 1 # will give the total chemical names with non-repeated multiple SMILEs
  }
  
}
if(no_rep_BOOl == 0) {
  SMILES_2 = SMILES_2[,1:2] #if all Canonical smiles were repeats, the 3rd column was removed (so only the first smile was kept)
}
# some smiles are not repeats
# some of these might be for cis or trans isomers 
# check the data to see if sometimes cis and trans isomers (or whatever) are specified separately or whether they are always given in the non-specifical multi-SMILES form
# if the data is always non-specific, pick one SMILES and set that to be the smiles for the named chemical compound.
# if the data is specific, the data needs to be checked and some isomers need to be added together to standardise the data.
# the random selection of the SMILES may affect the chemical emission matrices (alter chemical distances),  so retain both standardised and non-standardised datasets and ask Gita

unique_multi_SMILES = which(is.na(SMILES_2[,3]) == FALSE)
for(index2 in unique_multi_SMILES) {
  list_SMILES = vector(mode="list")
  SMILE1 = SMILES_2[index2,2]
  SMILE2 = SMILES_2[index2,3]
  SMILE1_com = c()
  SMILE2_com = c()
  for(i in 1:nrow(SMILES_2)) {
    CAN_row2 = SMILES_2[i,2]
    CAN_row3 = SMILES_2[i,3]
    if(is.na(CAN_row2) == FALSE & CAN_row2 == SMILE1 || is.na(CAN_row3) == FALSE & CAN_row3 == SMILE1) { # this code is wrong here 
      SMILE1_com = append(SMILE1_com,SMILES_2[i,1])
    }
    else if(is.na(CAN_row2) == FALSE & SMILES_2[i,2] == SMILE2 || is.na(CAN_row3) == FALSE & SMILES_2[i,3] == SMILE2) {
      SMILE2_com = append(SMILE2_com,SMILES_2[i,1])
    }
    compound_name = SMILES_2[index2,1]
    a_list = list(length = 2)
    a_list[[1]] = SMILE1_com
    a_list[[2]] = SMILE2_com
    assign(compound_name,a_list) # I'm not sure if this works
  }
}

# from inspection, the first SMILES given in SMILES_2[,2] are the unique values. The second column will be the default SMILE value
# Exceptions (1): 
## Cymol: Cymol and Thymol are the same in SMILES_2[,2]
cymol_data = c(which(SMILES_2[,1] == "Cymol"),which(SMILES_2[,1] == "Thymol"),which(SMILES_2[,1] == "thymol"))
for(c in cymol_data) {
  sample = plant_features$Profile_no[[c]]
  print(sample)
}
# cymol and thymol are only found in one dataset together. These points should be added together.
# I will do this manually after the SMILES data is added to the dataframe - it'll be easier
SMILES_2 = SMILES_2[,1:2]
write.csv(SMILES_2, "C:\\Users\\Jamila\\Documents\\compounds_norm.csv",row.names=FALSE)

# Inspecting compounds where SMILEs could not be obtained
no_SMILE = as.data.frame(unique(SMILES_2[,1][is.na(SMILES_2[,2])]))
no_SMILE_ind = which(is.na(SMILES_2[,2]))
write.csv(no_SMILE, "C:\\Users\\Jamila\\Documents\\compounds_noSMILE.csv",row.names=FALSE) # data will need to be input manually

# from manual sampling of SMILEs, SMILEs are now a mixture of isomeric and canonical SMILEs
# this is bad for analysis and plotting individual compounds
# RDKIT accepts these values fine, however, so it is ok for Morgan Circular Fingerprint generation

no_SMILE = read.csv("~/compounds_noSMILE.csv")
no_SMILE = no_SMILE[,1:3]
SMILES_2 = read.csv("~/compounds_norm.csv")
SMILES_2 = cbind(SMILES_2,rep(NA,length=nrow(SMILES_2)),rep(NA,length=nrow(SMILES_2)),rep(NA,length=nrow(SMILES_2)))


# 1) convert to isomeric SMILES
new_can_no_SMILE = as.data.frame(matrix(NA,nrow=nrow(no_SMILE),ncol=10))
chem_no2 = 1
while(chem_no2 < (nrow(no_SMILE)+1)) {
  chem1 = no_SMILE[chem_no2,1]
  chem2 = no_SMILE[chem_no2,2]
  chem3 = no_SMILE[chem_no2,3]
  chem_list1 = chem1
  chem_list2 = c(chem2,chem3)
  chem_list3 = NULL
  for(can in chem_list2) {
    if(is.na(can) == FALSE) {
      chem_list1 = append(chem_list1, can)
      can_smi = cir_query(identifier = can, representation = "smiles") #variable takes a named list form, name = common name (identifier), value = smiles (representation)
      can_smi = can_smi[[1]] 
      chem_list3 = append(chem_list3, can_smi)
    }
  }
  chem_list1 = append(chem_list1,chem_list3)
  new_can_no_SMILE[chem_no2,1:length(chem_list1)] = chem_list1 
  chem_no2 = chem_no2 + 1
}

# check to see if isomeric SMILES were not found 
cc = c()
for(i in 1:nrow(new_can_no_SMILE)) {
  if(is.na(new_can_no_SMILE[i,2]) == FALSE & is.na(new_can_no_SMILE[i,2]) == TRUE) {
    cc = append(cc,new_can_no_SMILE[i,1])
  }
}
# cc = NULL, all SMILES were found
# however, not all compounds have isomeric SMILEs so returned Canonical - this will need to be redone (see lines 420-465)
new_can_no_SMILE[,2] == no_SMILE[,2]

# add the new data to SMILES_3
new_can_no_SMILE[96,1] = "germacrene" # I changed the name without thinking 
SMILES_3 = SMILES_2

multi_smile_ind_2 = which(is.na(no_SMILE[,3]) == FALSE)
for(index in 1:nrow(new_can_no_SMILE)) {
  compound = new_can_no_SMILE[index,1]
  where_com = which(SMILES_3[,1] == compound)
  if(index %in% multi_smile_ind_2 == FALSE) {
    SMILES_3[where_com,2] = new_can_no_SMILE[index,3]
  } else {
    SMILES_3[where_com,2:3] = new_can_no_SMILE[index,4:5]
  }
}
write.csv(SMILES_3, "C:\\Users\\Jamila\\Documents\\compounds_norm.csv",row.names=FALSE)


###############################################################################################################
plant_features = read.csv("~/plant_features_database_no_Sphagneticola.csv")

################################ more normalisation - binary data conversion ######################################

# convert invasivity into 2 columns, 1 for invasive status, one for native/non-native
# col 1 = native status. If a plant is native (0) or not (1) to a region
# col 2 = invasive status. If a plant is invasive (1) or not (0) in a region

invasive_native = as.data.frame(matrix(NA,nrow=0,ncol=2))

for(r in 1:nrow(plant_features)) {
  status = plant_features$invasivity[[r]]
  if(status == "invasive") {
    val = c(1,1)
  } else if(status == "native") {
    val = c(0,0)
  } else {
    val = c(1,0)
  }
  invasive_native = rbind(invasive_native,val)
}
names(invasive_native) = c("native_status","invasive_status")

# convert wild/not_wild into a single binary column 
# wild = 1, not_wild = 0

wild_status = as.data.frame((matrix(NA,ncol=1,nrow=0)))
for(r in 1:nrow(plant_features)) {
  status = plant_features$wild.not_wild[[r]]
  if(is.na(status) == FALSE) {
    if(status == "wild") {
      val = 1
    } else {
      val = 0
    }
    wild_status = rbind(wild_status,val)
  } else if(is.na(status) == TRUE) {
    wild_status = rbind(wild_status,status)
  }
}
names(wild_status) = "wild.not_wild"

# convert plant part into either a bit vector or into many columns - do both I think

plant_part = unique(plant_features$Plant_Part)

# 1) bit vector 
ppart_bit = as.data.frame(matrix(NA,ncol=nrow(plant_features),nrow=1))

for(r in 1:nrow(plant_features)) {
  bit_vec = c(rep(0,length=length(plant_part)))
  part = plant_features$Plant_Part[[r]]
  part_ind = which(plant_part == part)
  bit_vec[part_ind] = 1
  ppart_bit[[r]] = list(c(bit_vec))
}

ppart_bit_trans = t(ppart_bit)
ppart_bit_trans = as.data.frame(ppart_bit_trans,row.names = FALSE)
names(ppart_bit_trans) = "Bit vector plant part"

# 2) many columns
ppart_col = as.data.frame(matrix(NA,ncol=1,nrow=0))

for(r in 1:nrow(plant_features)) {
  bit_vec = c(rep(0,length=length(plant_part)))
  part = plant_features$Plant_Part[[r]]
  part_ind = which(plant_part == part)
  bit_vec[part_ind] = 1
  ppart_col = rbind(ppart_col,bit_vec)
}
names(ppart_col) = plant_part  

##################################### Compiling the dataset

factors_df = read.csv("~/factors_df.csv")
SMILES_dat = read.csv("~/compounds_norm.csv")

plant_features = plant_features[,2:ncol(plant_features)]

final_df = cbind(plant_features[,1:2],SMILES_dat,ppart_bit_trans,plant_features[,6:8],invasive_native,plant_features[,9],wild_status)
final_df = subset(final_df, select = -c(V5))
colnames(final_df)[c(3,4,11)] = c("Compound", "SMILES", "Location")
final_df = cbind(final_df,replicate((ncol(factors_df)-8),rep(NA,length=nrow(final_df))))

# temp_df = as.data.frame(unique(plant_features$Location)) - check to make sure spellings match up
factors_df = subset(factors_df,select=-c(Agrionomic.land.use...agriculture.or.settlements.,Soil.Ammonium,Soil.Nitrate,Grazing.prone.,plantae_SR__records_under75km,plantae_SR_datasets_meanunder75km,
                                         plantae_SR_country,an_SR__records_under75km,an_SR_datasets_meanunder75km,
                                         an_SR_country,plantae_SR_20ydatasets_meanunder75km,an_SR_20ydatasets_meanunder75km))

for(ind in 1:nrow(factors_df)) {
  locat = factors_df$locations[[ind]]
  
  if(is.na(locat) == FALSE) {
    emission_ind = which(final_df$Location %in% locat)
    correct_data = factors_df[ind,5:ncol(factors_df)]
    correct_data = subset(correct_data,select=-c(wild.not_wild,Habitat))
    final_df[emission_ind,13:(length(correct_data)+13)] = correct_data
    
  } else { 
    emission_ind = which(final_df$source_article %in% factors_df$articles[[ind]])
    
  }
  
  correct_data = factors_df[ind,5:ncol(factors_df)]
  correct_data = subset(correct_data,select=-c(wild.not_wild,Habitat))
  final_df[emission_ind,13:(length(correct_data)+12)] = correct_data
  
} 

final_df = final_df[,1:34]
colnames(final_df)[c(13:34)] = c(colnames(factors_df[5]),colnames(factors_df[8:ncol(factors_df)]))

# save dataset
saveRDS(final_df,file="final_plant_features.rds")


### THIS IS FROM AFTER INITIAL_ANALYSIS.R AND MORE ANALYSIS FOR GITA (ANLYSIS 2).R: error - isomeric and Canonical Smiles have not been converted properly
# convert all molecules to Canonical Smiles, as all molecules have them (not all have isomeric)
# run initial_analysis.R and more analysis for gita (analysis 2).R after this


smi1 = parse.smiles("CC(C)[C@@H]1CC[C@](C)(O)[C@@H]2CCC(=C[C@@H]12)C")[[1]]
get.smiles(smi1,smiles.flavors("Canonical"))
smi2 = parse.smiles("CC(C)C1CCC(C)(O)C2CCC(=CC12)C")[[1]]
get.smiles(smi2,smiles.flavors("Canonical"))

comp_feat = readRDS("final_plant_features.rds")
init_SMILES = comp_feat$SMILES
new_CANSMI = as.data.frame(matrix(NA,nrow=length(init_SMILES),ncol=1))
rep_SMILES = NULL

for(smile in init_SMILES) {
  if(length(which(rep_SMILES == smile)) == 0) {
    smi_parse = parse.smiles(smile)[[1]]
    can_smi = get.smiles(smi_parse,smiles.flavors("Canonical"))
    smi_ind = which(init_SMILES %in% smile)
    new_CANSMI[smi_ind,1] = can_smi
    rep_SMILES = append(rep_SMILES,smile)
  }
}
# some NAs created - compare data 
new_CANSMI = cbind(init_SMILES,new_CANSMI)
for(ind2 in 1:nrow(new_CANSMI)) {
  smile = new_CANSMI[ind2,1]
  smile2 = new_CANSMI[ind2,2]
  if(is.na(smile) == FALSE & is.na(smile2) == TRUE) {
    smi_parse = parse.smiles(smile)[[1]]
    #print(smi_parse)
    can_smi = get.smiles(smi_parse,smiles.flavors("Canonical"))
    #print(can_smi)
    new_CANSMI[ind2,2] = can_smi
    
  }
}

# ^^  package function do not like NAs and disrupted querying for other functions, now fixed
# make a new comp_feat in case the conversion was wrong
final_df2 = comp_feat
final_df2$SMILES = new_CANSMI[,2]
# actually name comp_feat, but save as a different RDS file
comp_feat = final_df2
saveRDS(comp_feat,"final_plant_feature2.rds")


# removing bad samples
comp_feat = comp_feat[-1,] # remove first profile 

# remove a profile known to have only 2 compounds, name of the profile unknown
samples = unique(comp_feat$Profile_no)
for(s in samples) {
  ind = c(which(comp_feat$Profile_no %in% s))
  l_ind = length(ind)
  if(l_ind == 2) {
    print(s)
  }
}

# s=56
comp_feat = comp_feat[which(comp_feat$Profile_no != 56),]


####### create a csv file for the dataframe 
comp_feat_for_csv = comp_feat
comp_feat_for_csv$`Bit vector plant part`

null_data = as.data.frame(matrix(NA,ncol=6,nrow=nrow(comp_feat_for_csv)))
for(row in 1:length(comp_feat_for_csv$`Bit vector plant part`)) {
  data = comp_feat_for_csv$`Bit vector plant part`[[row]]
  null_data[row,] = data
}

names(null_data) = unique(plant_features$Plant_Part)

comp_feat_for_csv = cbind(comp_feat_for_csv[,1:5],null_data,comp_feat_for_csv[,7:ncol(comp_feat_for_csv)])
write.csv(comp_feat_for_csv,"C:\\Users\\Jamila\\Documents\\convertable_finalfeatures2.csv")



##################################### 

# problem - some compounds within the same profile have been given the same SMILES
# check data to see common name of each compound to see if they can be grouped together into 2 categories, and manually find the appropriate SMILES from PubChem
# smiles are still technically correct, so chemical distances should be ok. Rerun if we have time just in case.

rep_SMILES = NULL

for(pro in profiles) {
  
  profile_SMI = comp_feat_for_csv$SMILES[which(comp_feat_for_csv$Profile_no %in% pro)]
  dup_SMI = unique(profile_SMI[duplicated(profile_SMI)])
  
  for(dup in dup_SMI) {
    comp_name = comp_feat_for_csv$Compound[which(comp_feat_for_csv$Profile_no %in% pro & comp_feat_for_csv$SMILES %in% dup)]
    
    for(comp in comp_name) {
      subdata = c(pro,dup,comp)
      rep_SMILES = rbind(rep_SMILES,subdata)
      
    }
  }
}

rep_SMILES = as.data.frame(rep_SMILES)
names(rep_SMILES) = c("Profile_no","SMILES","Compound")

# remove repeats of the same spelling
comp_feat_csv_cleaned = comp_feat_for_csv
rep_SMILES2 = rep_SMILES

for(row in 1:nrow(rep_SMILES2)) {
  pro = rep_SMILES$Profile_no[[row]]
  comp = rep_SMILES$Compound[[row]]
  
  ind = which(comp_feat_csv_cleaned$Profile_no %in% pro & comp_feat_csv_cleaned$Compound %in% comp)
  ind = ind[-1]
  
  if(length(ind) > 0) {
    comp_feat_csv_cleaned = comp_feat_csv_cleaned[-c(ind),]
    rep_SMILES2 = rep_SMILES2[-c(row,row+1),]
  }
}

# removing tau-cadinol from 1 dataset
ind_tau = which(comp_feat_csv_cleaned$Profile_no %in% 4 & comp_feat_csv_cleaned$Compound %in% "tau-cadinol")
comp_feat_csv_cleaned = comp_feat_csv_cleaned[-c(ind_tau),]

# removing any NAs
rep_SMILES2 = rep_SMILES2[-c(which(is.na(rep_SMILES2$SMILES)==TRUE)),]

# obtaining unique compounds to manually find new smiles
comp_search = unique(rep_SMILES2$Compound)
comp_search = sort(comp_search)
comp_search = as.data.frame(comp_search)
write.csv(comp_search,"C:\\Users\\Jamila\\Documents\\smiles_search.csv")

# some isomers do not have unique SMILES - combine their rows if so.
# in fact, after manual inspection, all SMILES are non-unique in PubChem

 
################### Data Matrix Construction

# presence/absence
# version 1 - each compound has a column, presence(1)/absence(0) is filled in 
comp_feat_for_csv = read.csv("~/convertable_finalfeatures2.csv")

compounds = unique(comp_feat_for_csv$SMILES)
samples = unique(comp_feat_for_csv$Profile_no)


new_com_format = as.data.frame(matrix(NA,ncol=length(compounds)+1,nrow=length(samples)))
names(new_com_format) = c("Profile_no",compounds)
new_features = as.data.frame(matrix(NA,ncol=ncol(comp_feat_for_csv)-4,nrow=length(samples)))
names(new_features) = names(comp_feat_for_csv)[!names(comp_feat_for_csv) %in% c("Compound","Chemical_family","SMILES","percentage")]
comp_feat_no_chem = subset(comp_feat_for_csv,select=-c(Compound,Chemical_family,SMILES,percentage))

for(sam_ind in 1:length(samples)) {
  
  sam = samples[sam_ind]
  
  row_data = comp_feat_no_chem[which(comp_feat_no_chem$Profile_no %in% sam)[[1]],]
  new_features[sam_ind,] = row_data 
  
  new_com_format$Profile_no[[sam_ind]] = sam
  profile = comp_feat_for_csv$SMILES[which(comp_feat_for_csv$Profile_no %in% sam)]
  for(comp in compounds) {
    if(is.na(comp)==FALSE) {
      if(length(which(profile == comp)) == 0) {
        new_com_format[[comp]][sam_ind] = 0
      } else if(length(which(profile == comp)) > 0) {
        new_com_format[[comp]][sam_ind] = 1
      }
    }
  }
}

data_clean = c(which(colnames(new_com_format) == "O"), which(is.na(colnames(new_com_format)) == TRUE))
new_com_format = new_com_format[,-c(data_clean)]
matrix_data_feat = cbind(new_features[,2:3],new_com_format[,2:ncol(new_com_format)],new_features[,4:ncol(new_features)])
saveRDS(matrix_data_feat,"matrix_format_features.rds")
matrix_data_feat = readRDS("matrix_format_features.rds")



#### % amounts

# clean data a bit more to remove NAs
comp_feat_csv_cleaned = comp_feat_csv_cleaned[which(is.na(comp_feat_csv_cleaned$SMILES)==FALSE),]

# need to identify trace amounts - set to NA
trace_val = sort(decreasing=TRUE,unique(comp_feat_csv_cleaned$percentage))[1:4]
comp_feat_csv_cleaned$percentage[which(comp_feat_csv_cleaned$percentage %in% trace_val)] = NA
# convert data to numeric
comp_feat_csv_cleaned$percentage = as.numeric(comp_feat_csv_cleaned$percentage)
saveRDS(comp_feat_csv_cleaned, "comp_feat_csv_cleaned.rds")


## for non-unique SMILES in a profile, combine all rows with the same SMILES 
matrix_feat_percent = matrix_data_feat
compound_col = which(colnames(matrix_feat_percent) %in% compounds)

for(row in 1:nrow(matrix_feat_percent)) {
  pro = matrix_feat_percent$Profile_no[[row]] 
  percentage_smile = subset(comp_feat_csv_cleaned[which(comp_feat_csv_cleaned$Profile_no %in% pro),],select=c(SMILES,percentage))
  
  for(col in compound_col) {
    SMILE = colnames(matrix_feat_percent)[col]
    SMILE_bin = matrix_feat_percent[row,col]
    
    if(SMILE_bin == 1) {
      percent = percentage_smile$percentage[which(percentage_smile$SMILES %in% SMILE)]
      
      if(length(percent) > 1) {
        add_percent = sum(percent,na.rm=TRUE)
        matrix_feat_percent[row,col] = add_percent
        
      } else {
        matrix_feat_percent[row,col] = percent
        
      }
    }
  }
}

# some chemodiversity was lost, so do not do the same for previous analysis to make it less significant
saveRDS(matrix_data_feat,"matrix_feat_percent.rds")

################################################ Chemical distances ###########################################

# get chemical families from EssoilDB info_compound data files
info_compound = read.csv("~/info_compound_15th_may.csv")
SMILES = unique(comp_feat$SMILES)

comp_feat = subset(comp_feat,select=-c(Chemical_family))
comp_feat = add_column(comp_feat, d = rep(NA,length=nrow(comp_feat)), .after = "Compound")
colnames(comp_feat)[4] = "Chemical_family"

for(smi in SMILES) {
  com_names = comp_feat$Compound[which(comp_feat$SMILES %in% smi)]
  for(compound in com_names)  {
    fam_ind = which(info_compound$COMPOUND %in% compound)
    if(length(fam_ind) > 0) {
      chem_fam = info_compound$CHEMICAL_FAMILY[fam_ind][[1]]
      Smiles_ind = which(comp_feat$SMILES %in% smi)
      comp_feat$Chemical_family[Smiles_ind] = chem_fam
      break
    }
  }
}
saveRDS(comp_feat,"final_plant_features.rds")

######################### Chemical Distances Between SMILES
# function just does 1 - number across a dataframe/matrix/whatever
dis_sim_converter = function(measure) {
  new_measure = 1 - measure
  return(new_measure)
}

mols = parse.smiles(compounds) # parse smiles to get a molecule object
circ_fing = lapply(mols, get.fingerprint, type="circular") # fingerpints

# similarity matrix w. tanimoto metric again
fing_sim_matx = fingerprint::fp.sim.matrix(circ_fing,method="tanimoto")
heatmap(fing_sim_matx)

fing_sim_matx = as.data.frame(fing_sim_matx)
colnames(fing_sim_matx) = compounds
rownames(fing_sim_matx) = compounds

# distance matrix
fing_diss_matx = as.data.frame(dis_sim_converter(fing_sim_matx))

# save files 
saveRDS(fing_sim_matx,"fing_sim_matrx.rds")
saveRDS(fing_diss_matx,"fing_diss_matrx.rds")






