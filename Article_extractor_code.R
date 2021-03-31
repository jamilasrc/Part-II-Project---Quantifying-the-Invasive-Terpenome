library(XML)
library(rvest)
library(stringr)
library(dplyr)

######################################## Lam_cam #####################################################

##################################### Step 1 - convert an excel table with external hyperlinks to a file where R can read these hyperlinks ##############################################
# Needed for EssoilDB as hyperlinks are named not as the web address. 
# The hyperlinks in question are the "Details" hyperlink. This gives us a link to an html with the source article for each chemical compound. 
# We need to know the source data to find information about the plant's location and growing conditions to build the feature vectors. 
# Initial stages are done with L. camara 

lan.cam.excel = "lantana_camara_essoildb.xlsx"

# rename file to .zip
lan.cam.zipfile = sub("xlsx", "zip", lan.cam.excel)   # converts the .xslx file to a zip file
file.copy(from = lan.cam.excel, to = lan.cam.zipfile)  # I think this checks the data matches so that the conversion has been successful?

# unzip the file
unzip(lan.cam.zipfile) 

# unzipping produces a bunch of files which we can read using the XML package
# assume sheet1 has our data
xml_lancam = xmlParse("xl/worksheets/sheet1.xml")

# finally grab the hyperlinks
hyperlinks_lancam = xpathApply(xml_lancam, "//x:hyperlink/@display", namespaces="x")

############################################### Step 2 - create a function that extracts Article Name from the hyperlinked url ###########################################

article_extract = function(hyperlinks) { 
  
  article_titles = c()
  
  for(i in 1:length(hyperlinks)) {
    x = readLines(as.character(hyperlinks[[i]]))                                                # produces a list with html source details, which includes the Article title
    val = grep("Article Title",x)                                                               # Identifies the element no. of the list that contains "Article Title"
    val2 = gregexpr("Article Title",x[val])                                                     # Finds the beginning element of "Article Title" within the list element (the element is a string of sorts) to narrow down where the information is
    art_str = substring(x[val],as.numeric(val2),(as.numeric(val2)+200))                         # creates a string starting with the beginning of "Article Title" and ending 200 characters ahead of the start of the string. This should (hopefully) include the whole title + plus extra text
    art_str_red = str_match(art_str,"Article Title</b></td><td>\\s*(.*?)\\s*</td></tr><tr><td") # removes the surrounding text around the beginning and end of the article title. (\\s*(.*?)\\s* is the function code separating the boundaries of the text to be removed. Any characters within these boundaries is kept.)
    article_titles = append(article_titles,art_str_red[,2])                                     # adds the article title to a vector. (art_str_red[,2] needed as art_str_red[,2] produces a list of 2 & the 2nd element is the extracted string)
    print(length(article_titles))                                                               # keeps track of how far through the data the function is         
  }
  return(article_titles)
}

titles_lamcam = article_extract(hyperlinks_lancam)

########################################### Step 3 - produce a dataframe of chemical emissions and article titles - the beginning of the features dataframe ####################################

lam_cam = read.csv("~/lantana_camara_essoildb.csv")               # obtains a df taken directly from EssoilDB 
lam_cam = lam_cam[,1:5]                                           # removes "Details" column - contains hyperlink, but in the .csv confers no info
lam_cam = cbind(lam_cam[1:3017,],titles_lamcam)                   # adds the article names for each chemical compound to the df
write.csv(lam_cam,"C:\\Users\\Jamila\\Documents\\lam_cam.csv",row.names=FALSE)  # writing a new df saves the df (so the article_extract() function does not need to be run again - it has a long runtime)

lam_cam = read.csv("~/lam_cam.csv")                               # obtains the saved df for further editing. Here I added some countries based off the article titles manually 


############################################# Step 4 - Tidy data ################################ 
# Remove duplicate rows - inspection of the data indicates some data sources have been added multiple times
library(dplyr)
no_rep_lam_cam = distinct(lam_cam) # creates new df of only unique rows, so removes any repeat data
# 1/3 data was removed in the above step. 
names(no_rep_lam_cam)[names(no_rep_lam_cam)=="titles_lamcam"] = "source_article"   # renaming column with source papers as "source_articles"
write.csv(no_rep_lam_cam,"C:\\Users\\Jamila\\Documents\\lam_cam_no_repeats.csv",row.names=FALSE)  # export df to add info manually after reading through papers. Also identify sources of profiles where source_article == NA.

length(unique(no_rep_lam_cam$source_article)) 
# output = 31, 31 papers to look through
# create df of all the papers that need reading, and a column that indicates if they have been read or not
papers_to_read_lamcam = as.data.frame(unique(no_rep_lam_cam$source_article))
papers_to_read_lamcam = cbind(papers_to_read_lamcam,rep("no",times=length(papers_to_read_lamcam)))
names(papers_to_read_lamcam) = c("paper","read?")
write.csv(papers_to_read_lamcam,"C:\\Users\\Jamila\\Documents\\source_papers_lamcam.csv",row.names=FALSE)


################################################ M. quinquenervia #######################################################

mel.qui.excel = "mel_qui.xlsx"
mel.qui.zipfile = sub("xlsx", "zip", mel.qui.excel)   
file.copy(from = mel.qui.excel, to = mel.qui.zipfile)  
unzip(mel.qui.zipfile) 
xml_melqui = xmlParse("xl/worksheets/sheet1.xml")
hyperlinks_melqui = xpathApply(xml_melqui, "//x:hyperlink/@display", namespaces="x")

titles_melqui = article_extract(hyperlinks_melqui)

mel_qui = read.csv("~/mel_qui.csv")               
mel_qui = mel_qui[,1:5]                                           
mel_qui = cbind(mel_qui,titles_melqui)
names(mel_qui)[names(mel_qui)=="titles_melqui"] = "source_article"
write.csv(mel_qui,"C:\\Users\\Jamila\\Documents\\mel_qui.csv",row.names=FALSE)  

library(dplyr)
no_rep_mel_qui = distinct(mel_qui) # creates new df of only unique rows, so removes any repeat data
write.csv(no_rep_mel_qui,"C:\\Users\\Jamila\\Documents\\mel_qui_no_repeats.csv",row.names=FALSE)  # export df to add info manually after reading through papers. Also identify sources of profiles where source_article == NA.

papers_to_read_melqui = as.data.frame(unique(no_rep_mel_qui$source_article))
papers_to_read_melqui = cbind(papers_to_read_melqui,rep("no",times=length(papers_to_read_melqui)))
names(papers_to_read_melqui) = c("paper","read?")
write.csv(papers_to_read_melqui,"C:\\Users\\Jamila\\Documents\\source_papers_melqui.csv",row.names=FALSE)


###################################################################### P. cattleainum #########################################

psi.catt.excel = "psi_catt.xlsx"
psi.catt.zipfile = sub("xlsx", "zip", psi.catt.excel) 
file.copy(from = psi.catt.excel, to = psi.catt.zipfile)  
unzip(psi.catt.zipfile) 
xml_psicatt = xmlParse("xl/worksheets/sheet1.xml")
hyperlinks_psicatt = xpathApply(xml_psicatt, "//x:hyperlink/@display", namespaces="x")
titles_psicatt = article_extract(hyperlinks_psicatt)

psi_catt = read.csv("~/psi_catt.csv")               
psi_catt = psi_catt[,1:5]                                           
psi_catt = cbind(psi_catt,titles_psicatt)
names(psi_catt)[names(psi_catt)=="titles_psicatt"] = "source_article"
write.csv(psi_catt,"C:\\Users\\Jamila\\Documents\\psi_catt.csv",row.names=FALSE)  

no_rep_psi_catt = distinct(psi_catt) # same length as psi_catt, so there is no repeat data

papers_to_read_psicatt = as.data.frame(unique(psi_catt$source_article))
papers_to_read_psicatt = cbind(papers_to_read_psicatt,rep("no",times=length(papers_to_read_psicatt)))
names(papers_to_read_psicatt) = c("paper","read?")
write.csv(papers_to_read_psicatt,"C:\\Users\\Jamila\\Documents\\source_papers_psicatt.csv",row.names=FALSE)


################################## A, conyzoides ##################################################
# Completely different method - obtaining the article titles from the EssoilDB raw data
a_cony_info_plant = read.csv("~/a_cony_info_plant.csv")
a_cony_info_compound = read.csv("~/info_compound_15th_may.csv")

# dataframe extracts information that would be found in source articles - compounds present, plant features etc
a_conyzoides = NULL
for(row in 1:nrow(a_cony_info_plant)) {
  code = a_cony_info_plant$CODE[[row]]
  locat = a_cony_info_plant$LOCATION[[row]]
  article = a_cony_info_plant$SOURCE_ARTICLE[[row]]
  samp_date = as.character(a_cony_info_plant$DATE[[row]])
  code_data = a_cony_info_compound[which(a_cony_info_compound$CODE == code),]
  locat_add = rep(locat,times=nrow(code_data))
  art_add = rep(article,times=nrow(code_data))
  date_add = rep(samp_date,times=nrow(code_data))
  code_data = cbind(code_data,locat_add,art_add,date_add)
  a_conyzoides = rbind(a_conyzoides,code_data,make.row.names = FALSE)
}

a_conyzoides = subset(a_conyzoides, select=c("COMPOUND","PERCENTAGE","PLANT_PART","ELUTION_METH","CHEMICAL_FAMILY","locat_add","art_add","date_add"))
fill_NA = rep(NA,times=nrow(a_conyzoides))
a_con_spec = rep("Ageratum conyzoides",times=nrow(a_conyzoides))
a_conyzoides = cbind(a_con_spec,a_conyzoides$COMPOUND,a_conyzoides$PLANT_PART,
                     a_conyzoides$CHEMICAL_FAMILY,a_conyzoides$PERCENTAGE,a_conyzoides$art_add,
                     fill_NA,a_conyzoides$locat_add,fill_NA,fill_NA,fill_NA,a_conyzoides$date_add)

a_conyzoides = as.data.frame(a_conyzoides)
names(a_conyzoides) = c("Species","Chemical","Plant_Part","Chemical_family","percentage","source_article",
                        "Country","Location","Invasivity","wild.not_wild","Habitat","sampling_date")
write.csv(a_conyzoides,"C:\\Users\\Jamila\\Documents\\a_conyzoides.csv")


###################### editing a dataframe from all the source articles (created manually)
# source articles were collected from the species dataframes created in this script and compiled together with the article's DOI

comp_art = read.csv("~/Compiled Source Papers (1).csv") # L. camara, M. quinquenervia and P. cattleainum included (and another species since removed)
comp_art = comp_art[,1:4]
# add A. conyzoides to the dataframe
articles_to_add = unique(a_conyzoides$source_article)
DOI = rep(NA,times=length(articles_to_add))
read_stat = rep("no",times=length(articles_to_add)) # adds a check to see if the article has been read or not to get plant features and raw data (updated manually)
species = rep("Ageratum conyzoides",times=length(articles_to_add))
art_data_to_add = cbind(species,articles_to_add,DOI,read_stat)
art_data_to_add = as.data.frame(art_data_to_add[,1:4])
names(art_data_to_add) = names(comp_art)
comp_art = rbind(comp_art,art_data_to_add,make.row.names=FALSE)
write.csv(comp_art,"C:\\Users\\Jamila\\Documents\\Compiled Source Papers (1).csv")

