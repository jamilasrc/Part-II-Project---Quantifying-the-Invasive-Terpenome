####################### Packages #########################
library(ggplot2)
library(tibble)
library(vegan)
library(dplyr)
library(rio)
library(zip)
library(car)
library(MASS)
library(ggeffects)
library(viridis)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(swamp)
library(dendsort)
library(cluster)
library(ClusterR)
library(factoextra)
library(ggsignif)
library(rstatix)
library(wordspace)
library(ggforce)
library(ggalt)
library(chisq.posthoc.test)
library(proxy)
library(car)
library(sfsmisc)

################## Files #########################

matrix_data_feat = readRDS("matrix_format_features.rds")
matrix_feat_percent = readRDS("matrix_feat_percent.rds")
chemical_coordinates = readRDS("chemical_coordinates.rds")

##################### Functions ################################

# function converts categorical data to numeric - needed for clustering algorithms

cat_to_num = function(mat_in, variable) {
  var_vals = unique(mat_in[[variable]])
  var_numeric = as.numeric(as.factor(var_vals)) # convert variable values to numeric
  names(var_numeric) = var_vals # name variables with corresponding initial name
  
  var_replace = mat_in[[variable]]
  
  for(v in names(var_numeric)) {
    
    if(is.na(v) == FALSE) {
      var_replace[which(var_replace %in% v)] = unname(var_numeric[which(names(var_numeric) == v)]) # replace initial name with new numeric name
      
    }
  }
  mat_in[[variable]] = as.numeric(var_replace) # reconvert numbers to numeric dataform 
  
  return(mat_in)
}

# converts binary plant parts data to numeric factors in 1 column - needed for lm()
plantpart_factor = function(df,cols) {
  # df = input dataframe
  # cols = sequence of column numbers of the data being transformed
  
  factors = seq(1:length(cols))
  names(factors) = cols
  
  return_df = df[,-cols]
  new_col = NULL 
  
  for(row in 1:nrow(df)) {
    part = which(df[row,cols] == 1)
    new_col = append(new_col, unname(factors[which(names(factors) == part)]))
    
  }
  return_df = as.data.frame(cbind(new_col,return_df))
  colnames(return_df)[1] = "Plant_part"
  
  return(return_df)
}


# compare_clust creates clusters of explanatory and predictor variables to see if there is any correlation of sample clustering within these
compare_clust = function(df_feat, df_pred, clust_args, opclus_arg = NULL) { 
  # df_feat = raw dataframe of features/independent variables
  # df_pred = raw dataframe of predicted variable(s)
  
  # clust_args = named list where names are function arguments and values are arg. values for eclust(). 
  # Possible arguments:
  ## FUNcluster = clustering algorithm. Options: "kmeans", "pam", "clara", "fanny", "hclust", "agnes", "diana"
  ## k = no. of clusters. Should be set to null if the user wants the optimal cluster no. to be found, but can be an integer > 1 (possibly < 10)
  ## If FUNcluster = "hclust", "agnes", "diana", hc_metric and hc_method can be specified:
  ### hc_metric = distance method. Options = (Distance) "euclidean", "manhattan", "maximum", "canberra", "binary", "minkowski"; (Similarity) "pearson", "spearman" or "kendall"
  ### hc_method = agglomeration method. Options: Any from hclust() e.g. "pearson", "spearman" or "kendall"
  
  # opclus_args = named list of arguments needed for fviz_nbclust
  # see fviz_nbclust details for potential arguments
  
  set.seed(123) # keep data consistent
  
  out_list = list()
  
  k = unname(clust_args$k) # extrace no. of clusters
  clust_argsin = clust_args # define new cluster arguments that will be modfied in the loop, but allows the interations to rely on initial conditions (from clust_arg)
  
  df_listed_first = list(df_feat,df_pred)
  
  if(is.null(k) == TRUE) { # if k has not been defined, so the shared optimal no. of clusters between the 2 groups needs to be found
    op_clus = NULL
    
    for(i in 1:length(df_listed)) {
      
      df = df_listed_first[[i]]
      opclus_arg[["x"]] = df # adds dataframe to the arguments for fviz_nbclust(), not done in the compare_cluster() arguments so the dataframes don't need to be added to different function argument lists multiple times
      
      opclu_plot = do.call(fviz_nbclust, opclus_arg) # runs fviz_nbclust() with arg. values defined in opclus_args
      op_clu = opclu_plot$data$y
      
      op_clus = rbind(op_clus,op_clu)
      
    }
    find_op = colSums(op_clus)
    no_clu = which.max(find_op) # identifies opt. cluster no. as the cluster with the highest combined av. silhouette for both datasets
    
    
    clust_argsin$k = no_clu # new k now added to function arguments
    
  }
  
  vars = colnames(df_feat)
  var_count = 0
  
  while(var_count < length(vars) + 1) {
    
    if(var_count == 0) {
      df_featvar = df_feat
      
    } else if(length(vars) > 1) {
      df_featvar = subset(df_feat, select = vars[-var_count]) # removes a variable from the features dataframe before the cluster analysis. This acts as a formal of backwards stepwise elimination, identifying which variables are important in influencing clustering
      
    }
    
    df_listed = list(df_featvar,df_pred) # listing allows loop to iterate through dataframes
    
    for(i in 1:length(df_listed)) {
      df = df_listed[[i]]
      clust_argsin[["x"]] = df
      
      mat_cluster = do.call(eclust,clust_argsin) # runs a clustering on df with the (potentially modified) cluster args.
      
      assign(paste0("dend",i),mat_cluster)
      
    }
    clust_dif = cluster.stats(dist(df_pred),dend1$cluster,dend2$cluster) # computes correlations by multiple metrics between clustering (from clusters only, no other info) of the predicted and (expected) predictor variables
    clust_sim = clust_dif$corrected.rand # extracts a correlation metric corrected Rand index (-1/1)
    
    # names correlation according to whether the complete model was tested, or named after the variable that had been removed
    if(var_count == 0) {
      out_list[["full model"]] = clust_sim
      
    } else {
      var = vars[var_count]
      out_list[[var]] = clust_sim
      
    }
    var_count = var_count + 1
    
  }
  
  if(length(vars) == 1) {
    out_list = out_list[[1]] # if df_feat only had one variable, this ensures only the full model is returned. Otherwise rep. values would be returned
    
  } 
  final_list = list(clust_argsin$k,out_list)
  return(final_list)
  
}

# y = (tot. production compounds in cluster n) /( no. compounds in cluster n), per emission profiles
# for invasives only

cluster_production = function(df_feat,df_pred) {
  # assumes df_feat will already have the desired variables
  # and also only the desired species and features etc, use data from above only
  # 0 columns will have also been removed
  
  new_df = df_feat
  
  prod_per_clu = NULL
  
  for(row in 1:nrow(df_pred)) {
    comps = colnames(df_pred)[which(! df_pred[row,] %in% 0)] # remove data where compounds are absent per profile
    comps_prod = unlist(c(df_pred[row,comps])) # creates vector of production rates/numeric values from the columns where compounds were present
    clust = chemical_coordinates$cluster[which(chemical_coordinates$name %in% comps)] # identifies cluster no. of all compounds present 
    
    prod_clusters = NULL
    
    for(clu in c(1,2,3)) { # these are the 3 clusters from k-medoids clustering 
      ind_clu = which(clust %in% clu) 
      
      if(length(ind_clu) > 0) { # if compounds from the cluster are present 
        clu_prod = comps_prod[ind_clu] # gets production rates for compounds present in a profile in the cluster
        prod_for_clu = sum(clu_prod, na.rm = TRUE)/length(ind_clu) # calculates mean production rate per compound in the cluster
        
      } else {
        prod_for_clu = 0 # if a cluster is absent from a profile, set mean production rate as 0
        
      }
      prod_clusters = append(prod_clusters, prod_for_clu) # iterate for each cluster
      
    }
    prod_per_clu = rbind(prod_per_clu,prod_clusters) # iterate for each profile
    
  }
  new_df = cbind(new_df,prod_per_clu) # adds production rates to explanatory variables dataframe so all data is in one place
  new_df= as.data.frame(new_df)
  colnames(new_df)[(ncol(new_df)-2):ncol(new_df)] = c("clu1","clu2","clu3")
  rownames(new_df) = rownames(df_feat)
  
  return(new_df)
  
}



##################################################### Actual Analysis #######################################

# Step 1 - create dataframes per species of feature data and predicted variables, removing samples/variables as appropriate
# Step 2 - create any transfromed data necessary for my arguments
# Step 3 - implement appropriate function on the data sets and report the results

# 2 main analyses:


## Investigating what drives invasive evolution - studying variance within invasive phenotypes.
### This will identify how invasive emission profiles evolve after establishing in a non-invasive environment.
### Will include all data, prioritising including biotic and disturbance data in the model.
### Appropriate Algoriths: compare_clust, cluster_lm, compare_dist
### Species investigated: L. camara, M. quinquenervia (w. and w/o chemotype compounds), P. cattleainum

# I will probably make lists of the relevant data for a loop to iterate through

########### L. camara dataframe construction

# features dataframe
lc_feat = matrix_data_feat[which(matrix_data_feat$Species %in% species[1]),-compound_col] # removes predicted variable data
rownames(lc_feat) = lc_feat$Profile_no
lc_feat = lc_feat[-28,] # remove samples with too much NAs 
## get profile no.s of invasive, non-invasive, native and non-native profiles
lc_invpro = rownames(lc_feat)[which(lc_feat$invasive_status %in% 1)]
lc_ninvpro = rownames(lc_feat)[which(lc_feat$invasive_status %in% 0)]
lc_natpro = rownames(lc_feat)[which(lc_feat$native_status %in% 1)]
lc_nnatpro = rownames(lc_feat)[which(lc_feat$native_status %in% 0)]
##
lc_feat = lc_feat[,-c(1,2,8:10,12,13,16,33,35,37)] # removes variables not of interest
lc_feat = cat_to_num(lc_feat,"iso2c") # converts country code to numeric 
# note - we may need to remove country code - high collinearity with biotic data, which was obtained by country
lc_feat = lc_feat[,-c(6,9,16,17,19)]# remove samples with too many NAs
lc_feat[,c(11,15)] = lapply(lc_feat[,c(11,15)], as.numeric) # convert accidently converted character data back to numeric
lc_feat[2,19] = 1 # adds missing data
lc_feat_noNA = as.data.frame(lc_feat[,-c(6,14,15)]) # final df
lc_feat_facpp = as.data.frame(plantpart_factor(lc_feat_noNA,c(1:5))) # creates factorial plant part needed for cluster_lm


# predicted dataframes
## p/a data
lc_pred_pa = matrix_data_feat[which(matrix_data_feat$Species %in% species[1]),compound_col]
rownames(lc_pred_pa) = matrix_data_feat$Profile_no[which(matrix_data_feat$Species %in% species[1])]
lc_pred_pa = lc_pred_pa[-28,] # removes sample removed from feature data
lc_pred_pa = lc_pred_pa[,which(colSums(lc_pred_pa) != 0)] # remove compounds not present in any emission profiles
## amounts data
lc_pred_am = matrix_feat_percent[which(matrix_feat_percent$Species %in% species[1]),compound_col]
rownames(lc_pred_am) = matrix_feat_percent$Profile_no[which(matrix_feat_percent$Species %in% species[1])]
lc_pred_am = lc_pred_am[-28,] 
lc_pred_am = lc_pred_am[,which(colSums(lc_pred_am) != 0)] 
lc_pred_am = as.matrix(lc_pred_am)
lc_pred_am[which(is.na(lc_pred_am) == TRUE)] = 0
lc_pred_am = as.data.frame(lc_pred_am)

########### analysis 
# clustering comparison
clusterinv_lcam = compare_clust(lc_feat_noNA[which(rownames(lc_feat_noNA) %in% lc_invpro),],
                                lc_pred_am[which(rownames(lc_pred_am) %in% lc_invpro),],
                                clust_args = list(FUNcluster="hclust",k=NULL,hc_metric="euclidean",hc_method="ward.D2"),
                                opclus_arg = list(FUNcluster=get("hcut"),method="silhouette"))



# production per cluster association with explanatory variables 
# only done on L. camara and cluster 1 and 2 due to lack of data from other species and cluster 3

df_prediction = lc_pred_am[which(rownames(lc_pred_am) %in% lc_invpro),]
df_prediction = df_prediction[,which(colSums(df_prediction) != 0)]

lc_clu_prod = cluster_production(df_feat = lc_feat_facpp[which(rownames(lc_feat_facpp) %in% lc_invpro),],df_pred = df_prediction)

# plots to check for association
# cluster 1 
for(col in 1:(ncol(lc_clu_prod)-3)) {
  plot(x=lc_clu_prod[,col],y=lc_clu_prod$clu1,xlab=colnames(lc_clu_prod)[col])
}

for(col in 1:(ncol(lc_clu_prod)-3)) {
  plot(x=lc_clu_prod[,col],y=lc_clu_prod$clu2,xlab=colnames(lc_clu_prod)[col])
}
# data is highly distributed, but some varibles could show correlation
# Run a robust linear model (rlm) on the explantory variables and production rates
# with numeric data only (no factors or binary)
rlm_lc_clu1 = rlm(lc_clu_prod[,-c(1,2,10,12,15:17)],lc_clu_prod$clu1,maxit=40) # cluster 1 production
summary(rlm_lc_clu1)
rlm_lc_clu2 = rlm(lc_clu_prod[,-c(1,2,10,12,15:17)],lc_clu_prod$clu2,maxit=40) # cluster 2 production
summary(rlm_lc_clu2)

# inspect each variable coefficient to see if the variable signficantly predicts cluster production rate
# cluster 1 
sigvar_lc_1 = NULL
pval_11 = NULL
for(var in colnames(lc_clu_prod)[-c(1,2,10,12,15:17)]) {
  f_rob = f.robftest(rlm_lc_clu1, var = var) # this is a wald test, tests if the variables in rlm() make a signficant contribution to the model
  pval_11 = append(pval_11,f_rob$p.value) # gets all p values
  
  if(f_rob$p.value <= 0.1) {
    sigvar_lc_1 = rbind(sigvar_lc_1,c(var,f_rob$p.value)) # creates a list of variables with a p value < 0.1 (signficant and noteable) and the respective significance level
    
  }
}

# repeat process with cluster 2
sigvar_lc_2 = NULL
pval_12 = NULL
for(var in colnames(lc_clu_prod)[-c(1,2,10,12,15:17)]) {
  f_rob = f.robftest(rlm_lc_clu2, var = var)
  pval_12 = append(pval_12,f_rob$p.value) 
  
  if(f_rob$p.value <= 0.1) {
    sigvar_lc_2 = rbind(sigvar_lc_2,c(var,f_rob$p.value))
    
  }
}

# repeat with unique-to-invasive compounds only 
rem_nouniq = which(colSums(lc_pred_am[which(rownames(lc_pred_am) %in% lc_ninvpro),]) == 0) # compounds that are absent from the non-invasive populations (data already cleaned so compounds absent from both populations are removed)
df_prediction2 = lc_pred_am[which(rownames(lc_pred_am) %in% lc_invpro),rem_nouniq]

lc_clu_prod2 = cluster_production(df_feat = lc_feat_facpp[which(rownames(lc_feat_facpp) %in% lc_invpro),],df_pred = df_prediction2)

# plots to check for association
## cluster 1 
for(col in 1:(ncol(lc_clu_prod2)-3)) {
  plot(x=lc_clu_prod2[,col],y=lc_clu_prod2$clu1,xlab=colnames(lc_clu_prod2)[col])
}
## cluster 2
for(col in 1:(ncol(lc_clu_prod2)-3)) {
  plot(x=lc_clu_prod2[,col],y=lc_clu_prod$clu2,xlab=colnames(lc_clu_prod2)[col])
}

# rlm for cluster 1 and 2 
rlm_lc2_clu1 = rlm(lc_clu_prod2[,-c(1,2,10,12,15:17)],lc_clu_prod2$clu1,maxit=20)
summary(rlm_lc2_clu1)
rlm_lc2_clu2 = rlm(lc_clu_prod2[,-c(1,2,10,12,15:17)],lc_clu_prod2$clu2,maxit=20)
summary(rlm_lc2_clu2)

# repeat wald tests and signficance levels for unique compounds
# cluster 1
sigvar_lc2_1 = NULL
pval_21 = NULL
for(var in colnames(lc_clu_prod2)[-c(1,2,10,12,15:17)]) {
  f_rob = f.robftest(rlm_lc2_clu1, var = var) # this is a wald test, tests if the variables in lm() make a signficant contribution to the model
  pval_21 = append(pval_21,f_rob$p.value)
  
  if(f_rob$p.value <= 0.1) {
    sigvar_lc2_1 = rbind(sigvar_lc2_1,c(var,f_rob$p.value))
    
  }
}

# cluster 2
sigvar_lc2_2 = NULL
pval_22 = NULL
for(var in colnames(lc_clu_prod2)[-c(1,2,10,12,15:17)]) {
  f_rob = f.robftest(rlm_lc2_clu2, var = var)
  pval_22 = append(pval_22,f_rob$p.value)
  
  if(f_rob$p.value <= 0.1) {
    sigvar_lc2_2 = rbind(sigvar_lc2_2,c(var,f_rob$p.value))
    
  }
}


# plot model coefficients and signficance levels
## add model coefficients to a dataframe for each model 
rlm_coeff = NULL
rlm_coeff = cbind(rep("Cluster 1 (All)",length(rlm_lc_clu1$coefficients)),names(rlm_lc_clu1$coefficients),unname(rlm_lc_clu1$coefficients))
rlm_coeff = rbind(rlm_coeff,cbind(rep("Cluster 2 (All)",length(rlm_lc_clu2$coefficients)),
                                  names(rlm_lc_clu2$coefficients),unname(rlm_lc_clu2$coefficients)))
rlm_coeff = rbind(rlm_coeff,cbind(rep("Cluster 1 (Unique)",length(rlm_lc2_clu1$coefficients)),
                                  names(rlm_lc2_clu1$coefficients),unname(rlm_lc2_clu1$coefficients)))
rlm_coeff = rbind(rlm_coeff,cbind(rep("Cluster 2 (Unique)",length(rlm_lc2_clu2$coefficients)),
                                  names(rlm_lc2_clu2$coefficients),unname(rlm_lc2_clu2$coefficients)))
rlm_coeff = as.data.frame(rlm_coeff)
colnames(rlm_coeff) = c("RLM","Explanatory Variable","Model Coefficient") # rename variables nicely

# explantory variable names are messily written, so replace with shorter names for plotting
# dataframe colnames = original name, value = replacement
rename_coeff = as.data.frame(matrix(NA,ncol=length(rlm_lc_clu1$coefficients),nrow=1))
colnames(rename_coeff) = names(rlm_lc_clu1$coefficients)
rename_coeff$Preciptation..mm...annual = "Precipitation"
rename_coeff$Humidity..Av..Rel...monthly.or.annum..... = "Relative Humidity"
rename_coeff$Average.high.Temperature..annual.or.month...celcius. = "Temperature (High)"
rename_coeff$Average.low.Temperature..celcius. = "Temperature (Low)"
rename_coeff$Fire.Risk..VIIRS.alerts.2020. = "Wildfire Risk"
rename_coeff$Deforestation..total.tree.cover.loss.over.18yrs. = "Deforestation"
rename_coeff$plantae_20ySR_country = "Plant SR"
rename_coeff$an_20ySR_country = "Animal SR"

# find cells with a given explantory variable in the coefficients dataframe and replace with the new name
for(col in 1:ncol(rename_coeff)) {
  
  if(is.na(rename_coeff[,col]) == FALSE) {
    rlm_coeff$`Explanatory Variable`[which(rlm_coeff$`Explanatory Variable` %in% colnames(rename_coeff)[col])] = rename_coeff[,col]
    
  }
}

# add p value stars as signficance levels
pval_rlm = c(pval_11,pval_12,pval_21,pval_22)
pval_labels_rlm = NULL

for(p in pval_rlm[which(pval_rlm <= 0.05)]) { # only p-values <= 0.05
  
  if(p <= 0.0001) {
    p_star = "****"
  } else if (p <= 0.001) {
    p_star = "***"
  } else if(p <= 0.01) {
    p_star = "**"
  } else if(p <= 0.05) {
    p_star = "*"
  } 
  pval_labels_rlm = append(pval_labels_rlm,p_star)
}
rlm_coeff[,4] = pval_rlm # adding new variable allows for easy plotting of stars
colnames(rlm_coeff)[4] = "p"

rlm_coeff$`Model Coefficient` = as.numeric(rlm_coeff$`Model Coefficient`) # convert coefficients back to numeric

# plots - coefficent points given as lines, starred where signficant
# line at 0 indicates where a coefficient is signficant 
coeff_plot = ggplot(rlm_coeff,aes(y=`Explanatory Variable`,x=`Model Coefficient`,color=RLM)) + geom_point(shape=108,size=10) +
  geom_vline(xintercept = 0,linetype = "dashed") + theme(panel.background =  element_rect(fill = "white",size = 2, linetype = "solid"),
                                                         axis.line = element_line(colour = "black"),legend.title = element_text(size=13), legend.text = element_text(size=11.5), axis.text = element_text(size=13), axis.title = element_text(size=15), plot.title = element_text(size=16))  
# add p_value stars 
coeff_plot = coeff_plot + geom_text(data=rlm_coeff[which(rlm_coeff$p <= 0.05),],aes(x=`Model Coefficient`,y=`Explanatory Variable`),label=pval_labels_rlm,color="black",size=6.5, position = position_nudge(y=0.3))


