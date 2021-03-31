############################### Packages used ##########################################
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
library(usedist)
library(dendextend)
library(flextable)
library(ggpubr)
library(XNomial)
library(emstreeR)

######################################### Files ##############################################################

comp_feat_csv_cleaned = readRDS("comp_feat_csv_cleaned.rds")
matrix_data_feat = readRDS("matrix_format_features.rds")
matrix_data_feat = readRDS("matrix_feat_percent.rds")

######################################## Stats Analysis and Plots ###############################################
species = unique(comp_feat$Species) # identify unique species
compound_col = which(!colnames(matrix_data_feat) %in% colnames(comp_feat_csv_cleaned)) # identify column indices containing compound data

####################### Profile plots (presence and absence)  ############################
# plots are generated as heatmaps to look like a BLAST plot
# colouring indicates presence of a compounds
# profiles are coloured according to whether they are invasive (purple) /non-invasive (green) for all species but A. conyzoides, which are coloured by native and non-native status

colbreak_proplot = c(0,0.99,1.99,2) # split colouring so that values 0 = black, 1 = green and 2 = purple
cols = c("black","#66CC00","#CC66CC")

for(spec in species) {
  
  if(spec == species[4]) {
    var_newcolor = "native_status" # colour according to native status
  } else {
    var_newcolor = "invasive_status" # colour according to invasivity
  }
  submat = matrix_data_feat[which(matrix_data_feat$Species %in% spec),compound_col] # gets compound data for each species
  rownames(submat) = matrix_data_feat$Profile_no[which(matrix_data_feat$Species %in% spec)] # rownames = profile no. for the df
  submat = submat[,which(colSums(submat) != 0)] # removes any absent compounds from the df
  
  # current presence/absence data is binary
  # non-invasive/native no. stays 1, presence of compounds in invasive/non-native profiles is converted to 2 for plotting
  recolour = matrix_data_feat$Profile_no[which(matrix_data_feat[[var_newcolor]] %in% 1 & matrix_data_feat$Species %in% spec)] # identify all invasive/non-native profiles
  for(pro in recolour) {
    row = which(rownames(submat) == pro)
    submat[row,which(submat[row,] %in% 1)] = 2 # set all present compounds to equal 2 
  }
  
  # reorder data so that all invasive/non-native profiles are grouped together and all native/non-invasive profiles are together
  rowinds = which(rownames(submat) %in% recolour) # identifies invasive/non-native profile indices 
  submat = as.data.frame(rbind(submat[rowinds,],submat[-c(rowinds),])) # invasive goes first, non-invasive after
  # plot heatmap
  pheatmap(t(submat),cluster_rows=FALSE,cluster_cols=FALSE,color=cols,breaks=colbreak_proplot,border_color = "black",show_rownames=FALSE,show_colnames=FALSE)
}


############################# Profile Size
profile_sizes = NULL

for(spec in species) {

  ## Profile size by invasive_status
  spec_data = comp_feat_csv_cleaned[which(comp_feat_csv_cleaned$Species %in% spec),] # get species data 
  sampl_spec = unique(spec_data$Profile_no) # get profile no.s 
  
  for(sam in sampl_spec) {
    
    profile_ind = which(spec_data$Profile_no %in% sam) # get row indices for the profile 
    prof_size = length(profile_ind) # get number of compounds per profile = profile size
    inv_stat = as.numeric(as.character((spec_data$invasive_status[[profile_ind[1]]]))) # add profile invasivity
    nat_stat = as.numeric(as.character((spec_data$native_status[[profile_ind[1]]]))) # add profile invasive status
    size_data = c(spec,prof_size,inv_stat,nat_stat)
    profile_sizes = rbind(profile_sizes,size_data)
    
  }
}

profile_sizes = as.data.frame(profile_sizes)
colnames(profile_sizes) = c("Species","Profile_size","invasive_status","native_status")

# reconvert data back to numeric where necessary 
profile_sizes$invasive_status = as.numeric(profile_sizes$invasive_status)
profile_sizes$native_status = as.numeric(profile_sizes$invasive_status)
profile_sizes$Profile_size = as.numeric(profile_sizes$Profile_size)

# A. conyzoides had to be re-done with native and non-native datasets 
profile_size_newac = profile_sizes
# change A. conyzoides invasivity to native status
profile_size_newac$invasive_status[which(profile_size_newac$Species %in% species[4])] = 
  profile_size_newac$native_status[which(profile_size_newac$Species %in% species[4])]

# get number of profiles for each invasive status or native status per species and profile size summary statistics (for initial tables)

for(spec in species) {
  
  for(i in c(0,1)) { # 0 = native/non-invasive, 1 = non-native/invasive
    
    if(spec_i == 4) {
      prosize = profile_size_newac$Profile_size[which(profile_size_newac$Species %in% spec & profile_size_newac$native_status %in% i)] # gets list of profile sizes of a given native status and species
      
    } else {
      prosize = profile_size_newac$Profile_size[which(profile_size_newac$Species %in% spec & profile_size_newac$invasive_status %in% i)] # gets list of profile sizes of a given invasive status and species
      
    }
    print(length(prosize)) # gets no. profiles belonging to each species and invasivity/native status cateogory
    print(summary(prosize)) # gets summary statistics of profile size belonging to each species and invasivity/native status cateogory
    
  }
}


# boxplot of profile sizes for each species, comparing invasive or native status
profile_size_plot = ggplot(profile_size_newac,aes(x=Species,y=Profile_size,fill=invasive_status)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..),hide.ns=TRUE) + ylim(0,100) +
  theme(panel.background =  element_rect(fill = "white",size = 2, linetype = "solid"),axis.line = element_line(colour = "black"),
        legend.title = element_text(size=12), legend.text = element_text(size=12), axis.text = element_text(size=12), axis.title = element_text(size=13)) + 
  xlab("Species") + scale_fill_manual(name = "Invasivity Phenotype",labels = c("Non-Invasive","Invasive"),
                                      values=c("green","red")) + ggtitle("Emission Profiles Sizes of Invasive and Non-invasive Plant Phenotypes") + ylab("Profile Size")

# K-W tests for comparing emission profile size between and within groups
kruskal.test(Profile_size~invasive_status,profile_size_newac)
kruskal.test(Profile_size~Species,profile_size_newac)
for(spec in species) {
  print(kruskal.test(Profile_size~invasive_status,profile_size_newac[which(profile_size_newac$Species %in% spec),]))
}


############## Bray-Curtis Dissimilarity 
bray_curtis = vegdist(subset(matrix_data_feat, select = c(compounds)),method="bray") # calculates pairwise BC creates distance matrix
bray_curtis2 = as.matrix(unlist(as.matrix(bray_curtis))) # finalise symmetrical distance matrix

bray_compare = NULL

for(spec in species) {
  
  pros = matrix_data_feat$Profile_no[which(matrix_data_feat$Species %in% spec)] # gets profile no. of the species
  inv_stat = matrix_data_feat$invasive_status[which(matrix_data_feat$Species %in% spec)] # gets all invasive status
  names(inv_stat) = as.character(pros) # assigns profile no. to invasivity
  
  for(pro1 in pros) {
    
    for(pro2 in pros) {
      
      if(pro1 < pro2) { # makes sure repeat comparisons aren't added to the data table
        
        dist_col = which(colnames(bray_curtis2) == pro1) # gets column number of profile 1 
        dist_row = which(rownames(bray_curtis2) == pro2) # gets row number of profile 2
        bc_dist = bray_curtis2[dist_row,dist_col] # gets pairwise BC dissimiarlity for the 2 profiles
        
        # add invasivity and invasivity combination in the pairwise comparison
        inv1 = unname(inv_stat[which(names(inv_stat) == as.character(pro1))]) # gets invasivity of profile 1
        inv2 = unname(inv_stat[which(names(inv_stat) == as.character(pro2))])
        
        # invasivity combo
        if(inv1 == inv2) { # produces 0/0 or 1/1
          inv_comb = paste0(inv1,"/",inv2)
          
        } else { # makes sure non-invasive/invasive and invasive/non-invasive combinations have the same labelling
          inv_comb = "0/1"
        }
        
        data_add = c(pro1,pro2,spec,spec,inv1,inv2,inv_comb,bc_dist)
        bray_compare = rbind(bray_compare,data_add)
        
      }
    }
  }
}

bray_compare = as.data.frame(bray_compare)
colnames(bray_compare) = c("Profile_A","Profile_B","Species_A","Species_B","Inv_A","Inv_B","inv_comb","dissimilarity")

# convert back to numeric data 
bray_compare$Profile_A = as.numeric(bray_compare$Profile_A)
bray_compare$Profile_B = as.numeric(bray_compare$Profile_B)

# convert A. conyzoides to proper values
bray_compare_newac = ac_data_converter(bray_compare,"Profile_A","Profile_B","Species_A","Inv_A","Inv_B","inv_comb")

# boxplot of pairwise comparisons for each species
ggplot(bray_compare_newac,aes(x=Species_A,y=dissimilarity,fill=inv_comb)) + geom_boxplot() + stat_compare_means(aes(label = ..p.signif..),hide.ns=TRUE) + ylim(0,1) + 
  theme(panel.background =  element_rect(fill = "white",size = 2, linetype = "solid"),axis.line = element_line(colour = "black"),
        legend.title = element_text(size=14), legend.text = element_text(size=14), axis.text = element_text(size=15), 
        axis.title = element_text(size=15),axis.text.x = element_text(face="italic"), plot.title = element_text(size=16)) + 
  xlab("Species") + ylab("Dissimilarity") + scale_fill_manual(name = "",labels = c("Non-Invasive/Non-Invasive","Non-Invasive/Invasive","Invasive/Invasive"),
                                                              values=c("#66CC00","yellow","#CC66CC")) + ggtitle("Pairwise Bray-Curtis Dissimilarity Between Emission Profiles")

# averages and stats test
for(spec in species) {
  print(median(bray_compare_newac$dissimilarity[which(bray_compare_newac$Species_A %in% spec)])) # get median B-C dissimiarity per species
  print(kruskal.test(dissimilarity~inv_comb,bray_compare_newac[which(bray_compare_newac$Species_A %in% spec),])) # Kruskal-Wallis test on pairwise dissimilarity between invasivity groups, done for each species
  print(dunn_test(data=bray_compare_newac,formula=dissimilarity~inv_comb)) # post-hoc for Krusk-Wallis
  
}


###################################################################### Clustering of Profiles ############################

# Overall rationale 
## all clustering will be split into 2 groups/clusters, as this is predicted by our hypothesis that profiles 
##  should segregate into invasive and non-invasive groups, or native and non-native groups. 
## The true clustering into 2 groups is validated against expected clustering

################################################## Functions #############################

# splits dataset into quantiles
quantile_breaks = function(df, n = 10) {
  breaks = quantile(df, probs = seq(0, 1, length.out = n),na.rm=TRUE)
  breaks = breaks[!duplicated(breaks)]
  return(breaks)
}

# creates breaks for a colour scale where data is split into 0 + 10% quanitles (removing 0 values in case the other quantiles = 0)
quantile_colours = function(matrix) {
  mat_col_scale = matrix
  matrix_col_vec = NULL
  
  for(r in 1:nrow(mat_col_scale)) {
    
    for(c in 1:ncol(mat_col_scale)) {
      
      if(mat_col_scale[r,c] == 0) {
        mat_col_scale[r,c] = NA
        
      } 
      matrix_col_vec = append(matrix_col_vec,mat_col_scale[r,c])
      
    }
  }
  colour_scale = quantile_breaks(matrix_col_vec, n = 11)
  final_colour = c(0,colour_scale)
  return(final_colour)
}



# pheatmap plotter

pheatmap_cluster = function(mat,df,group_split=NULL,annotate_vars=NULL,annotate_rename=NULL,new_colour=NULL,col_pal,datatype=NULL) { 
  # mat = matrix with rows - variables to be clustered and columns - samples
  # df = original full dataset of features (from which mat was derived)
  # group_split = named list of variable name and variable values
  # annotate_vars = max length 2, variables to annotate cluster w. on a featmap
  # col_pal = call for a colour palette from a various R package
  
  # mat = matrix_clu_per2
  #df = matrix_feat_percent
  # group_split = list(species)
  # names(group_split) = "Species"
  # annotate_vars = "invasive_status"
  # new_colour = "Y"
  # datatype = "percentage amount"
  #col_pal = RColorBrewer::brewer.pal(11, "RdYlBu")
  
  #annotate_vars = "invasive_status"
  
  
  
  #ann_col = NA
  #ann_row = NA
  
  if(is.null(annotate_vars)==FALSE) {
    annotator = as.data.frame(matrix(NA,ncol=length(annotate_vars),nrow=ncol(mat)))
    for(ind in 1:length(annotate_vars)) {
      
      anno_var = annotate_vars[ind]
      anno = as.data.frame(unlist(lapply(df[[anno_var]], as.character)))
      annotator[,ind] = anno
    }
    rownames(annotator) = colnames(mat)
    if(is.null(annotate_rename)==FALSE) {
      colnames(annotator) = annotate_rename
    } else {
      colnames(annotator) = annotate_vars
    }
  } 
  
  if(is.null(group_split)==FALSE) {
    group_list = list()
    var = names(group_split)
    group_split2 = unlist(group_split)
    
    for(gr_val in group_split2) {
      
      
      rows_gr = which(df[[var]] %in% gr_val)
      mat_sub = mat[,rows_gr]
      mat_sub = mat_sub[rowSums(mat_sub!=0)>0,]
      
      anno_gr = subset(annotator,rownames(annotator) %in% colnames(mat_sub))
      anno_gr = as.data.frame(anno_gr)
      
      
      colnames(anno_gr) = colnames(annotator)
      rownames(anno_gr) = colnames(mat_sub)
      
      
      
      if(is.null(new_colour)==FALSE) {
        final_col = quantile_colours(mat_sub)
      } else if(datatype == "binary") {
        final_col = c(0,0.99,1)
      } else {
        final_col = NA
      }
      
      mat_cluster_cols = hclust(dist(t(mat_sub)))
      mat_cluster_cols = sort_hclust(mat_cluster_cols)
      #mat_cluster_rows = sort_hclust(hclust(dist(mat_sub)))
      
      pheatmap = pheatmap(mat_sub,show_rownames=FALSE,show_colnames=FALSE,annotation_col=anno_gr,
                          scale = "none",clustering_method="ward.D2",
                          clustering_distance_cols="euclidean",color=col_pal,breaks=final_col,border_color = "black",
                          cluster_rows=FALSE,cluster_cols=mat_cluster_cols,main=bquote(italic(.(gr_val))),annotation_names_col=FALSE,fontsize=13)
      
      list_ind = length(group_list)+1
      group_list[[list_ind]] = list(mat_sub,pheatmap)
      
    }
    names(group_list) = group_split
    return(group_list)
    
  } else {
    
    
    if(is.null(new_colour)==FALSE) {
      final_col = quantile_colours(mat)
    } else if(datatype == "binary") {
      final_col = c(0,0.99,1)
    } else {
      final_col = NA
      
    }
    print(final_col)
    
    mat_cluster_cols = hclust(dist(t(mat)))
    mat_cluster_cols = sort_hclust(mat_cluster_cols)
    # mat_cluster_rows = sort_hclust(hclust(dist(mat)))
    
    pheatmap = pheatmap(mat,show_rownames=FALSE,show_colnames=TRUE,annotation_col=annotator,
                        scale = "none",clustering_method="ward.D2", 
                        clustering_distance_cols="euclidean",color=col_pal,breaks=final_col,border_color = "black",
                        cluster_rows=FALSE,cluster_cols=mat_cluster_cols,main="Hierarchal Cluster Pheatmap",annotation_names_col=FALSE,fontsize=13)
    
    return(pheatmap) 
    
  }
}

# Cluster validation function
cluster_valid = function(mat,df=NULL,group_split=NULL,cluster_no,expect_cluster) {
  # mat,df and group_split are the same as in pheatmap plotter
  # cluster_no = no. of clusters that should be clustered into
  # expected cluster needs to be numeric and in the same order as the matrix or matrix split (I think change binary to 1,2)
  
  
  if(is.null(group_split)==FALSE) {
    group_list = list()
    var = names(group_split)
    group_split2 = unlist(group_split)
    
    for(gr_val in group_split2) {
      
      rows_gr = which(df[[var]] %in% gr_val)
      mat_sub = mat[rows_gr,]
      mat_sub = mat_sub[,which(colSums(mat_sub) != 0)]
      exp_clu_sub = expect_cluster[rows_gr]
      
      mat_cluster = eclust(mat_sub,"hclust",k=cluster_no,hc_metric = "euclidean", hc_method = "ward.D2")
      clust_dif = cluster.stats(dist(mat_sub),exp_clu_sub,mat_cluster$cluster) # calculates how similar the expected cluster and cluster by the algorithm is 
      clust_sim = clust_dif$corrected.rand # varies between -1 (no similarity) to 1 (complete similarity)
      
      list_ind = length(group_list)+1
      group_list[[list_ind]] = clust_sim
      
    }
    names(group_list) = group_split2
    
    return(group_list)
    
  } else {
    mat_cluster = eclust(mat,"hclust",k=cluster_no,hc_metric = "euclidean", hc_method = "ward.D2")
    clust_dif = cluster.stats(dist(mat),expect_cluster,mat_cluster$cluster) 
    clust_sim = clust_dif$corrected.rand 
    
    return(clust_sim)
    
  }
}


##### Clustering of Profiles from presence/absence compound data

# create dendograms
for(spec in species) {
  
  
  submat = matrix_data_feat[which(matrix_data_feat$Species %in% spec),compound_col] # matrix of all profiles of the specified species, data = compounds only
  rownames(submat) = matrix_data_feat$Profile_no[which(matrix_data_feat$Species %in% spec)] # set rownames as profile number
  submat = submat[,which(colSums(submat) != 0)] # remove any compounds that are not present in any chemical profile of this subset
  
  hc_cut = eclust(submat,"hclust",k=2,hc_metric = "euclidean", hc_method = "ward.D2") # cluster matrix by hierachical clustering using Ward Clustering and Euclidean distances, cut into 2 groups  
  
  if(spec == species[4]) {
    annotate = matrix_data_feat$native_status[which(matrix_data_feat$Species %in% spec)] # for A. conyzoides, native and non-native groups are annotated on the dendrogram
    
  } else {
    annotate = matrix_data_feat$invasive_status[which(matrix_data_feat$Species %in% spec)] # for all other species invasive and non-invasve groups are annotated
    
  }
  # convert annotated groups, which can take binary values 0 and 1, to 1 and 2 so that it can be directly compared with cluster number
  #   note this means if the first cluster is set as 1, but mainly contains invasive/non-native profiles, this will result in -ve clustering
  #   this is ok as we only have 2 cluster groups, but note the modulus of the clustering should be taken into account, rather than whether it is +ve/-ve
  ind1 = which(annotate %in% 1) 
  ind0 = which(annotate %in% 0)
  annotate[ind1] = 1
  annotate[ind0] = 2
  
  # create vector of annotation colours 
  anno_cols = vector(length=length(annotate))
  anno_cols[ind1] = "#CC66CC"  # invasive/non-native = purple
  anno_cols[ind0] = "#66CC00"  # non-invasive/native = green
  
  hc_plot = as.dendrogram(hc_cut) # produce dendrogram object from hierachal clustering
  hc_plot = colour_branches(hc_plot,k=2) # colour branches according to true clustering into 2 groups
  labels_colors(hc_plot) = "white" # removes profile no. labels
  plot(hc_plot, main = bquote(italic(.(spec))), cex.main = 1.5) # plots the dendrogram 
  colored_bars(anno_cols,hc_plot,rowLabels = "") # adds annotation bar to the plot 
  
 
}

# validate clusters 
group_splitter = list(species)
names(group_splitter) = "Species"

expect_cluster_all = matrix_feat_percent$invasive_status
expect_cluster_all[which(matrix_feat_percent$invasive_status %in% 0)] = 1
expect_cluster_all[which(matrix_feat_percent$invasive_status %in% 1)] = 2

bin_cluster = matrix_data_feat[,compound_col]
bin_cluster = as.data.frame(bin_cluster)
rownames(bin_cluster) = matrix_data_feat$Profile_no

pre_abs_cluster = cluster_valid(mat=data.matrix(bin_cluster),df=matrix_data_feat,group_split=group_splitter,
                                cluster_no=2,expect_cluster=expect_cluster_all) 

# with native/non-native A. conyzoides
expect_ac = matrix_data_feat$native_status[which(matrix_data_feat$Species %in% species[4])]
expect_ac_func = expect_ac
expect_ac_func[which(expect_ac %in% 0)] =  1
expect_ac_func[which(expect_ac %in% 1)] =  2

ac_preabs_cluster = cluster_valid(mat=data.matrix(bin_cluster[which(matrix_data_feat$Species %in% species[4]),]),
                                  cluster_no=2,expect_cluster=expect_ac_func)

########################### clustering with amounts 

matrix_am = matrix_feat_percent[,compound_col]
rownames(matrix_am) = matrix_feat_percent$Profile_no
matrix_am = as.data.frame(t(matrix_am))

for(c in 1:ncol(matrix_am)) {
  for(r in 1:nrow(matrix_am)) {
    
    if(is.na(matrix_am[r,c]) == TRUE) {
      matrix_am[r,c] = 0
    }
  }
}


 
# plot pheatmaps

final_plots = pheatmap_cluster(mat=matrix_am,df=matrix_feat_percent,
                               annotate_vars="invasive_status",annotate_rename="Invasivity",group_split=group_splitter,
                               new_colour="Y",col_pal=c("black",RColorBrewer::brewer.pal(10, "RdYlBu")))

# remove chemotype-biasing compounds from M. quinquenervia
# repeat cluster analysis for melaleuca quinquenervia removing chemotype
# 1,8-cineole, viridiflorol and E-nerolidol need to be removed
rem_comp = c("1,8-cineole","viridiflorol","E-nerolidol")
SMI_rem = NULL
for(comp in rem_comp) {
  SMI = comp_feat_csv_cleaned$SMILES[which(comp_feat_csv_cleaned$Compound %in% comp)[[1]]]
  SMI_rem = append(SMI_rem,SMI)
}


data_melrem = pheatmap_cluster(mat=subset(matrix_am[-c(which(rownames(matrix_am) %in% SMI_rem)),]),df=matrix_feat_percent,
                               annotate_vars="invasive_status",annotate_rename="Invasivity",group_split=group_splitter,
                               new_colour="Y",col_pal=c("black",RColorBrewer::brewer.pal(10, "RdYlBu")))

# A. conyzoides native/non-native 
data_ac_nat = pheatmap_cluster(mat=matrix_am,df=matrix_feat_percent,
                               annotate_vars="native_status",annotate_rename="Native or Non-Native",group_split=group_splitter,
                               new_colour="Y",col_pal=c("black",RColorBrewer::brewer.pal(10, "RdYlBu")))


# validate clusters 

# remove chemotype-biasing compounds from M. quinquenervia
mat_rem = matrix_am
mat_rem = subset(mat_rem, select = (colnames(mat_rem) %in% matrix_feat_percent$Profile_no[which(matrix_feat_percent$Species %in% species[2])]))
mat_rem = mat_rem[-c(which(rownames(mat_rem) %in% SMI_rem)),]
mat_rem = as.data.frame(t(mat_rem))

expect_mq = matrix_data_feat$invasive_status[which(matrix_data_feat$Species %in% species[2])]
expect_mq_func = expect_mq
expect_mq_func[which(expect_mq %in% 0)] =  1
expect_mq_func[which(expect_mq %in% 1)] =  2


am_clu_mel = cluster_valid(mat=data.matrix(mat_rem),cluster_no=2,expect_cluster=expect_mq_func)


# A. conyzoides native/non-native 
mat_acam = matrix_am[,which(matrix_feat_percent$Species %in% species[4])]
mat_acam = as.data.frame(t(mat_acam))
mat_acam = mat_acam[,which(colSums(mat_acam) != 0)]

ac_am_cluster = cluster_valid(mat=mat_acam,cluster_no=2,expect_cluster=expect_ac_func)

###################### chemical similarity clustering 

dist_between_centroids(fing_diss_matx,colnames(fing_diss_matx)[2:5],rownames(fing_diss_matx)[2:5])

# calculate centroid distance from tanimoto's distance matrix of compounds, where each profile = group of compounds belonging to it

dist_profiles_chemsim = function(df,distmat,profilesub) {
  # df = dataframe of each profile with p/a data
  # distmat = distance matrix containing items that will be in the profile
  # profilesub = list of profile numbers for which the average distance is going to be calculated
  
  out_distmat = as.data.frame(matrix(NA,nrow=length(profilesub),ncol=length(profilesub))) # initialise the outputted distance matrix
  colnames(out_distmat) = profilesub
  rownames(out_distmat) = profilesub
  
  for(pro1 in profilesub) {
    
    for(pro2 in profilesub) {
      
      cell_val = out_distmat[which(rownames(out_distmat) == pro2),which(colnames(out_distmat) == pro1)]
      
      if(is.na(cell_val) == TRUE) {
        right_row1 = which(rownames(df) == pro1)
        comps_1 = colnames(df)[which(df[right_row1,] %in% 1)]
        right_row2 = which(rownames(df) == pro2)
        comps_2 = colnames(df)[which(df[right_row2,] %in% 1)]
        
        dist_bet_group = dist_between_centroids(distmat,comps_1,comps_2)
        
        out_distmat[which(rownames(out_distmat) == pro2),which(colnames(out_distmat) == pro1)] = dist_bet_group
        out_distmat[which(rownames(out_distmat) == pro1),which(colnames(out_distmat) == pro2)] = dist_bet_group
        
      }
    }
  }
  return(out_distmat)
}

for(spec in species) {
  
  spec_distcomp = dist_profiles_chemsim(as.data.frame(matrix_data_feat[,compound_col],row.names=matrix_data_feat$Profile_no),
                                        fing_diss_matx,matrix_data_feat$Profile_no[which(matrix_data_feat$Species %in% spec)])
  
  for(i in 1:nrow(spec_distcomp)) {
    spec_distcomp[i,i] = 0
  }
  
  spec_distcomp = as.dist(spec_distcomp)
  hc_clust = hclust(spec_distcomp)
  hc_cut = hcut(spec_distcomp,k=2)
  
  if(spec == species[4]) {
    annotate = matrix_data_feat$native_status[which(matrix_data_feat$Species %in% spec)]
    
  } else {
    annotate = matrix_data_feat$invasive_status[which(matrix_data_feat$Species %in% spec)]
    
  }
  ind1 = which(annotate %in% 1)
  ind0 = which(annotate %in% 0)
  annotate[ind1] = 1
  annotate[ind0] = 2
  
  anno_cols = vector(length=length(annotate))
  anno_cols[ind1] = "#CC66CC" 
  anno_cols[ind0] = "#66CC00"
  
  hc_plot = as.dendrogram(hc_clust)
  hc_plot = colour_branches(hc_plot,k=2)
  labels_colors(hc_plot) = "white"
  plot(hc_plot, main = bquote(italic(.(spec))), cex.main = 1.5)
  
  if(spec == species[4]) {
    colored_bars(anno_cols,hc_plot,rowLabels = "")
  } else {
    colored_bars(anno_cols,hc_plot,rowLabels = "")
  }
  
  
  clusters = hc_cut$cluster
  clu_stat = cluster.stats(spec_distcomp,clustering = clusters, alt.clustering = annotate)
  print(clu_stat$corrected.rand)
  
}

################################## Compound Clustering ##############################################################


# Identify optimal cluster number 
Optimal_Clusters_Medoids(fing_diss_matx,max_clusters = 20,distance_metric = "euclidean",criterion = "silhouette")
# 3 is optimal 


chem_clust = pam(fing_diss_matx,k=3) # cluster compounds according to k-medoids clustering w. optim. cluster no.  
# plot clustering 
clust_plot = fviz_cluster(chem_clust,geom="point",ggtheme = theme_bw(),main="Hierarchal Cluster Plot of Chemical Distances of Compounds Present in the Emission Profiles")
chemical_coordinates = data.frame(clust_plot$data[c("name","x","y","cluster")],row.names=NULL) #extract coordinates generated from the distance matrix for future plotting 

##################################### MST plots for clustering ########################################################
# these plots represent total diversity per population (invasive or non-invasive, non-native or native)
# plotted from coordinates from clustering and compounds are coloured according to clustering 
# edges = min. distance network, but generally of not much signficance, just to plot a network

########### functions
MST_clusters = function(nodes_df,data_df,split_var,var_rename=NULL,clust_colours) {
  # nodes_df = df of coordinates that the MST will be generated from.
  #  This should contain the clustering of each node, and the x and y coordinates should be in the first 2 columns
  # data_df = df of data containg information on which nodes should be selected for clustering. Named x, y, cluster
  #  e.g. only 1 species if 1 species should be clustered
  # split_var = variable the different comparative trees will be generated from
  # var_rename = changes title according to variable value, vector named w. real variable values and values = renaming
  # clust_colurs = named vector, names = as.character(clust no.), values = desired colours
  
  out_list = list()
  
  split_vals = unique(data_df[[split_var]]) # gets values of varible data is being split by (here will be invasive_status or native_status)
  
  for(val in split_vals) {
    
    small_data = data_df[which(data_df[[split_var]] %in% val),] # get data for particular variable value
    small_data = small_data[,compound_col] # get compound data only
    comp_keep = unname(which(colSums(small_data) != 0)) # get indices of which compounds are present for the variable value
    
    comp_df = nodes_df[comp_keep,]
    
    print(c(val,nrow(comp_df))) # returns number of compounds present in the variable value group (for tables)
    
    MST = ComputeMST(comp_df[,1:2]) # create a minimum spanning tree from the coordinates from initial clustering, only for compounds present for the variable value group
    # add cluster number to the data for plotting
    MST = as.data.frame(cbind(MST,comp_df[,3])) 
    colnames(MST)[ncol(MST)] = "cluster"
    
    colouring = unname(clust_colours[which(names(clust_colours) %in% as.character(unique(MST$cluster)))]) # get the cluster colours for for clusters present in the variable value group
    
    # title the plot with the grouping variable value
    if(is.null(var_rename) == FALSE) {
      title = unname(var_rename[which(names(var_rename) == as.character(val))])
    } else {
      title = as.character(var)
    }
    
   # generate network plot
     plot = ggplot(MST, aes(x = x, y = y, from = from, to = to, colour = cluster)) + geom_point(size=4) + 
      stat_MST(linetype="dashed",colour="black") + theme_void() + ggtitle(title) + scale_colour_manual(values=colouring) + xlim(-15,20) + ylim(-20,15)
    
    out_list[[length(out_list)+1]] = plot # add plot to a list of plots for all variable values
    
    
  }
  return(out_list)
  
}

# create plots for each species, split/grouped according to invasive_status
Lc_MST = MST_clusters(nodes_df = chemical_coordinates[,2:4],data_df = matrix_data_feat[which(matrix_data_feat$Species %in% species[1]),],
                      split_var = "invasive_status",var_rename = c("1" = "Invasive", "0" = "Non-Invasive"),
                      clust_colours = c("1" = "#FF3399", "2" = "#00CCFF", "3" = "#9966FF"))
Lc_MST[[1]]
Lc_MST[[2]]

Mq_MST = MST_clusters(nodes_df = chemical_coordinates[,2:4],data_df = matrix_data_feat[which(matrix_data_feat$Species %in% species[2]),],
                      split_var = "invasive_status",var_rename = c("1" = "Invasive", "0" = "Non-Invasive"),
                      clust_colours = c("1" = "#FF3399", "2" = "#00CCFF", "3" = "#9966FF"))
Mq_MST[[1]]
Mq_MST[[2]]

Pc_MST = MST_clusters(nodes_df = chemical_coordinates[,2:4],data_df = matrix_data_feat[which(matrix_data_feat$Species %in% species[3]),],
                      split_var = "invasive_status",var_rename = c("1" = "Invasive", "0" = "Non-Invasive"),
                      clust_colours = c("1" = "#FF3399", "2" = "#00CCFF", "3" = "#9966FF"))
Pc_MST[[1]]
Pc_MST[[2]]

# grouped according to native_status for A. conyzoides
Ac_MST = MST_clusters(nodes_df = chemical_coordinates[,2:4],data_df = matrix_data_feat[which(matrix_data_feat$Species %in% species[4]),],
                      split_var = "native_status",var_rename = c("1" = "Non-Native", "0" = "Native"),
                      clust_colours = c("1" = "#FF3399", "2" = "#00CCFF", "3" = "#9966FF"))

# plots will be grouped together and finished off in word 

################################    Compound Production Levels from Clustering  ######################################

#### functions ##########

# extract legend from a ggplot for plotting ggplots together in grid.arrange()
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#####################

# initialise data matrix containing chemical coordiantes production rates (0s included and removed), 
#  p-values where there is a signficant difference between production in invasive and non-invasive populations
# chemical coordiantes from clustering and compound SMILES are included immediately
cords_spec_ac = as.data.frame(cbind(rep(NA,length=nrow(chemical_coordinates)),
                                    chemical_coordinates,rep(NA,length=nrow(chemical_coordinates)),rep(NA,length=nrow(chemical_coordinates))))

# add species list, where rep. of species = length of all compounds
species_add = NULL
for(spec in species) {
  sub_add = rep(spec,length=nrow(chemical_coordinates))
  species_add = append(species_add,sub_add)
}

# duplicate data matrix for all species
cords_spec_ac = do.call(rbind, replicate(4, cords_spec_ac, simplify=FALSE))
colnames(cords_spec_ac)[c(1,6,7)] = c("Species","dif_pro_produce","dif_amount")

cords_spec_ac[,1] = species_add # add repeated species names to the matrix 
cords_spec_ac = as.data.frame(cbind(cords_spec_ac,rep(NA,length=nrow(cords_spec_ac)),rep(NA,length=nrow(cords_spec_ac))))
colnames(cords_spec_ac)[8:9] = c("p","mean_dif_no0") # add columns for p-values and the difference in production when 0s are removed from invasive and non-invasive populations (for plotting purposes only)


for(spec in species) {
  
  matrix_spec = matrix_feat_percent[which(matrix_feat_percent$Species %in% spec),] # intialise matrix with data for one species only 
  matrix_spec = matrix_spec[,-c(which(colSums(matrix_spec[,compound_col] != 0) == 0))] # removes absent compounds from the species dataframe
  spec_comp = unique(comp_feat_csv_cleaned$SMILES[which(comp_feat_csv_cleaned$Species %in% spec)]) # identifies distinct compounds found in the species
  
  for(comp in spec_comp) { # for each compound
    comp_ind_spec = NULL
    comp_ind_spec = which(cords_spec_ac$name %in% comp & cords_spec_ac$Species %in% spec) # identifies index for the compound and species
    
    if(is.null(comp_ind_spec) == FALSE) {
      inv_mean = NULL
      ninv_mean = NULL
      
      if(spec == species[4]) { # uses native status 
        inv_vec = c(matrix_spec[[comp]][which(matrix_spec$native_status %in% 1)]) # get production rates of the compound in non-native profiles 
        ninv_vec = c(matrix_spec[[comp]][which(matrix_spec$native_status %in% 0)]) # get production rates of the compound in native profiles 
        
      } else { # uses invasivity 
        inv_vec = c(matrix_spec[[comp]][which(matrix_spec$invasive_status %in% 1)]) # get production rates of the compound in invasive profiles 
        ninv_vec = c(matrix_spec[[comp]][which(matrix_spec$invasive_status %in% 0)]) # get production rates of the compound in non-invasive profiles 
        
      }
      # remove NA values
      inv_vec = inv_vec[!is.na(inv_vec)] 
      ninv_vec = ninv_vec[!is.na(ninv_vec)] 
      
      # remove 0 values from production rates (for plotting purposes)
      inv_vec_n0 = inv_vec[which(inv_vec != 0)]
      ninv_vec_n0 = ninv_vec[which(ninv_vec != 0)]
      
      dif = NULL # mean production rate of the compound in invasive/non-native populations (no 0 values)
      dif2 = NULL # mean production rate of the compound in non-invasive/native populations (no 0 values)
      dif_w0 = NULL # mean production rate of the compound in invasive/non-native populations (with 0 values)
      dif2_w0 = NULL # mean production rate of the compound in non-invasive/native populations (with 0 values)
      MW1 = NULL # is it possible to calculate mann-whitney test considering invasive/non-native data? ("Y" or NULL)
      MW2 = NULL # is it possible to calculate mann-whitney test considering non-invasive/native data? ("Y" or NULL)
      
      # invasive/non-native profiles
      if(length(inv_vec_n0) > 1) { # if the compound is produced in more than 1 invasive/non-native profile - can calculate M-W test
        dif = mean(inv_vec_n0,na.rm=TRUE) # mean production rate (no 0 values/where it is produced)
        dif_w0 = mean(inv_vec,na.rm=TRUE) # mean production rate (with 0 values)
        MW1 = "Y"
        
      } else if(length(inv_vec_n0) == 1) { # if the compound is produced in only 1 invasive profile - cannot calculate M-W test
        dif = inv_vec_n0
        dif_w0 = mean(inv_vec,na.rm=TRUE)
        
      } else if(length(inv_vec_n0) == 0) { # if the compound is absent from invasve profiles - add "only in non-invasive" to all columns
        subdata = "only non-invasive"
        cords_spec_ac[comp_ind_spec,6:9] = rep(subdata,length=4)
      }
      
      # repeat process for non-invasive/native profiles 
      if(length(ninv_vec_n0) > 1) {
        dif2 = mean(ninv_vec_n0,na.rm=TRUE)
        dif2_w0 = mean(ninv_vec,na.rm=TRUE)
        MW2 = "Y"
        
      } else if(length(ninv_vec_n0) == 1) {
        dif2 = ninv_vec_n0
        dif2_w0 = mean(ninv_vec,na.rm=TRUE)
        
      } else if(length(ninv_vec_n0) == 0) {
        subdata = "only invasive"
        cords_spec_ac[comp_ind_spec,6:9] = rep(subdata,length=4)
      }
      
      # calculating differences in production rate between invasive/non-native and non-invasive/native populations
      # calculated proportionally/in terms of amplification 
      
      if(is.null(dif)==FALSE & is.null(dif2)==FALSE) { # if the compound was present in both populations
        dif_n0 = dif/dif2 # difference in production (no 0s) - plotting
        dif_w0 = dif_w0/dif2_w0 # difference in production (with 0s) - stats tests
        
        cords_spec_ac$dif_amount[[comp_ind_spec]] = dif_w0
        cords_spec_ac$mean_dif_no0[[comp_ind_spec]] = dif_n0
        
      }
      
      if(length(inv_vec_n0) > 0 & length(ninv_vec_n0) > 0) { # if the compound was present in both populations
        # calculates the difference between the proportion of profiles producing the compound in invasive/nnon-native populations vs non-invasive/native populations
        #  for plotting
        dif_samp = (length(inv_vec_n0)/length(inv_vec))/(length(ninv_vec_n0)/length(ninv_vec))  
        cords_spec_ac$dif_pro_produce[[comp_ind_spec]] = dif_samp
        
      }
      
      if(is.null(MW1) == FALSE & is.null(MW2) == FALSE) { # if a M-W test can be calculated
        MW = wilcox.test(inv_vec,ninv_vec) # M-W test between the 2 populations
        cords_spec_ac$p[[comp_ind_spec]] = MW$p.value # add significance level to the dataframe 
        
      }
    }
  }
}

cords_spec_ac = cords_spec_ac[!is.na(cords_spec_ac$dif_pro_produce),] # removes compound that were absent from a species
cords_spec_1ac = cords_spec_ac

cords_spec_1ac[,5:9] = lapply(cords_spec_1ac[,5:9], as.numeric) # reconverts data back to numeric, converting any compounds only present in one population to NA
cords_spec_1ac = cords_spec_1ac[!is.na(cords_spec_1ac$dif_pro_produce),] # removes compounds only present in 1 population (by removing NAs)


# create colour and size scales from data for all species - this means the colouring and legend labelling will be consistent between species plots

# create a colour scale for the plot based on production levels in profiles producing the compound
# based on quantile breaks in production rates
cords_spec1ac_colscale = c(quantile_breaks(cords_spec_1ac$mean_dif_no0,n=11))
colour_labels_ac = NULL
for(col in 2:length(unname(cords_spec1ac_colscale))) {
  col_1 = signif(unname(cords_spec1ac_colscale)[col-1],digits=3) # production breaks rounded to 3.s.f
  col_2  = signif(unname(cords_spec1ac_colscale)[col],digits=3)
  
  if(col == 2) {
    colname = paste0("[",as.character(col_1),",",as.character(col_2),"]") # create legend labelling
  } else {
    colname = paste0("(",as.character(col_1),",",as.character(col_2),"]")
  }
  
  colour_labels_ac = append(colour_labels_ac,colname)
}


col_breaks_ac = c(rev(brewer.pal(11,"RdYlGn"))) # colour scale green-red (green = produced at high levels in non-invasive populations, low in invasive and visa versa for red)
names(col_breaks_ac) = colour_labels_ac

# point size is based on relative proportion of profiles producing the compound per population (to account for the colour scale removing the 0s)
# also based on quantile breaks
size_breakac = unname(quantile_breaks(cords_spec_1ac$dif_amount,n=11))
size_breakac = c(size_breakac)
size_break_labsac = NULL
for(siz in 2:length(size_breakac)) {
  siz_1 = round(size_breakac[siz-1],digits=3)
  siz_2  = round(size_breakac[siz],digits=3)
  sizname = paste0(as.character(siz_1),"-",as.character(siz_2))
  size_break_labsac = append(size_break_labsac,sizname)
}

# create plots and compute exact test between observed and expected proportions of compounds upregulated and downregualted in invasive populations
#  these are computed per cluster 
exact_test = list() # vector for exact tests of all species

for(spec in species) {
  cords_spec_2ac = cords_spec_1ac[which(cords_spec_1ac$Species %in% spec),]
  
  exp_val = NULL # vector of expected proportions for each cluster 
  dat_val = NULL # vector of compounds upregulated per cluster in the invasive
  dat_val2 = NULL # vector of compounds downregulated per cluster in the invasive
  
  for(clu in unique(chemical_coordinates$cluster)) { # for each cluster of compounds (1, 2 and 3, defined earlier)
    
    comp_spec = unique(comp_feat_csv_cleaned$SMILES[which(comp_feat_csv_cleaned$Species %in% spec)]) # get compounds present per species
    cord_sub = chemical_coordinates[which(chemical_coordinates$name %in% comp_spec),]
    tot_exp = length(which(cord_sub$cluster %in% clu)) # get the number of compounds belonging to the cluster in a species
    prop_exp = tot_exp/nrow(cord_sub) # identifies the proportion of total compound belonging to a cluster (in a given species) - acts as an expected value 
    exp_val = append(exp_val,prop_exp)
    
    # calculate the number of compounds upregulated (difference in production > 1) in the invasive belonging to the cluster
    tot_dat = length(which(cords_spec_2ac$cluster %in% clu & cords_spec_2ac$dif_amount > 1))
    dat_val = append(dat_val,tot_dat)
    
    # calculate the number of compounds dowregulated (difference in production < 1) in the invasive belonging to the cluster
    tot_dat2 = length(which(cords_spec_2ac$cluster %in% clu & cords_spec_2ac$dif_amount < 1))
    dat_val2 = append(dat_val2,tot_dat2)
    
  }

  # calculate the number of compounds upregulated and downregulated per species
  length1 = sum(dat_val)
  length2 = sum(dat_val2)
  print(cbind(dat_val,dat_val2)) # print numbers of compounds upregulated and downregulated per cluster (for tables)
  
  # calculate the expected number of compounds upregulated and downregulated in the invasive, based on how many compounds have altered regulation (assuming it does not differ from the expected proportion of compounds in each cluster in all species)
  expect_data1 = exp_val * length1 
  expect_data2 = exp_val * length2
  
  # multinomial goodness-of-fit tests between expected and observed proportions of clusters
  # add to the exact test list 
  exact1 = xmulti(dat_val,expect_data1,detail=2)
  exact_test[[paste0("up_",spec)]] = exact1
  exact2 = xmulti(dat_val2,expect_data2,detail=2)
  exact_test[[paste0("down_",spec)]] = exact2
  
 # post-hoc tests - binomial tests between observed and expected proportions for each cluster
   for(i in 1:length(dat_val)) {
    
    bin_test1 = binom.test(dat_val[i],length1,exp_val[i],alternative="two.sided") # two-sided has less optimistic p-values)
    exact_test[[paste0("post_hoc_up_",i,"_",spec)]] = bin_test1
    bin_test2 = binom.test(dat_val2[i],length2,exp_val[i],alternative="two.sided")
    exact_test[[paste0("post_hoc_down_",i,"_",spec)]] = bin_test2
  }
  
  # Plotting - plot the coordinates in 2D space from the compound clustering with production levels
  # consistent colour scale
  cords_spec_2ac = transform(cords_spec_2ac, discrete=cut(mean_dif_no0, cords_spec1ac_colscale, include.lowest=T)) # cut data into the quantiles derived for the colourscale
  data_colbreakac = sort(unique(cords_spec_2ac$discrete),decreasing=FALSE) # great the breaks in increasing order
  data_colbreakac = as.character(data_colbreakac)
  cols_usedac = unname(col_breaks_ac[which(names(col_breaks_ac) %in% data_colbreakac)]) # select which colours are used based on the range of quantiles the data falls into 
  
  # get the right sizes for proportion of profiles produced in data, so sizing is consistent between plots
  max_sizeac = round(max(cords_spec_2ac$dif_amount),digits=3)
  
  if(length(which(round(size_breakac,digits=3) == max_sizeac))>0) {
    cut_off = which(round(size_breakac,digits=3) == max_sizeac)
  } else {
    vec_temp = sort(c(max_sizeac,size_breakac),decreasing=FALSE)
    ind = which(vec_temp == max_sizeac)
    cut_off = which(size_breakac == vec_temp[ind+1])
  }
  
  # if a compound has signficantly different regulation in the invasive and non-invasive populations, 
  #   label the corresponding point with the signficance levels using *s
  p_val_labels = NULL
  if(length(which(cords_spec_2ac$p <= 0.05))>0) { # select compounds with a p-value <= 0.05
    for(p in cords_spec_2ac$p[which(cords_spec_2ac$p <= 0.05)]) {
      
      if(p <= 0.0001) {
        p_star = "****"
      } else if (p <= 0.001) {
        p_star = "***"
      } else if(p <= 0.01) {
        p_star = "**"
      } else if(p <= 0.05) {
        p_star = "*"
      }
      p_val_labels = append(p_val_labels,p_star)
    }
  }
  
  if(spec == species[1]) {
    spec_plot = ggplot(cords_spec_2ac,aes(x=x,y=y,color=discrete)) + geom_point(aes(size=dif_amount)) + 
      scale_color_manual(name = "Amplification of Chemical Production in the Invasive",labels=data_colbreakac,
                         values = cols_usedac) + 
      scale_size_continuous(name="Proportion Produced by Invasive/Proportion in non-invasive",
                            breaks=round(size_breakac,digits=3), range = c(1,cut_off)) +
      theme(panel.background =  element_rect(fill = "white",size = 2, linetype = "solid"),axis.line = element_line(colour = "black"),
            legend.title = element_text(size=15), legend.text = element_text(size=11.5), axis.text = element_text(size=13), axis.title = element_text(size=15), 
            plot.title = element_text(size=16)) + ggtitle(bquote(italic(.(spec)))) + ylim(-20,10) +xlim(-15,20) +
      xlab("Dim1 (28.9%)") + ylab("Dim2 (9.8%)")  
    
    
    legend_exp_1 = get_legend(spec_plot) # get the legend for plotting all the plots together 
    
    spec_plot = spec_plot + theme(legend.position = "none")
    
    if(is.null(p_val_labels)==FALSE) {
      spec_plot = spec_plot + geom_text(data=cords_spec_2ac[which(cords_spec_2ac$p <= 0.05),],aes(x=x,y=(y+0.15)),label=p_val_labels,color="black",size=6.5)
    }
    
  } else {
    
    spec_plot = ggplot(cords_spec_2ac,aes(x=x,y=y,color=discrete)) + geom_point(aes(size=dif_amount)) +
      scale_color_manual(name = "Amplification of Chemical Production in the Invasive",labels=data_colbreakac,
                         values = cols_usedac) + 
      scale_size_continuous(name="Proportion Produced by Invasive/Proportion in non-invasive",
                            breaks=round(size_breakac,digits=3), range = c(1,cut_off)) +
      theme(panel.background =  element_rect(fill = "white",size = 2, linetype = "solid"),axis.line = element_line(colour = "black"),legend.position = "none",
            axis.text = element_text(size=13), axis.title = element_text(size=15), plot.title = element_text(size=16)) + ggtitle(bquote(italic(.(spec)))) + ylim(-20,10) +xlim(-15,20) +
      xlab("Dim1 (28.9%)") + ylab("Dim2 (9.8%)") 
    
    if(is.null(p_val_labels)==FALSE) {
      spec_plot = spec_plot + geom_text(data=cords_spec_2ac[which(cords_spec_2ac$p <= 0.05),],aes(x=x,y=(y+0.19)),label=p_val_labels,color="black",size=6.5)
    }
    
  }
  
  assign(paste(spec,"expression_amount_1_ac"),spec_plot)
  
}

exact_test_updown = exact_test # rename exact test to a relevant name

# plot production rate plots together with 1 legend 
grid.arrange(arrangeGrob(`Ageratum conyzoides expression_amount_1_ac`,`Lantana camara expression_amount_1_ac`,
                         `Melaleuca quinquenervia expression_amount_1_ac`,`Psidium cattleianum expression_amount_1_ac`,ncol=2,nrow=2))
grid.arrange(arrangeGrob(legend_exp_1,ncol=1))



################################### Unique Compounds from Clustering ###########################################

# initialise dataframe for compounds unique to populations 
unique_compounds_ac = NULL

for(spec in species) {
  
  for(comp in compound_col) { # for each compound present in the entire dataset
    compound = colnames(matrix_data_feat)[comp] # get the SMILES
    x = cords_spec_ac$x[which(cords_spec_ac$name %in% compound)][[1]] # get the x coordinate from k-medoids clustering
    y = cords_spec_ac$y[which(cords_spec_ac$name %in% compound)][[1]] # get the y coordinate
    clu = cords_spec_ac$cluster[which(cords_spec_ac$name %in% compound)][[1]] # get the cluster the compound belongs to 
    
    # identify which profiles contain the compound
    if(spec == species[4]) { # presence in native and non-native profiles
      inv = c(matrix_data_feat[[compound]][which(matrix_data_feat$Species %in% spec & matrix_data_feat$native_status %in% 1)]) 
      ninv = c(matrix_data_feat[[compound]][which(matrix_data_feat$Species %in% spec & matrix_data_feat$native_status %in% 0)])
      
    } else { # presence in non-invasive and invasive profiles
      inv = c(matrix_data_feat[[compound]][which(matrix_data_feat$Species %in% spec & matrix_data_feat$invasive_status %in% 1)])
      ninv = c(matrix_data_feat[[compound]][which(matrix_data_feat$Species %in% spec & matrix_data_feat$invasive_status %in% 0)])
      
    }
     
    # remove 0 values to check if compounds are truly present
    inv_n0 = inv[which(inv != 0)]
    ninv_n0 = ninv[which(ninv != 0)]
    
    if(length(inv_n0)!=0 | length(ninv_n0)!=0) { #  if the compound is present in either population
      
      if(length(inv_n0)>0 & length(ninv_n0)==0) { # if the compound is found in invasive/non-native populations ONLY
        inv_stat = 1 # add which population it is unique to the dataframe
        no = length(inv_n0)
        tot = length(inv)
        prop = no/tot # proportion of profiles producing the compound
        subdata = c(spec,compound,x,y,clu,inv_stat,no,tot,prop) # add species, coordinates from clustering, cluster number, number of profiles produced in, total profiles in the population and proportion of profiles produced in
        unique_compounds_ac = rbind(unique_compounds_ac,subdata)
        
        
      } else if(length(inv_n0)==0 & length(ninv_n0)>0) { # if the compound is found in non-invasive/native populations ONLY
        inv_stat = 0
        no = length(ninv_n0)
        tot = length(ninv)
        prop = no/tot
        subdata = c(spec,compound,x,y,clu,inv_stat,no,tot,prop)
        unique_compounds_ac = rbind(unique_compounds_ac,subdata)
        
      }
    }
  }
}

unique_compounds_ac = as.data.frame(unique_compounds_ac)
names(unique_compounds_ac) = c("Species","SMILES","x","y","cluster","Invasivity",
                               "No. Profiles Present","Total Profiles of the Phenotype","Proportion Produced In")
# reconvert data back to numeric
unique_compounds_ac$x = as.numeric(unique_compounds_ac$x)
unique_compounds_ac$y = as.numeric(unique_compounds_ac$y)
unique_compounds_ac$`No. Profiles Present` = as.numeric(unique_compounds_ac$`No. Profiles Present`)

# Exact Binomial tests comparing number of compounds present in the invasive compared to non-invasive population per species
exact_test_unique = list() # initialise exact tests lists

for(spec in species) {
  
  chi_unique_ac = as.data.frame(matrix(NA,ncol=2,nrow=3)) # called this because chi-square tests were originally used
  
  for(clu_i in 1:length(unique(unique_compounds_ac$cluster))) { # index for cluster no.
    clu = unique(unique_compounds_ac$cluster)[clu_i] # get cluster 
    count_i = length(which(unique_compounds_ac$cluster %in% clu & unique_compounds_ac$Invasivity %in% 1 & unique_compounds_ac$Species %in% spec)) # calculate no. of unique-to-invasive population clusters
    count_ni = length(which(unique_compounds_ac$cluster %in% clu & unique_compounds_ac$Invasivity %in% 0 & unique_compounds_ac$Species %in% spec)) # calculate no. of unique-to-non-invasive population clusters
    chi_unique_ac[clu_i,] = c(count_i,count_ni)
    
  }
  colnames(chi_unique_ac) = c("Invasive","Non-invasive")
  rownames(chi_unique_ac) = unique(unique_compounds_ac$cluster)
  print(chi_unique_ac)
  
  for(row in 1:nrow(chi_unique_ac)) {
    bin_test = binom.test(chi_unique_ac[row,1],sum(chi_unique_ac[row,]),0.5,alternative="two.sided") # test to see if number of compounds is the same between populations (expected same numbers)
    exact_test_unique[[paste0(rownames(chi_unique_ac)[row],spec)]] = bin_test
    
  }
  
}

# repeat for all species using chi-square rather than exact tests (due to larger sample size)
chi_unique = as.data.frame(matrix(NA,ncol=2,nrow=3))

for(clu_i in 1:length(unique(unique_compounds_ac$cluster))) {
  clu = unique(unique_compounds_ac$cluster)[clu_i]
  count_i = length(which(unique_compounds_ac$cluster %in% clu & unique_compounds_ac$Invasivity %in% 1 ))
  count_ni = length(which(unique_compounds_ac$cluster %in% clu & unique_compounds_ac$Invasivity %in% 0))
  chi_unique[clu_i,] = c(count_i,count_ni)
  
}
colnames(chi_unique) = c("Invasive","Non-invasive")
rownames(chi_unique) = unique(unique_compounds_ac$cluster)

chisq.test(chi_unique)
chisq.posthoc.test(chi_unique) # post-hoc analysis per cluster






