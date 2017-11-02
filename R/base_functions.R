# Includes functions for reading and plotting



#
# General variables
#
.fill_vals = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(12, "Set3")[-9],
              RColorBrewer::brewer.pal(8, "Pastel1"), RColorBrewer::brewer.pal(7, "Pastel2"), RColorBrewer::brewer.pal(7, "Dark2"),
              RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(7, "Accent"),
              RColorBrewer::brewer.pal(11, "Spectral"))
.shade_cols = c("grey9", "grey29", "grey55", "grey78", "grey91")


#
# Process recombinants file
#
## Adapted from: https://stackoverflow.com/questions/12626637/reading-a-text-file-in-r-line-by-line
processRecomb = function(filepath){
  con = file(filepath, "r")
  cells_dic = list()
  while(TRUE){
    line = readLines(con, n = 1)
    if(length(line) == 0){ #EOF
      break
    } else if(grepl("\t", line) & !grepl("cell_name", line)){
      line = unlist(strsplit(gsub("\n", "", line), "\t"))
      cell = line[1]
      if(!cell %in% names(cells_dic)){
        cells_dic[[cell]] = c("notAssigned", "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No")
      }
      if(startsWith(line[2], "No seqs")){
        cells_dic[[cell]] = c("NoTCR", "NA", "NA", "NA", "NA",
                           "NA", "NA", "NA", "NA",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No")
      } else if(line[2]=="A"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")

        if(cells_dic[[cell]][2]=="No"){
          cells_dic[[cell]][2] = seq
          cells_dic[[cell]][10] = prod
        } else{
          cells_dic[[cell]][3] = seq
          cells_dic[[cell]][11] = prod
        }
      } else if(line[2]=="B"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")

        if(cells_dic[[cell]][4]=="No"){
          cells_dic[[cell]][4] = seq
          cells_dic[[cell]][12] = prod
        } else{
          cells_dic[[cell]][5] = seq
          cells_dic[[cell]][13] = prod
        }
      } else if(line[2]=="G"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")

        if(cells_dic[[cell]][6]=="No"){
          cells_dic[[cell]][6] = seq
          cells_dic[[cell]][14] = prod
        } else{
          cells_dic[[cell]][7] = seq
          cells_dic[[cell]][15] = prod
        }
      } else if(line[2]=="D"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")

        if(cells_dic[[cell]][8]=="No"){
          cells_dic[[cell]][8] = seq
          cells_dic[[cell]][16] = prod
        } else{
          cells_dic[[cell]][9] = seq
          cells_dic[[cell]][17] = prod
        }
      }
    }
  }
  close(con)
  return(cells_dic)
}


#
# Process TCR summary file
#
## Adapted from: https://stackoverflow.com/questions/12626637/reading-a-text-file-in-r-line-by-line
processTCR = function(filepath, cells_dic){
  con = file(filepath, "r")
  bool_vec = c("multi" = F, "iNKT" = F, "MAIT" = F, "clon" = F)
  cl_n = 0
  while(TRUE){
    line = readLines(con, n = 1)
    if(length(line) == 0){
      break
    }

    if(grepl("Cells with", line)){
      bool_vec["multi"] = T
    } else if(grepl("iNKT", line)){
      bool_vec["multi"] = F
      bool_vec["iNKT"] = T
    } else if(grepl("MAIT", line)){
      bool_vec["iNKT"] = F
      bool_vec["MAIT"] = T
    } else if(grepl("Clonotype", line)){
      bool_vec["iNKT"] = F
      bool_vec["MAIT"] = F
      bool_vec["clon"] = T
    }

    if(startsWith(line, "###")){
      cell = gsub(" ###", "", gsub("### ", "", line))
      if(!cell %in% names(cells_dic)){
        cells_dic[[cell]] = c("notAssigned", "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No")
      }
      if(bool_vec["multi"]){
        cells_dic[[cell]][1] = "Multi_recomb"
      }else if(bool_vec["iNKT"]){
        cells_dic[[cell]][1] = "iNKT"
      }else if(bool_vec["MAIT"]){
        cells_dic[[cell]][1] = "MAIT"
      }
    }
    if(bool_vec["clon"] & grepl(",", line)){
      cl_n = cl_n + 1
      line = unlist(strsplit(gsub("\n", "", gsub(" ", "", line)), ","))
      for(cell in line){
        cells_dic[[cell]][1] = paste0("cl", cl_n)
      }
    }
  }

  close(con)
  return(cells_dic)
}


#
# Obtain matrix of matching chains between pairs of cells reported as sharing chains
#
matchingChainsMatrix <- function(tracerData){
  cl_tracer_data = tracerData[grepl("cl", tracerData$tcr_info),]

  # Matrix of chains shared
  labs_chains = rep(c("A", "B", "G", "D"), each = 2)
  mat_chains = matrix("", nrow(cl_tracer_data), nrow(cl_tracer_data))
  rownames(mat_chains) = rownames(cl_tracer_data)
  colnames(mat_chains) = rownames(cl_tracer_data)
  cell_comb = combn(1:nrow(cl_tracer_data),2)
  for(comb_col in 1:ncol(cell_comb)){
    i = cell_comb[1,comb_col]
    j = cell_comb[2,comb_col]
    for(k in 2:9){ # Go through detected chains for each cell
      if(as.character(cl_tracer_data[i,k]) %in% as.character(unlist(cl_tracer_data[j,2:9])) & cl_tracer_data[i,k]!="No"){
        newchain = labs_chains[k-1]
        posj = which(as.character(unlist(cl_tracer_data[j,2:9]))==as.character(cl_tracer_data[i,k]))+9
        mat_chains[i,j] = paste0(mat_chains[i,j],
                                 ifelse(cl_tracer_data[i,k+8]=="No" | cl_tracer_data[j,posj]=="No",
                                        tolower(newchain), newchain))
        mat_chains[j,i] = paste0(mat_chains[j,i],
                                 ifelse(cl_tracer_data[j,posj]=="No" | cl_tracer_data[i,k+8]=="No",
                                        tolower(newchain), newchain))
      }
    }
  }
  return(mat_chains)
}


#
# Define Clonotypes based on specific criteria
#
defineClonotypes <- function(tracerData, matChains, criteriaList, nameVar = "custom_cl"){
  chains_tcr = c("A", "B", "G", "D")

  # Reformat matrix
  mat_melt = reshape2::melt(matChains)
  mat_melt = mat_melt[mat_melt[,1]!=mat_melt[,2] & mat_melt[,3]!="",]

  # Apply filters
  cond_list = list()
  for(c_n in 1:length(criteriaList)){ # each filter
    cond_list[[c_n]] = T
    crit = criteriaList[[c_n]]
    for(chain in 1:length(crit[[1]])){ # each chain
      if(crit[[1]][chain]){ # chain is required
        if(crit[[2]][chain]){ # chain has to be productive
          cond_list[[c_n]] = cond_list[[c_n]] & grepl(chains_tcr[chain], mat_melt[,3])
        }else{
          cond_list[[c_n]] = cond_list[[c_n]] & (grepl(chains_tcr[chain], mat_melt[,3]) |
                                                   grepl(tolower(chains_tcr[chain]), mat_melt[,3]))
        }
      }
    }
  }
  condition = F
  for(i in cond_list){
    condition = condition | i
  }
  mat_melt = mat_melt[condition,]

  # regroup new clonotypes
  clusters <- igraph::clusters(igraph::graph.data.frame(mat_melt[,1:2]))$membership
  clusters = data.frame(row.names = names(clusters), nameVar = paste0("cl",clusters), stringsAsFactors = F)
  colnames(clusters) = nameVar

  tracerData = merge(tracerData, clusters, by = 0, all = T)
  rownames(tracerData) = tracerData[,1]
  tracerData = tracerData[,-1]

  # labels for cells that are not in a clonotype
  ## cells that were in a clonotype but are not anymore will be assigned as "notAssigned"
  tracerData[is.na(tracerData[,nameVar]) & grepl("cl", tracerData[,"tcr_info"]), nameVar] = "notAssigned"
  tracerData[is.na(tracerData[,nameVar]), nameVar] = as.character(tracerData[is.na(tracerData[,nameVar]),1])

  return(tracerData)
}


#
# Read TraCeR results
#
readTracer <- function(summaryPath,
                       recombinants = "recombinants.txt", tcrSum = "TCR_summary.txt",
                       clonotypes = NULL, nameVar = "custom_cl",
                       get_matching_chains = F) {

  # Process recombinants files (in function argument) and TCR info (iNKT, clonotypes, ...)
  processed_files = processTCR(paste0(summaryPath, tcrSum),
                               processRecomb(paste0(summaryPath, recombinants)))

  # Reformat as a data frame
  tracer_data = t(data.frame(processed_files))
  colnames(tracer_data) = paste0("tcr_", c("info", "A_1", "A_2", "B_1", "B_2",
                            "G_1", "G_2", "D_1", "D_2",
                            "pA_1", "pA_2", "pB_1", "pB_2",
                            "pG_1", "pG_2", "pD_1", "pD_2"))
  tracer_data = data.frame(tracer_data)

  # Add clonotype information according to the criteria defined
  # Examples
  # criteriaList = list(list(c(T, T, F, F), c(T, T, F, F))) - only A and B productive
  # criteriaList = list(list(c(T, F, F, F), c(T, F, F, F)), list(c(F, T, F, F), c(F, T, F, F))) - either A or B productive
  if(get_matching_chains | !is.null(clonotypes)){
    chains_mat = matchingChainsMatrix(tracer_data)
    if(!is.null(clonotypes)){
      tracer_data = defineClonotypes(tracer_data, chains_mat,
                                     clonotypes, nameVar = nameVar)
    }
  }


  # Define final object to be returned
  result = list("tracer_metadata" = tracer_data)
  ## add matching chains matrix
  result$matching_chains = if(get_matching_chains | !is.null(clonotypes)) chains_mat else NULL
  ## add VDJ segments table slot
  result$vdj_segments = NULL

  return(result)
}


#
# Project clonotypes in dim red plot
#
plotProjection <- function(tracer_data, pheno_data,
                           dimensions = c("dim1", "dim2"), clonotypes = "tcr_info",
                           additional_pheno = NULL, plot_out = T){
  if(!is.null(additional_pheno) & length(additional_pheno)>1){
    stop("Please supply only one variable as additional information.")
  }
  if(length(unique(pheno_data[,additional_pheno]))>10){
    stop("No more than 11 classes can be shown as shapes.")
  }

  points_df = data.frame(row.names = rownames(tracer_data$tracer_metadata),
                         "tcr_info" = tracer_data$tracer_metadata[,clonotypes])
  points_df = merge(pheno_data, points_df, by = 0, all.x = T)
  points_df = points_df[,c("tcr_info", dimensions, additional_pheno)]
  points_df$tcr_info = as.character(points_df$tcr_info)
  # Make sure additional_pheno is a factor
  if(!is.null(additional_pheno) & !is.factor(points_df[,additional_pheno])){
    points_df[,additional_pheno] = as.factor(points_df[,additional_pheno])
  }
  # label cells that do not have a TraCeR output
  if(sum(is.na(points_df$tcr_info))>1){
    warning("Not all cells have a TraCeR output.")
    points_df$tcr_info[is.na(points_df$tcr_info)] = "NO TraCeR"
  }

  lvl_cl = unique(points_df$tcr_info)
  non_cl = c("notAssigned", "iNKT", "MAIT", "NoTCR", "NO TraCeR")
  points_df$tcr_info = factor(points_df$tcr_info, levels = c(lvl_cl[!lvl_cl %in% non_cl], non_cl))
  points_df$iscl = grepl("cl", points_df$tcr_info)
  points_df$iscl = factor(points_df$iscl, levels = c(F, T))

  # DF to draw lines - only for clonotypes
  lines_df = points_df[grepl("cl", points_df$tcr_info),]

  # define colours and shapes for plot
  if(length(fill_vals)<length(unique(lines_df$tcr_info))){
    fill_vals_plt = c(fill_vals, viridisLite::viridis(length(unique(lines_df$tcr_info))-length(fill_vals)),
                  shade_cols)
  } else{
    fill_vals_plt = c(fill_vals[seq(1, length(unique(lines_df$tcr_info)))],
                      shade_cols)
  }
  shape_vals = c(15, 19:17, 0:2, 5, 6, 8)

  colour_l_vals = fill_vals_plt

  # Make plot
  plot_proj = ggplot() +
    geom_line(data = lines_df, mapping = aes_string(x = dimensions[1], y = dimensions[2],
                                                    colour = "tcr_info"), size = 1)
  ## are additional features defined?
  plot_proj = if(!is.null(additional_pheno)){
    plot_proj + geom_point(data = points_df, mapping = aes_string(x = dimensions[1], y = dimensions[2],
                                                                  colour = "tcr_info", shape = additional_pheno,
                                                                  size = "iscl", alpha = "iscl"))
  } else{
    plot_proj + geom_point(data = points_df, mapping = aes_string(x = dimensions[1], y = dimensions[2],
                                                                  colour = "tcr_info", size = "iscl", alpha = "iscl"))
  }
  ## cont
  plot_proj = plot_proj+
    # scales
    scale_colour_manual(values = colour_l_vals, drop=FALSE)+
    scale_shape_manual(values = shape_vals, drop=FALSE)+
    scale_size_manual(values = c(1.2, 2), guide = "none")+
    scale_alpha_manual(values = c(0.75, 1), guide = "none")+
    guides(shape = guide_legend(order = 1),
           colour = guide_legend(order = 2, title = clonotypes))+
    theme_classic()+
    theme(legend.text = element_text(size = 8, colour = "black"),
          legend.title = element_text(size = 9, colour = "black"),
          axis.title = element_text(size = 9, colour = "black"),
          axis.text = element_text(size = 7.7, colour = "black"),
          legend.position="right",
          legend.key.size = unit(0.3, "cm"),
          legend.box.spacing = unit(0.03, "cm"))

  if(plot_out){
    print(plot_proj)
  }
  return(plot_proj)
}

#
# Count clonotypes by category
#
plotCounts <- function(tracer_data, pheno_data, category,
                       clonotypes = "tcr_info", plot_out = T){
  if(length(category)>2){
    stop("Choose at most 2 categories for plotting.")
  }

  clon_df = data.frame(row.names = rownames(tracer_data$tracer_metadata),
                       "tcr_info" = tracer_data$tracer_metadata[,clonotypes])
  clon_df = merge(pheno_data, clon_df, by = 0, all.x = T)
  clon_df = clon_df[,c("tcr_info", category)]
  #clon_df = clon_df[grepl("cl", clon_df$tcr_info),]
  colnames(clon_df) = if(length(category)==2) c("tcr_info", "cat1", "cat2") else c("tcr_info", "cat1")

  lvl_cl = levels(clon_df$tcr_info)
  non_cl = c("notAssigned", "iNKT", "MAIT", "NoTCR", "NO TraCeR")
  clon_df$tcr_info = factor(clon_df$tcr_info, levels = c(lvl_cl[!lvl_cl %in% non_cl], non_cl))

  # define colours and shapes for plot
  if(length(fill_vals)<length(unique(clon_df$tcr_info))){
    fill_vals_plt = c(fill_vals, viridisLite::viridis(length(unique(clon_df$tcr_info))-length(fill_vals)),
                  shade_cols)
  } else{
    fill_vals_plt = c(fill_vals[seq(1, length(unique(clon_df$tcr_info[grepl("cl", clon_df$tcr_info)])))],
                  shade_cols)
  }

  # Make plot
  plot_count = ggplot(clon_df, aes(x = cat1, fill = tcr_info))
  plot_count = if(length(category)==2){
    plot_count + facet_wrap(~ cat2)
  } else {plot_count}

  plot_count = plot_count + geom_bar()+
    guides(fill = guide_legend(title = clonotypes))+
    # scales
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(name = category[1])+
    scale_fill_manual(values = fill_vals_plt)+
    theme_classic()+
    theme(legend.text = element_text(size = 8, colour = "black"),
          legend.title = element_text(size = 9, colour = "black"),
          axis.title = element_text(size = 9, colour = "black"),
          axis.text = element_text(size = 7.7, colour = "black"),
          legend.position="right",
          legend.key.size = unit(0.3, "cm"),
          legend.box.spacing = unit(0.03, "cm"))

  if(plot_out){
    print(plot_count)
  }
  return(plot_count)
}



#
# Matrix plot showing shared chains
#


#
# Sankey diagramme of clonotype sharing
#


#
# Sankey diagramme of VDJ usage
#


