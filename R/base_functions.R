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
# Define cell chain code
#
.cellCode <- function(tracer_data){
  chains_ref = rep(c("A","B","G","D"), each = 2)

  chains_col = c()
  for(i in 1:nrow(tracer_data)){
    if(tracer_data[i, "tcr_info"]=="notAssigned" |
       grepl("cl", tracer_data[i, "tcr_info"])){
      chains = tracer_data[i,2:9]
      chains_code = chains_ref[which(chains!="No")]
      prod = tracer_data[i,10:17][which(chains!="No")]
      chains_code[prod=="No"] = tolower(chains_code[prod=="No"])
      chains_col = c(chains_col, paste0(chains_code, collapse=""))
    } else{
      chains_col = c(chains_col, as.character(tracer_data[i, "tcr_info"]))
    }
  }

  return(chains_col)
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

  message("Obtaining matching chains for clonotypes...")
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
# Cell labeler - used to label cells as gamma-delta
#
cellLabeler <- function(tracerData, ref_col = "tcr_info", new_col = "tcr_info_mod",
                        presentChains = list(c("G", "g"), c("D", "d")), # each group 'or', between 'and'
                        absentChains = list(c("A", "a"), c("B", "b")),
                        new_label = "gdCell"){
  # check if new column exists - issue warning
  if(ref_col==new_col | new_col %in% colnames(tracerData$tracer_metadata)){
    warning("You are overwriting an already existing column in tracer_metadata.")
  }

  # both conditions can't be NULL
  if(is.null(presentChains) & is.null(absentChains)){
    stop("Both conditions can not be NULL.")
  }

  # chain is present
  if(!is.null(presentChains)){
    chain_codes = tracerData$tracer_metadata$tcrCode
    present_cond = list()
    for(gg in 1:length(presentChains)){
      start_cond = rep(F, length(chain_codes))
      for(cc in presentChains[[gg]]){
        start_cond = start_cond | grepl(cc, chain_codes)
        chain_codes = sub(cc, "", chain_codes)
      }
      present_cond[[gg]] = start_cond
    }
    p_c = T
    for(pc in present_cond){ p_c = p_c & pc }
  }

  # chain is absent
  if(!is.null(absentChains)){
    chain_codes = tracerData$tracer_metadata$tcrCode
    absent_cond = list()
    for(gg in 1:length(absentChains)){
      start_cond = rep(F, length(chain_codes))
      for(cc in absentChains[[gg]]){
        start_cond = start_cond | !grepl(cc, chain_codes)
        chain_codes = sub(cc, "", chain_codes)
      }
      absent_cond[[gg]] = start_cond
    }
    a_c = T
    for(ac in absent_cond){ a_c = a_c & ac }
  }

  # merge info and classify
  final_cond = if(!is.null(presentChains) & !is.null(absentChains)){
    p_c & a_c
  } else if(!is.null(presentChains)){
    p_c
  } else if(!is.null(absentChains)){
    a_c
  }
  newcol = as.character(tracerData$tracer_metadata[,ref_col])
  newcol[final_cond] = new_label

  # add column and return tracerData
  tracerData$tracer_metadata[,new_col] = newcol
  return(tracerData)

}


# CHANGE
# Obtain matrix of matching chains between pairs of cells reported as sharing chains
#
segmentMatrix <- function(tracerData){
  chain_cols = c("tcr_A_1", "tcr_A_2", "tcr_B_1", "tcr_B_2",
                 "tcr_G_1", "tcr_G_2", "tcr_D_1", "tcr_D_2")
  prod_cols = c("tcr_pA_1", "tcr_pA_2", "tcr_pB_1", "tcr_pB_2",
                "tcr_pG_1", "tcr_pG_2", "tcr_pD_1", "tcr_pD_2")
  seg_df = tracerData$tracer_metadata[,-1]
  seg_df$cell_name = rownames(seg_df)

  # transform matrix
  seg_df_c = reshape2::melt(seg_df[,c(1:8, 17)], id.vars = "cell_name")
  seg_df_p = reshape2::melt(seg_df[,9:17], id.vars = "cell_name")

  l_segs = strsplit(seg_df_c$value, "_")
  l_segs = lapply(l_segs, function(x) if(length(x)<3) rep(x, 3) else x) #CLARIFY HOW THE SEGMENTS WORK

  seg_df_c$productive = seg_df_p$value
}


#
# Define Clonotypes based on specific criteria
#
defineClonotypes <- function(tracerData, matChains, criteriaList,
                             nameVar = "custom_cl"){
  chains_tcr = c("A", "B", "G", "D")

  # Reformat matrix
  mat_melt = reshape2::melt(matChains)
  mat_melt = mat_melt[mat_melt[,1]!=mat_melt[,2] & mat_melt[,3]!="",]

  # Apply filters
  chain_codes = mat_melt$value
  present_cond = list()
  for(gg in 1:length(criteriaList)){
    start_cond = rep(F, length(chain_codes))
    for(cc in criteriaList[[gg]]){
      start_cond = start_cond | grepl(cc, chain_codes)
      chain_codes = sub(cc, "", chain_codes)
    }
    present_cond[[gg]] = start_cond
  }
  p_c = T
  for(pc in present_cond){ p_c = p_c & pc }
  mat_melt = mat_melt[p_c,]

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
# Correct clonotypes
#
.correctClonotypes <- function(tcr_info_column){
  if(any(table(tcr_info_column)==1)){
    res = tcr_info_column
    res[res %in% names(table(res)[table(res)==1])] = "notAssigned"
    message('Some cells with detected clonotypes had no match in the current selection and are masked as "notAssigned"')
    return(res)
  }
}


#
# Read TraCeR results
#
readTracer <- function(summaryPath,
                       recombinants = "recombinants.txt", tcrSum = "TCR_summary.txt",
                       clonotypes = NULL, nameVar = "custom_cl",
                       get_matching_chains = F,
                       get_individual_segments = F) {

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
  tracer_data$tcrCode = .cellCode(tracer_data)

  # Add clonotype information according to the criteria defined
  # Examples
  # criteriaList = list(c("A"), c("B")) - only A and B productive
  # criteriaList = list(c("A", "B")) - either A or B productive
  if(get_matching_chains | !is.null(clonotypes)){
    chains_mat = matchingChainsMatrix(tracer_data)
    if(!is.null(clonotypes)){
      tracer_data = defineClonotypes(tracer_data, chains_mat,
                                     clonotypes, nameVar = nameVar)
    }
  }

  # Add individual segment information
  if(get_individual_segments){

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
  # annotate single clonotypes as "notAssigned"
  points_df$tcr_info = .correctClonotypes(points_df$tcr_info)

  lvl_cl = unique(points_df$tcr_info)
  non_cl = c("notAssigned", "iNKT", "MAIT", "NoTCR", "NO TraCeR")
  points_df$tcr_info = factor(points_df$tcr_info, levels = c(lvl_cl[!lvl_cl %in% non_cl], non_cl))
  points_df$iscl = grepl("cl", points_df$tcr_info)
  points_df$iscl = factor(points_df$iscl, levels = c(F, T))

  # DF to draw lines - only for clonotypes
  lines_df = points_df[grepl("cl", points_df$tcr_info),]

  # define colours and shapes for plot
  if(length(.fill_vals)<length(unique(lines_df$tcr_info))){
    fill_vals_plt = c(.fill_vals, viridisLite::viridis(length(unique(lines_df$tcr_info))-length(.fill_vals)),
                  .shade_cols)
  } else{
    fill_vals_plt = c(.fill_vals[seq(1, length(unique(lines_df$tcr_info)))],
                      .shade_cols)
  }
  shape_vals = c(15, 19:17, 0:2, 5, 6, 8)

  colour_l_vals = fill_vals_plt

  # Make plot
  plot_proj = ggplot2::ggplot() +
    ggplot2::geom_line(data = lines_df, mapping = ggplot2::aes_string(x = dimensions[1], y = dimensions[2],
                                                    colour = "tcr_info"), size = 1)
  ## are additional features defined?
  plot_proj = if(!is.null(additional_pheno)){
    plot_proj + ggplot2::geom_point(data = points_df, mapping = ggplot2::aes_string(x = dimensions[1], y = dimensions[2],
                                                                  colour = "tcr_info", shape = additional_pheno,
                                                                  size = "iscl", alpha = "iscl"))
  } else{
    plot_proj + ggplot2::geom_point(data = points_df, mapping = ggplot2::aes_string(x = dimensions[1], y = dimensions[2],
                                                                  colour = "tcr_info", size = "iscl", alpha = "iscl"))
  }
  ## cont
  plot_proj = plot_proj+
    # scales
    ggplot2::scale_colour_manual(values = colour_l_vals, drop=FALSE)+
    ggplot2::scale_shape_manual(values = shape_vals, drop=FALSE)+
    ggplot2::scale_size_manual(values = c(1.2, 2), guide = "none")+
    ggplot2::scale_alpha_manual(values = c(0.75, 1), guide = "none")+
    ggplot2::guides(shape = ggplot2::guide_legend(order = 1),
           colour = ggplot2::guide_legend(order = 2, title = clonotypes))+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8, colour = "black"),
          legend.title = ggplot2::element_text(size = 9, colour = "black"),
          axis.title = ggplot2::element_text(size = 9, colour = "black"),
          axis.text = ggplot2::element_text(size = 7.7, colour = "black"),
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
  colnames(clon_df) = if(length(category)==2) c("tcr_info", "cat1", "cat2") else c("tcr_info", "cat1")
  # annotate single clonotypes as "notAssigned"
  clon_df$tcr_info = .correctClonotypes(clon_df$tcr_info)

  lvl_cl = levels(clon_df$tcr_info)
  non_cl = c("notAssigned", "iNKT", "MAIT", "NoTCR", "NO TraCeR")
  clon_df$tcr_info = factor(clon_df$tcr_info, levels = c(lvl_cl[!lvl_cl %in% non_cl], non_cl))

  # define colours and shapes for plot
  if(length(.fill_vals)<length(unique(clon_df$tcr_info))){
    fill_vals_plt = c(.fill_vals, viridisLite::viridis(length(unique(clon_df$tcr_info))-length(.fill_vals)),
                  .shade_cols)
  } else{
    fill_vals_plt = c(.fill_vals[seq(1, length(unique(clon_df$tcr_info[grepl("cl", clon_df$tcr_info)])))],
                  .shade_cols)
  }

  # Make plot
  plot_count = ggplot2::ggplot(clon_df, ggplot2::aes(x = cat1, fill = tcr_info))
  plot_count = if(length(category)==2){
    plot_count + ggplot2::facet_wrap(~ cat2)
  } else {plot_count}

  plot_count = plot_count + ggplot2::geom_bar()+
    ggplot2::guides(fill = ggplot2::guide_legend(title = clonotypes))+
    # scales
    ggplot2::scale_y_continuous(expand = c(0,0))+
    ggplot2::scale_x_discrete(name = category[1])+
    ggplot2::scale_fill_manual(values = fill_vals_plt)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8, colour = "black"),
          legend.title = ggplot2::element_text(size = 9, colour = "black"),
          axis.title = ggplot2::element_text(size = 9, colour = "black"),
          axis.text = ggplot2::element_text(size = 7.7, colour = "black"),
          legend.position="right",
          legend.key.size = ggplot2::unit(0.3, "cm"),
          legend.box.spacing = ggplot2::unit(0.03, "cm"))

  if(plot_out){
    print(plot_count)
  }
  return(plot_count)
}


#
# Matrix plot showing shared chains
#
plotSharedClones <- function(tracer_data, pheno_data, category,
                             clonotypes = "tcr_info", plot_out = T){
  shared_list = list()

  clon_df = data.frame(row.names = rownames(tracer_data$tracer_metadata),
                       "cl1" = tracer_data$tracer_metadata[,clonotypes[1]])
  colnames(clon_df)[1] = clonotypes[1]
  clon_df = merge(pheno_data, clon_df, by = 0, all.x = T)
  clon_df = clon_df[,c(clonotypes, category)]
  # annotate single clonotypes as "notAssigned"
  clon_df$tcr_info = .correctClonotypes(clon_df$tcr_info)

  sub_clon_df = clon_df
  colnames(sub_clon_df)[which(grepl(clonotypes, x = colnames(sub_clon_df)))] = "classes"
  colnames(sub_clon_df)[which(grepl(category, x = colnames(sub_clon_df)))] = "cond"
  sub_clon_df$classes = factor(sub_clon_df$classes,
                              levels = unique(sub_clon_df$classes)[order(unique(as.character(sub_clon_df$classes)))])
  sub_clon_df$cond = factor(sub_clon_df$cond)
  sub_clon_df = sub_clon_df[grepl("cl", sub_clon_df$classes),]

  matrix_counts = matrix(rep(0, length(unique(sub_clon_df$cond))**2),
                         length(unique(sub_clon_df$cond)), length(unique(sub_clon_df$cond)))
  rownames(matrix_counts) = unique(sub_clon_df$cond)
  colnames(matrix_counts) = unique(sub_clon_df$cond)
  for(cl in unique(sub_clon_df$classes)){
    sub_pheno = sub_clon_df[sub_clon_df$classes==cl, "cond"]

    comb_pheno = t(combn(sub_pheno, 2))
    comb_pheno.sort = t(apply(comb_pheno, 1, sort))
    cor_mat = comb_pheno[!duplicated(comb_pheno.sort),]
    if(!is.matrix(cor_mat)){
      cor_mat = matrix(cor_mat, 1, 2)
    }
    for(i in 1:nrow(cor_mat)){
      matrix_counts[cor_mat[i,1],cor_mat[i,2]] = matrix_counts[cor_mat[i,1],cor_mat[i,2]] +1
      if(cor_mat[i,2]!=cor_mat[i,1]){
        matrix_counts[cor_mat[i,2],cor_mat[i,1]] = matrix_counts[cor_mat[i,2],cor_mat[i,1]] +1
      }
    }
  }
  matrix_list = matrix_counts

  inter_pop_cl = matrix_list
  lab_mat = matrix(as.character(inter_pop_cl), nrow(inter_pop_cl), ncol(inter_pop_cl))
  rownames(lab_mat) = rownames(inter_pop_cl)
  colnames(lab_mat) = colnames(inter_pop_cl)

  inter_pop_cl = reshape2::melt(inter_pop_cl)
  lab_mat = reshape2::melt(lab_mat)

  plot_shared = ggplot2::ggplot() +
    ggplot2::geom_point(data = inter_pop_cl,
               mapping = ggplot2::aes(x = as.character(Var1), y = as.character(Var2), size = as.integer(value)))+
    ggplot2::geom_text(data = lab_mat,
              mapping = ggplot2::aes(x = as.character(Var1), y = as.character(Var2), label = value),
              colour = "white", size = 4.3, fontface = "bold")+
    ggplot2::scale_x_discrete(name = category)+
    ggplot2::scale_y_discrete(name = category)+
    ggplot2::scale_size_continuous(range = c(4.9, 11))+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 10, colour = "black"),
          axis.text = ggplot2::element_text(size = 9, colour = "black"),
          axis.line = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "grey72"),
          panel.grid.major = ggplot2::element_line(colour = "grey90", linetype = "dashed"),
          legend.position="none")

  if(plot_out){
    print(plot_shared)
  }
  return(plot_shared)

}


#
# Table for matching chains REMOVE CELLS WITHOUT MATCH AFTER FILTERING
#
plotMatchingChains <- function(tracer_data, pheno_data = NULL,
                               categories = NULL, cell_names = T,
                               plot_out = T){
  if(!is.null(pheno_data)){
    rnpheno = rownames(pheno_data)
    tcr_info = cbind(rnpheno, as.character(tracer_data$tracer_metadata[rnpheno,"tcr_info"])) # use this to filter
    tcr_info = tcr_info[grepl("cl", tcr_info[,2]),] # only clonotypes (all)
    ti_tab = table(tcr_info[,2]) # get cells that after filtering have clonotype
    rnpheno = tcr_info[which(tcr_info[,2] %in% names(ti_tab[ti_tab>1])),1] # only cells in clonotypes (after filtering)
    tracer_data$matching_chains = tracer_data$matching_chains[rownames(tracer_data$matching_chains) %in% rnpheno,
                                                              colnames(tracer_data$matching_chains) %in% rnpheno]
  }
  mc_count = nchar(tracer_data$matching_chains)
  mc_text = tracer_data$matching_chains
  # cluster shared chains
  hc = hclust(dist(mc_count,
                   method = "manhattan"),
              method = "complete")
  mc_count = reshape2::melt(mc_count[hc$order,hc$order])
  mc_count[mc_count==0] = NA
  mc_text =  reshape2::melt(tracer_data$matching_chains[hc$order,hc$order])

  # plot
  mc_plot = ggplot2::ggplot()+
    ggplot2::geom_tile(data = mc_count, mapping = ggplot2::aes(x = Var1, y = Var2, fill = value))+
    ggplot2::geom_text(data = mc_text,
              mapping = ggplot2::aes(x = Var1, y = Var2, label = value),
              colour = "white", size = 1.5)+
    ggplot2::scale_fill_continuous(low = "#078ff2", high = "#241ece", na.value = "white")+
    ggplot2::scale_y_discrete(name = "Cells")+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 9, colour = "black"),
          axis.text.y = ggplot2::element_text(size = 7.5, colour = "black"),
          axis.text.x = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(l = 0.07, r = 0.07),
          panel.background = ggplot2::element_rect(fill = "grey72"),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA),
          legend.position="none")

  if(cell_names==FALSE){
    mc_plot = mc_plot + ggplot2::theme(axis.text = ggplot2::element_blank(),
                              axis.title = ggplot2::element_blank(),
                              axis.ticks = ggplot2::element_blank())
  }

  if(!is.null(pheno_data) & !is.null(categories)){
    info_data = data.frame(pheno_data[as.character(unique(mc_count$Var1)),categories])
    info_data$Cells = unique(mc_count$Var1)
    colnames(info_data) = c(categories, "Cells")
    info_data = reshape2::melt(info_data, id.var = "Cells")

    sideplot = ggplot2::ggplot()+
      ggplot2::geom_tile(data = info_data,
                mapping = ggplot2::aes(x = variable, y = Cells, fill = value))+
      ggplot2::scale_x_discrete(expand = c(0,0))+
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Values"))+
      ggplot2::theme_classic()+
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9),
            legend.text = element_text(size = 7.1),
            legend.key.size = ggplot2::unit(0.31, "cm"),
            legend.box.margin = ggplot2::margin(l = 0.05, r = 0.05),
            legend.margin = ggplot2::margin(l = 0.05, r = 0.05),
            legend.box.spacing = ggplot2::unit(0.19, "cm"),
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(colour = "black", size = 7.8,
                                       angle = 30, vjust = 0.9, hjust = 0.9))

    result = cowplot::plot_grid(mc_plot, sideplot, ncol = 2,
                                rel_widths = c(1, 0.21), align = "h")
  if(plot_out) print(result)
    return(result)

  } else{
    if(plot_out) print(mc_plot)
    return(mc_plot)
  }

}


#
# Clonotype sizes
#
plotClonotypeSize <- function(tracer_data, pheno_data = NULL, clonotypes = "tcr_info",
                          category = NULL, plot_out = T){
  if(!is.null(pheno_data) & !is.null(category)){
    plot_df = merge(tracer_data$tracer_metadata, pheno_data, by = 0)[,c(category, clonotypes)]
    plot_df[,clonotypes] = .correctClonotypes(plot_df[,clonotypes])

    plot_df_cl = plot_df[grepl("cl", as.character(plot_df[,clonotypes])),]
    plot_df_cl = reshape2::melt(lapply(tapply(plot_df_cl[,clonotypes], plot_df_cl[,category],
                                              function(x) table(table(x))), cbind))[,-2]
    plot_df_cl = plot_df_cl[plot_df_cl[,1]!=0,]

    plot_df_na = plot_df[grepl("notAssigned", as.character(plot_df[,clonotypes])),]
    plot_df_na = reshape2::melt(lapply(tapply(plot_df_na[,clonotypes], plot_df_na[,category],
                                              function(x) table(table(x))), cbind))[,-2]
    plot_df_na = plot_df_na[plot_df_na[,1]!=0,]
    plot_df_na[,2] = plot_df_na[,1]
    plot_df_na[,1] = 1

    plot_df = rbind(plot_df_na, plot_df_cl)
    colnames(plot_df) = c("size", "Number", "cat")
  } else{
    tab_all = table(tracer_data$tracer_metadata[,clonotypes])
    na_n = tab_all["notAssigned"]
    tab_all = tab_all[grepl("cl", names(tab_all))]

    plot_df = data.frame("size" = as.numeric(c(1, names(table(tab_all)))),
                         "Number" = c(na_n, table(tab_all)))
  }

  plot_cls = ggplot2::ggplot(plot_df, ggplot2::aes(x = size, y = Number))+
    ggplot2::geom_bar(stat = "identity", fill = "grey2")+
    ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, max(plot_df$Number)+max(plot_df$Number)/20))+
    ggplot2::scale_x_continuous(name = "Clonotype Size")
  if(!is.null(pheno_data) & !is.null(category)){
    plot_cls = plot_cls + ggplot2::facet_wrap(~ cat) +
      ggplot2::ggtitle(category)
  }
  plot_cls = plot_cls + ggplot2::theme_classic()+
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8, colour = "black"),
          legend.title = ggplot2::element_text(size = 9, colour = "black"),
          axis.title = ggplot2::element_text(size = 9, colour = "black"),
          axis.text = ggplot2::element_text(size = 7.7, colour = "black"),
          legend.position="right",
          legend.key.size = ggplot2::unit(0.3, "cm"),
          legend.box.spacing = ggplot2::unit(0.03, "cm"),
          plot.title = ggplot2::element_text(hjust = 0.5))

  if(plot_out){
    print(plot_cls)
  }
  return(plot_cls)
}


#
# Sankey diagramme of clonotype sharing
#


#
# Sankey diagramme of VDJ usage
#


