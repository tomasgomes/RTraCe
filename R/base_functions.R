#######################################################################################
#                                                                                     #
#                                                                                     #
#                                   Base Functions                                    #
#                                                                                     #
#                                                                                     #
#######################################################################################


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
                           # is productive?
                           "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           # len and CDR3 for each chain
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No")
      }
      if(startsWith(line[2], "No seqs")){
        cells_dic[[cell]] = c("NoTCR", "NA", "NA", "NA", "NA",
                           "NA", "NA", "NA", "NA",
                           # is productive?
                           "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           # len and CDR3 for each chain
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No")
      } else if(line[2]=="A"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")
        rec_len = line[5]
        CDR3aa = line[6]
        CDR3nt = line[7]

        if(cells_dic[[cell]][2]=="No"){
          cells_dic[[cell]][2] = seq
          cells_dic[[cell]][10] = prod
          cells_dic[[cell]][18] = rec_len
          cells_dic[[cell]][19] = CDR3aa
          cells_dic[[cell]][20] = CDR3nt
        } else{
          cells_dic[[cell]][3] = seq
          cells_dic[[cell]][11] = prod
          cells_dic[[cell]][21] = rec_len
          cells_dic[[cell]][22] = CDR3aa
          cells_dic[[cell]][23] = CDR3nt
        }
        
      } else if(line[2]=="B"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")
        rec_len = line[5]
        CDR3aa = line[6]
        CDR3nt = line[7]

        if(cells_dic[[cell]][4]=="No"){
          cells_dic[[cell]][4] = seq
          cells_dic[[cell]][12] = prod
          cells_dic[[cell]][24] = rec_len
          cells_dic[[cell]][25] = CDR3aa
          cells_dic[[cell]][26] = CDR3nt
        } else{
          cells_dic[[cell]][5] = seq
          cells_dic[[cell]][13] = prod
          cells_dic[[cell]][27] = rec_len
          cells_dic[[cell]][28] = CDR3aa
          cells_dic[[cell]][29] = CDR3nt
        }
        
      } else if(line[2]=="G"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")
        rec_len = line[5]
        CDR3aa = line[6]
        CDR3nt = line[7]

        if(cells_dic[[cell]][6]=="No"){
          cells_dic[[cell]][6] = seq
          cells_dic[[cell]][14] = prod
          cells_dic[[cell]][30] = rec_len
          cells_dic[[cell]][31] = CDR3aa
          cells_dic[[cell]][32] = CDR3nt
        } else{
          cells_dic[[cell]][7] = seq
          cells_dic[[cell]][15] = prod
          cells_dic[[cell]][33] = rec_len
          cells_dic[[cell]][34] = CDR3aa
          cells_dic[[cell]][35] = CDR3nt
        }
        
      } else if(line[2]=="D"){
        seq = line[3]
        prod = ifelse(line[4]=="True", "Yes", "No")
        rec_len = line[5]
        CDR3aa = line[6]
        CDR3nt = line[7]

        if(cells_dic[[cell]][8]=="No"){
          cells_dic[[cell]][8] = seq
          cells_dic[[cell]][16] = prod
          cells_dic[[cell]][36] = rec_len
          cells_dic[[cell]][37] = CDR3aa
          cells_dic[[cell]][38] = CDR3nt
        } else{
          cells_dic[[cell]][9] = seq
          cells_dic[[cell]][17] = prod
          cells_dic[[cell]][39] = rec_len
          cells_dic[[cell]][40] = CDR3aa
          cells_dic[[cell]][41] = CDR3nt
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
                           # is productive?
                           "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           # len and CDR3 for each chain
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No", 
                           "No", "No", "No", "No", "No", "No")
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
  cl_tracer_data = tracerData$tracer_metadata[grepl("cl", tracerData$tracer_metadata$tcr_info),]

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

  tracerData$matching_chains = mat_chains
  return(tracerData)
}


#
# Obtain a list of TCR segments for each cell
#
segmentMatrix <- function(tracerData){
  chain_cols = c("tcr_A_1", "tcr_A_2", "tcr_B_1", "tcr_B_2",
                 "tcr_G_1", "tcr_G_2", "tcr_D_1", "tcr_D_2")
  prod_cols = c("tcr_pA_1", "tcr_pA_2", "tcr_pB_1", "tcr_pB_2",
                "tcr_pG_1", "tcr_pG_2", "tcr_pD_1", "tcr_pD_2")
  cdr3_cols = c("A_1_recLen", "A_1_CDR3aa", "A_1_CDR3nt", 
                "A_2_recLen", "A_2_CDR3aa", "A_2_CDR3nt", 
                "B_1_recLen", "B_1_CDR3aa", "B_1_CDR3nt", 
                "B_2_recLen", "B_2_CDR3aa", "B_2_CDR3nt", 
                "G_1_recLen", "G_1_CDR3aa", "G_1_CDR3nt", 
                "G_2_recLen", "G_2_CDR3aa", "G_2_CDR3nt", 
                "D_1_recLen", "D_1_CDR3aa", "D_1_CDR3nt", 
                "D_2_recLen", "D_2_CDR3aa", "D_2_CDR3nt")
  seg_df = tracerData$tracer_metadata[,-1]
  seg_df$cell_name = rownames(seg_df)
  seg_df = data.frame(lapply(seg_df, as.character), stringsAsFactors = F)

  # transform matrix
  seg_df_c = reshape2::melt(seg_df[,c(grep("tcr_._[12]",
                                           colnames(seg_df), value = T), "cell_name")],
                            id.vars = "cell_name")
  seg_df_p = reshape2::melt(seg_df[,c(grep("tcr_p._[12]",
                                           colnames(seg_df), value = T), "cell_name")],
                            id.vars = "cell_name")
  seg_df_len = reshape2::melt(seg_df[,c(cdr3_cols[c(1,4,7,10, 13, 16, 19, 22)], "cell_name")],
                              id.vars = "cell_name")
  seg_df_aa = reshape2::melt(seg_df[,c(cdr3_cols[c(2,5,8,11,14,17,20,23)], "cell_name")],
                              id.vars = "cell_name")
  seg_df_nt = reshape2::melt(seg_df[,c(cdr3_cols[c(3,6,9,12,15,18,21,24)], "cell_name")],
                              id.vars = "cell_name")

  # format matrix
  l_segs = strsplit(seg_df_c$value, "_")
  l_segs = lapply(l_segs, function(x) if(length(x)>3) c(paste0(x[1], "_", x[2]), x[3:4]) else x)
  l_segs = lapply(l_segs, function(x) if(length(x)<3) rep(x, 3) else x)
  l_segs = data.frame(matrix(unlist(l_segs), ncol = 3, byrow = T), stringsAsFactors = F)
  l_segs$cell_name = seg_df_c$cell_name
  colnames(l_segs)[1:3] = c("variable", "diversity_link", "joining")
  l_segs = l_segs[,c(4,1:3)]
  l_segs$chain = gsub("TR", "", substr(l_segs$variable, 1, 3))
  l_segs$chain[l_segs$chain %in% c("NA", "No")] = NA
  l_segs[l_segs=="NA"] = NA
                  
  l_segs$productive = seg_df_p$value
  l_segs$productive[rowSums(is.na(l_segs))>1] = NA
                  
  l_segs$reconstructed_length = seg_df_len$value
  l_segs$reconstructed_length[rowSums(is.na(l_segs))>1] = NA
  l_segs$CDR3aa = seg_df_aa$value
  l_segs$CDR3aa[rowSums(is.na(l_segs))>1] = NA
  l_segs$CDR3nt = seg_df_nt$value
  l_segs$CDR3nt[rowSums(is.na(l_segs))>1] = NA
                  
  l_segs = merge(l_segs, tracerData$tracer_metadata[,"tcr_info"], by.x = 1, by.y = 0)
  tracerData$tracer_metadata = tracerData$tracer_metadata[,!(colnames(tracerData$tracer_metadata) %in% cdr3_cols)]
  tracerData$vdj_segments = l_segs
  return(tracerData)
}


#
# Define Clonotypes based on specific criteria
#
defineClonotypes <- function(tracerData, criteriaList,
                             nameVar = "custom_cl"){
  chains_tcr = c("A", "B", "G", "D")
  tracer_meta = tracerData$tracer_metadata
  matChains = if(!is.null(tracerData$matching_chains)){
    tracerData$matching_chains
  } else stop("Matching chains for cells have not been defined.")

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

  tracer_meta = merge(tracer_meta, clusters, by = 0, all = T)
  rownames(tracer_meta) = tracer_meta[,1]
  tracer_meta = tracer_meta[,-1]

  # labels for cells that are not in a clonotype
  ## cells that were in a clonotype but are not anymore will be assigned as "notAssigned"
  tracer_meta[is.na(tracer_meta[,nameVar]) & grepl("cl", tracer_meta[,"tcr_info"]), nameVar] = "notAssigned"
  tracer_meta[is.na(tracer_meta[,nameVar]), nameVar] = as.character(tracer_meta[is.na(tracer_meta[,nameVar]),1])

  tracerData$tracer_metadata = tracer_meta
  return(tracerData)
}


#
# Read TraCeR results
#
readTracer <- function(summaryPath,
                       recombinants = "recombinants.txt", tcrSum = "TCR_summary.txt",
                       clonotypes = NULL, nameVar = "custom_cl",
                       get_matching_chains = F,
                       get_vdj_segments = F) {

  # Process recombinants files (in function argument) and TCR info (iNKT, clonotypes, ...)
  processed_files = processTCR(paste0(summaryPath, tcrSum),
                               processRecomb(paste0(summaryPath, recombinants)))

  # Reformat as a data frame
  tracer_data = t(data.frame(processed_files))
  colnames(tracer_data) = paste0("tcr_", c("info", "A_1", "A_2", "B_1", "B_2",
                            "G_1", "G_2", "D_1", "D_2",
                            "pA_1", "pA_2", "pB_1", "pB_2",
                            "pG_1", "pG_2", "pD_1", "pD_2",
                            "A_1_recLen", "A_1_CDR3aa", "A_1_CDR3nt", 
                            "A_2_recLen", "A_2_CDR3aa", "A_2_CDR3nt", 
                            "B_1_recLen", "B_1_CDR3aa", "B_1_CDR3nt", 
                            "B_2_recLen", "B_2_CDR3aa", "B_2_CDR3nt", 
                            "G_1_recLen", "G_1_CDR3aa", "G_1_CDR3nt", 
                            "G_2_recLen", "G_2_CDR3aa", "G_2_CDR3nt", 
                            "D_1_recLen", "D_1_CDR3aa", "D_1_CDR3nt", 
                            "D_2_recLen", "D_2_CDR3aa", "D_2_CDR3nt"))
  tracer_data = data.frame(tracer_data)
  tracer_data$tcrCode = .cellCode(tracer_data)

  # Define final object to be returned
  result = list("tracer_metadata" = tracer_data)

  # Add clonotype information according to the criteria defined
  # Examples
  # criteriaList = list(c("A"), c("B")) - only A and B productive
  # criteriaList = list(c("A", "B")) - either A or B productive
  if(get_matching_chains | !is.null(clonotypes)){
    result = matchingChainsMatrix(result)
    if(!is.null(clonotypes)){
      result = defineClonotypes(result, clonotypes,
                                nameVar = nameVar)
    }
  }

  # Add individual segment information
  if(get_vdj_segments){ result = segmentMatrix(result) }

  return(result)
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


#
# Filter clonotypes if they're shared between variable classes
#
filterClonotypes <- function(tracerData, pheno_data, category,
                             clonotypes = "tcr_info", new_col = NULL,
                             label = "_public", append = T){
  # match with pheno_data
  clon_df = data.frame(row.names = rownames(tracerData$tracer_metadata),
                       "tcr_info" = tracerData$tracer_metadata[,clonotypes])
  colnames(clon_df)[1] = clonotypes[1]
  clon_df = merge(pheno_data, clon_df, by = 0, all.x = T)
  clon_df = clon_df[,c(clonotypes, category)]
  # annotate single clonotypes as "notAssigned"
  clon_df$tcr_info = .correctClonotypes(clon_df$tcr_info)
  clon_df$tcr_info = as.character(clon_df$tcr_info)
  clon_df[,category] = as.character(clon_df[,category])
  clon_df = clon_df[grepl("cl", clon_df$tcr_info),]

  # identify clonotypes spanning more than one category
  isSharedCat = rowSums(table(clon_df)>0)>1

  # edit names
  tcr_info = as.character(tracerData$tracer_metadata[,clonotypes])
  if(append){
    tcr_info[tcr_info %in% names(isSharedCat)[isSharedCat]] = paste0(tcr_info[tcr_info %in% names(isSharedCat)[isSharedCat]],
                                                                     label)
  } else{
    tcr_info[tcr_info %in% names(isSharedCat)[isSharedCat]] = label
  }

  # include in tracer_metadata
  if(!is.null(new_col)){
    tracerData$tracer_metadata[,new_col] = tcr_info
  } else{
    warning(paste0("TraCeR metadata column ", clonotypes, " will be edited."))
    tracerData$tracer_metadata[,clonotypes] = tcr_info
  }
  return(tracerData)
}




