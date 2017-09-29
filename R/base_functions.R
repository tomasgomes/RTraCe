# Includes functions for reading and plotting

# Process recombinants file
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



# Process TCR summary file
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



# Obtain matrix of matching chains between pairs of cells reported as sharing chains
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



# Define Clonotypes based on specific criteria
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

  clusters <- igraph::clusters(igraph::graph.data.frame(mat_melt[,1:2]))$membership
  clusters = data.frame(row.names = names(clusters), nameVar = paste0("cl",clusters), stringsAsFactors = F)
  colnames(clusters) = nameVar

  tracerData = merge(tracerData, clusters, by = 0, all = T)
  rownames(tracerData) = tracerData[,1]
  tracerData = tracerData[,-1]
  tracerData[is.na(tracerData[,nameVar]),nameVar] = as.character(tracerData[is.na(tracerData[,nameVar]),1])

  return(tracerData)
}



# Read TraCeR results
readTracer <- function(summaryPath,
                       recombinants = "recombinants.txt", tcrSum = "TCR_summary.txt",
                       clonotypes = NULL, nameVar = "custom_cl") {

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
  if(!is.null(clonotypes)){
    chains_mat = matchingChainsMatrix(tracer_data)
    tracer_data = defineClonotypes(tracer_data, chains_mat,
                                   clonotypes, nameVar = nameVar)
  }

  return(tracer_data)
}



# Project clonotypes in tSNE plot



# Count clonotypes by category



# Matrix plot showing shared chains



# Sankey diagramme of clonotype sharing



# Sankey diagramme of VDJ usage



