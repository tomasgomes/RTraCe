# Includes functions for reading and plotting

# import packages
require(scater)
require(ggplot2)
require(googleVis)


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
        cells_dic[cell] = c("NoTCR", "NA", "NA", "NA", "NA",
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
  while(TRUE){
    line = readLines(con, n = 1)
    if(length(line) == 0){
      break
    }
    bool_vec = c("multi" = F, "iNKT" = F, "MAIT" = F, "clon" = F)
    cl_n = 0
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

    if(startsWith("###", line)){

      cell = gsub(" ###\n", "", gsub("### ", "", line))
      if(!cell %in% names(cells_dic)){
        cells_dic[[cell]] = c("notAssigned", "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No",
                           "No", "No", "No", "No")
      }
      if(bool_vec["multi"]){
        cells_dic[[cell]][0] = "Multi_recomb"
      }else if(bool_vec["iNKT"]){
        cells_dic[[cell]][0] = "iNKT"
      }else if(bool_vec["MAIT"]){
        cells_dic[[cell]][0] = "MAIT"
      }
    }
    if(bool_vec["clon"] & grepl(",", line)){
      cl_n = cl_n + 1
      line = unlist(strsplit(gsub("\n", "", line), ","))
      for(cell in line){
        cells_dic[[cell]][0] = paste0("cl", cl_n)
      }
    }
  }

  close(con)
  return(cells_dic)
}


# Read TraCeR results
readTracer <- function(summaryPath,
                       recombinants = "recombinants.txt", tcrSum = "TCR_summary.txt") {

  processed_files = processTCR(paste(summary_path, tcrSum),
                               processRecomb(paste(summary_path, recombinants)))

  tracer_data = t(data.frame(processed_files))
  colnames(tracer_data) = c("clonotype_simp", "A_1", "A_2", "B_1", "B_2",
                            "G_1", "G_2", "D_1", "D_2",
                            "pA_1", "pA_2", "pB_1", "pB_2",
                            "pG_1", "pG_2", "pD_1", "pD_2")
  tracer_data = data.frame(tracer_data)

}


# Project clonotypes in tSNE plot


# Count clonotypes by category


# Matrix plot showing shared chains


# Sankey diagramme of clonotype sharing


# Sankey diagramme of VDJ usage



