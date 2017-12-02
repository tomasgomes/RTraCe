#######################################################################################
#                                                                                     #
#                                                                                     #
#                                   Support Functions                                 #
#                                                                                     #
#                                                                                     #
#######################################################################################


#
# General variables
#
.fill_vals = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(12, "Set3")[-9],
               RColorBrewer::brewer.pal(8, "Pastel1"), RColorBrewer::brewer.pal(7, "Pastel2"), RColorBrewer::brewer.pal(7, "Dark2"),
               RColorBrewer::brewer.pal(12, "Paired"), RColorBrewer::brewer.pal(7, "Accent"),
               RColorBrewer::brewer.pal(11, "Spectral"))
.shade_cols = c("grey9", "grey29", "grey55", "grey78", "grey91")


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
