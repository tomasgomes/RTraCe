Add AB clonotypes

```{r}
for(n in names(tpmMatrix_corrected_list)){
  info_tab = infoTable_corrected_list[[n]]
  info_tab = info_tab[info_tab$clonotypes!="notAssigned",
                      grepl("tracer", colnames(info_tab)) | colnames(info_tab) %in% c("tissue_cell",
                                                                                      "t_c_cond",
                                                                                      "clonotypes")]

  labs_chains = rep(c("A", "B", "G", "D"), each = 2)
  mat_chains = matrix("", nrow(info_tab), nrow(info_tab))
  rownames(mat_chains) = rownames(info_tab)
  colnames(mat_chains) = rownames(info_tab)
  for(i in 1:nrow(info_tab)){
    for(j in 1:nrow(info_tab)){
      for(k in 4:11){
        if(as.character(info_tab[i,k]) %in% as.character(unlist(info_tab[j,4:11])) &
           info_tab[i,k]!="No" & info_tab[i,k+8]=="Yes"){
          mat_chains[i,j] = paste0(mat_chains[i,j], labs_chains[k-3])
        }
      }
    }
  }

  mat_melt = melt(mat_chains)
  mat_melt = mat_melt[mat_melt[,1]!=mat_melt[,2] & grepl("AB", mat_melt[,3]),]

  clABlist = list()
  for(i in 1:nrow(mat_melt)){
    if(any(unlist(lapply(clABlist, function(x) mat_melt[i,1] %in% x)))){
      clABlist[[which(unlist(lapply(clABlist, function(x) mat_melt[i,1] %in% x)))]] = c(clABlist[[which(unlist(lapply(clABlist, function(x) mat_melt[i,1] %in% x)))]], as.character(unlist(mat_melt[i,1:2])))
    } else if(any(unlist(lapply(clABlist, function(x) mat_melt[i,2] %in% x)))){
      clABlist[[which(unlist(lapply(clABlist, function(x) mat_melt[i,2] %in% x)))]] = c(clABlist[[which(unlist(lapply(clABlist, function(x) mat_melt[i,2] %in% x)))]], as.character(unlist(mat_melt[i,1:2])))
    } else{
      clABlist[[length(clABlist)+1]] = as.character(unlist(mat_melt[i,1:2]))
    }
  }
  clABlist = lapply(clABlist, unique)
  names(clABlist) = paste0("clAB", 1:length(clABlist))
  clABlist = melt(clABlist)
  colnames(clABlist) = c("cell_name", "clAB")

  infoTable_corrected_list[[n]] = merge(infoTable_corrected_list[[n]], clABlist,
                                        by.x = 0, by.y = "cell_name", all.x = T)
  rownames(infoTable_corrected_list[[n]]) = infoTable_corrected_list[[n]][,1]
  infoTable_corrected_list[[n]] = infoTable_corrected_list[[n]][,-1]
  infoTable_corrected_list[[n]]$clAB[is.na(infoTable_corrected_list[[n]]$clAB)] = "notAssigned"
}
```



tSNE clonotypes

```{r}
cl_tsne_list = list()
pdf("../FINAL_PLOTS/tSNE_clonotypes_final_link_2.pdf", width = 6.5, height = 6.9, useDingbats = F)
for(dataset in names(tpmMatrix_corrected_list)){
  for(cc in c("clonotypes", "clAB")){

    tsne_df = infoTable_corrected_list[[dataset]]
    colnames(tsne_df)[which(grepl(cc, x = colnames(tsne_df)))] = "classes"
    colnames(tsne_df)[which(grepl("tissue_cell", x = colnames(tsne_df)))] = "cond"
    tsne_df$classes = factor(tsne_df$classes, levels = unique(tsne_df$classes)[order(unique(as.character(tsne_df$classes)))])
    tsne_df$cl_size = factor(ifelse(tsne_df$classes=="notAssigned", "no", "yes"), levels = c("no", "yes"))

    tsne_df$Tumour = ifelse(tsne_df$condition=="Tumour", "Tumour", "Control")

    cl_df = tsne_df[tsne_df$classes!="notAssigned", c("dim1", "dim2", "classes","cond")]

    n_cl1 = length(unique(tsne_df$classes))-1
    col_vec = c(colours_clon[1:n_cl1], "#7e7e7e")

    if(grepl("human", dataset)){
      shape_vals = c(10, 1, 19, 12, 0, 15, 11, 2, 17)
    } else{
      shape_vals = c(1, 19, 0, 15, 2, 17)
    }

    plot_tsne = ggplot(data = tsne_df, aes(x = dim1, y = dim2, shape = cond))
    if(grepl("mel", dataset)){
      plot_tsne = plot_tsne + geom_point(aes(colour = classes, fill = classes, size = Tumour), alpha = 0.8)+
        geom_line(data = cl_df, aes(x = dim1, y = dim2, colour = classes, group = classes), size = 1.06, show.legend = F)+
        scale_size_manual(values = c(1, 2.25), guide = T)
    } else{
      plot_tsne = plot_tsne + geom_point(aes(colour = classes, fill = classes, size = cl_size), alpha = 0.8)+
        geom_line(data = cl_df, aes(x = dim1, y = dim2, colour = classes, group = classes), size = 1.06, show.legend = F)+
        scale_size_manual(values = c(1, 2.02), guide = F)
    }

    plot_tsne = plot_tsne +
      #ggtitle(paste0("t-SNE ", dataset, " n_cells: ", nrow(tsne_df)))+
      scale_shape_manual(values = shape_vals)+
      guides(shape = guide_legend(nrow = 3))+
      scale_colour_manual(values = col_vec)+
      theme_tsne+
      theme(legend.text = element_text(size = 6),
            legend.title = element_text(size = 7),
            axis.title = element_text(size = 7),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position="right",
            legend.key.size = unit(0.3, "cm"),
            legend.box.spacing = unit(0.03, "cm"))

    if(F){
      if(!grepl("_mel", dataset)) plot_tsne = plot_tsne + guides(colour = guide_legend(ncol = 4), size = "none")
      if(grepl("_mel", dataset)) plot_tsne = plot_tsne + guides(colour = guide_legend(nrow = 3), size = guide_legend(nrow = 3))
      if(!grepl("human", dataset)) plot_tsne = plot_tsne + guides(shape = guide_legend(ncol = 2, byrow = T))
    } else{
      if(grepl("_mel", dataset)){
        plot_tsne = plot_tsne + guides(colour = "none", size = guide_legend(nrow = 3), fill = "none",
                                       shape = guide_legend(title = "Tissue/Cell", ncol = 1, byrow = T))
      } else if(grepl("human", dataset)){
        plot_tsne = plot_tsne + guides(shape = guide_legend(title = "Tissue/Cell", ncol = 1, byrow = T),
                                       colour = "none", fill = "none", size = "none")
      } else{
        plot_tsne = plot_tsne + guides(colour = "none", size = "none", fill = "none",
                                       shape = guide_legend(title = "Tissue/Cell", ncol = 1, byrow = T))
      }
    }

    cl_tsne_list[[paste0(dataset, "_", cc)]] = plot_tsne
    print(plot_tsne)
  }
}
dev.off()
save(cl_tsne_list, file = "./plot_repository/cl_tsne_list_2.RData")
```



Table for clonotype correspondence

```{r}
cl_tab_list = list()
for(n in names(infoTable_corrected_list)){
  for(cc in c("clonotypes", "clAB")){
    info_tab = infoTable_corrected_list[[n]]
    info_tab = info_tab[info_tab[,cc]!="notAssigned",
                        grepl("tracer", colnames(info_tab)) | colnames(info_tab) %in% c("tissue_cell","t_c_cond", cc)]

    labs_chains = rep(c("A", "B", "G", "D"), each = 2)
    mat_chains = matrix("", nrow(info_tab), nrow(info_tab))
    rownames(mat_chains) = rownames(info_tab)
    colnames(mat_chains) = rownames(info_tab)
    for(i in 1:nrow(info_tab)){
      for(j in 1:nrow(info_tab)){
        for(k in 4:11){
          if(as.character(info_tab[i,k]) %in% as.character(unlist(info_tab[j,4:11])) &
             info_tab[i,k]!="No"){
            mat_chains[i,j] = paste0(mat_chains[i,j], labs_chains[k-3])
          }
        }
      }
    }
    n_chains = matrix(sapply(mat_chains, nchar), nrow(info_tab), nrow(info_tab))
    rownames(n_chains) = rownames(info_tab)
    colnames(n_chains) = rownames(info_tab)

    annot_df = data.frame(row.names = rownames(info_tab),
                          "TissueCell" = info_tab$tissue_cell)
    annot_cl = data.frame(row.names = rownames(info_tab),
                          "Clonotype" = info_tab[,cc])

    col_vec = if(grepl("mouse", n)) colours_mouse else colours_h_tc
    names(col_vec) = levels(infoTable_corrected_list[[n]]$tissue_cell)

    n_cl1 = length(unique(info_tab[,cc]))
    cl_col = colours_clon[1:n_cl1]

    names(cl_col) = levels(factor(info_tab[,cc]))

    ann_col = list("TissueCell" = col_vec,
                   "Clonotype" = cl_col)

    if(n=="mouse_colon"){
      #gaps = cumsum(c(2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,2,3))
      gaps = length(cl_col)

      hc = hclust(dist(n_chains))

      cl_tab_list[[paste0(n, "_", cc)]] = pheatmap(n_chains[hc$order,hc$order], show_rownames = F, show_colnames = F,
                                                   annotation_row = annot_cl, annotation_col = annot_df,
                                                   display_numbers = mat_chains[hc$order,hc$order],
                                                   annotation_names_col = F, annotation_names_row = F,
                                                   cluster_rows = F, cluster_cols = F,
                                                   annotation_colors = ann_col, color = c("white", viridis::viridis(200)[76:200]),
                                                   treeheight_row = 0, treeheight_col = 0, gaps_col = gaps,
                                                   clustering_method = "ward.D2", fontsize = 5, fontsize_number = 3.5)$gtable

    } else{
      gaps = length(cl_col)

      cl_tab_list[[paste0(n, "_", cc)]] = pheatmap(n_chains, show_rownames = F, show_colnames = F,
                                                   annotation_row = annot_cl, annotation_col = annot_df, display_numbers = mat_chains,
                                                   annotation_names_col = F, annotation_names_row = F,
                                                   annotation_colors = ann_col, color = c("white", viridis::viridis(200)[76:200]),
                                                   treeheight_row = 0, treeheight_col = 0, cutree_cols = gaps,
                                                   clustering_method = "ward.D2", fontsize = 5, fontsize_number = 3.5)$gtable
    }

  }
}
save(cl_tab_list, file = "./plot_repository/cl_tab_list.RData")
```



Counts within and between tissues

```{r, fig.width=4, fig.height=4}
cl_tab_match_list = list()
for(n in names(infoTable_corrected_list)){
  matrix_list = list()
  for(cc in c("clonotypes", "clAB")){
    tsne_df = infoTable_corrected_list[[n]]
    colnames(tsne_df)[which(grepl(cc, x = colnames(tsne_df)))] = "classes"
    colnames(tsne_df)[which(grepl("tissue_cell", x = colnames(tsne_df)))] = "cond"
    tsne_df$classes = factor(tsne_df$classes, levels = unique(tsne_df$classes)[order(unique(as.character(tsne_df$classes)))])

    tsne_df = tsne_df[tsne_df$classes!="notAssigned",]

    matrix_counts = matrix(rep(0, length(levels(tsne_df$cond))**2),
                           length(levels(tsne_df$cond)), length(levels(tsne_df$cond)))
    rownames(matrix_counts) = levels(tsne_df$cond)
    colnames(matrix_counts) = levels(tsne_df$cond)
    for(cl in unique(tsne_df$classes)){
      sub_tsne_df = tsne_df[tsne_df$classes==cl,c("cond")]

      x = t(combn(sub_tsne_df, 2))
      x.sort = t(apply(x, 1, sort))
      cor_mat = x[!duplicated(x.sort),]
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
    matrix_list[[cc]] = matrix_counts
  }

  intra_pop_cl = data.frame("All" = diag(matrix_list[[1]]), "AB" = diag(matrix_list[[2]]))
  inter_pop_cl = matrix_list[[1]]
  gdata::lowerTriangle(inter_pop_cl) = gdata::lowerTriangle(matrix_list[[2]])
  diag(inter_pop_cl) = rep(max(inter_pop_cl)+1, length(diag(inter_pop_cl)))
  lab_mat = matrix(as.character(inter_pop_cl), nrow(inter_pop_cl), ncol(inter_pop_cl))
  diag(lab_mat) = ""

  diag(inter_pop_cl) = -1
  cl_tab_match_list[[paste0(n, "_inter")]] = pheatmap(inter_pop_cl,  show_rownames = T, show_colnames = T, display_numbers = lab_mat,
                                                      cluster_rows = F, cluster_cols = F, annotation_colors = ann_col,
                                                      color = c("white", viridis::viridis(200)[76:200]), fontsize = 7.5, fontsize_number = 6.5,
                                                      legend = F, number_color = "black", border_color = "white")$gtable
  cl_tab_match_list[[paste0(n, "_intra")]] = pheatmap(intra_pop_cl, show_rownames = T, show_colnames = T, display_numbers = round(intra_pop_cl, 0),
                                                      cluster_rows = F, cluster_cols = F, annotation_colors = ann_col,
                                                      color = c(viridis::viridis(200)[76:200]), fontsize = 7.5, fontsize_number = 6.5,
                                                      legend = F, number_color = "black", border_color = "white")$gtable
}
save(cl_tab_match_list, file = "./plot_repository/cl_tab_match_list.RData")
```
