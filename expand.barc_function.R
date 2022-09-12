expand.barc <- function(dataset_TCR = bs833_d0_TCR_clono){
  clonenames <- names(sort(table(as.data.frame(dataset_TCR)[,15]),decreasing = T))
  
  ###
  clone_limbot <- clonenames[sort(table(as.data.frame(dataset_TCR)[,15]),
                                  decreasing = T)>50]
  high_exp <- clone_limbot

  ###
  clone_limbot <- clonenames[sort(table(as.data.frame(dataset_TCR)[,15]),
                                  decreasing = T)>5]
  clone_limtop <- clonenames[sort(table(as.data.frame(dataset_TCR)[,15]),
                                  decreasing = T)<51]
  low2_exp <- duplicate(c(clone_limbot,clone_limtop))
  
  ###
  clone_limbot <- clonenames[sort(table(as.data.frame(dataset_TCR)[,15]),
                                  decreasing = T)>1]
  clone_limtop <- clonenames[sort(table(as.data.frame(dataset_TCR)[,15]),
                                  decreasing = T)<6]
  low1_exp <- duplicate(c(clone_limbot,clone_limtop))
  
  ###
  clone_limbot <- clonenames[sort(table(as.data.frame(dataset_TCR)[,15]),
                                  decreasing = T)==1]
  un_exp <- clone_limbot
  
  high_exp_barcode <- dataset_TCR$barcode[dataset_TCR$clonotype_new %in% high_exp]
  low2_exp_barcode <- dataset_TCR$barcode[dataset_TCR$clonotype_new %in% low2_exp]
  low1_exp_barcode <- dataset_TCR$barcode[dataset_TCR$clonotype_new %in% low1_exp]
  un_exp_barcode <- dataset_TCR$barcode[dataset_TCR$clonotype_new %in% un_exp]

  return(list(high_exp_barcode, 
              low2_exp_barcode, low1_exp_barcode, 
              un_exp_barcode))
}
