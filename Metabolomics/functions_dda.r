

parseScores = function(x) {
  nx = length(x)
  ss = rep(NA, nx)
  for(i in 1:nx) {
    xx = strsplit(x[i], ", ")[[1]]
    if(any(xx != "")) ss[i] = max(as.numeric(xx))   
  }
  ss
}


Detection_Freq = function(data, matrix){
  x = apply(sapply(X = c("NASH", matrix, "RT"), FUN = grepl, colnames(data)), 1, FUN = all)
  y = apply(sapply(X = c(matrix, "RT"), FUN = grepl, colnames(data)), 1, FUN = all)
  data = data[,x|y]
  
  n = ncol(data)
  freq = apply(data, 1, function(x)(length(which(x != 0))/n))
  return(freq)
}



remove_isotope = function(data){
  # identify isotope pairs ( Same RT, mass difference 1.003 Da )
  iso_tol = 0.001
  RT_tol = 0.05
  
  mz_list = data$`feature_m/z` #mz_list isotopic ion
  RT_list = data$Median
  row = rep(NA, nrow(data))
  pair_ind = rep("", nrow(data))
  n = 0
  
  for (i in 1:length(mz_list)){
    mz_diff = mz_list[i] - data$`feature_m/z`
    RT_diff = data$Median - RT_list[i]
    
    s = which(mz_diff < (1.003 + iso_tol) & mz_diff > (1.003 - iso_tol) & abs(RT_diff) < RT_tol)
    
    if (!identical(s, integer(0))){
      n = n + 1  # pair n
      row[i] = i # isotopic ion
      pair_ind[i] = paste0(paste0("isotope_", n), 
                           ifelse(pair_ind[i] != "", paste0(" | ", pair_ind[i]), ""))
      row[s] = s # parent ion
      pair_ind[s] = paste0(paste0("parent_", n), 
                           ifelse(pair_ind[s] != "", paste0(" | ", pair_ind[s]), ""))
    }
  }
  
  df = cbind(row, pair_ind)
  df = na.omit(df) # potential isotopic pairs
  
  
  # Annotate isotopic ratio : carbon * 1.1 % ?
  s = grepl("NASH", colnames(data)) & !grepl("RT_|score_|CE_", colnames(data)) 
  peak = data[,s]
  
  peak[is.na(peak)] = 0 # select a single value of peak area 
  for(k in 1:ncol(peak))peak[,k] = parseScores(as.character(peak[,k])) 
  peak[peak == 0] = NA
  
  ratio_median = isotopes = rep(NA, nrow(peak))
  pair_ind = strsplit(df[,2], " | ")
  
  for (i in 1:n){
    x = which(unlist(lapply(pair_ind, function(x)paste0("parent_", i) %in% x)))
    y = which(unlist(lapply(pair_ind, function(x)paste0("isotope_", i) %in% x)))
    
    p = as.numeric(df[x, 1])
    iso = as.numeric(df[y, 1])
    
    isotopes[p] = df[x, 2]
    isotopes[iso] = df[y, 2]
    
    d = rbind(peak[p,], peak[iso,])
    complete = apply(d, 2, function(x) !any(is.na(x)))
    
    d = data.frame(d[,complete])
    ratio = as.numeric(d[2,]) / as.numeric(d[1,]) / 0.011 
    ratio_median[iso] = median(ratio)
  }
  
  # IDs with very high ratio are less prone to be isotopes
  Non_iso = which(ratio_median > data$`feature_m/z`/12)
  x = grep("isotope", isotopes)
  x = x[!x %in% Non_iso]
  
  data = data[-x, ] # remove isotopes
  return(data)
}
