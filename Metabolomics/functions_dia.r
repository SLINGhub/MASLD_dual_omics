library(stringr)
library(stringi)
library(gsubfn)

# Get QC quant data for the matrix
QC_area = function(matrix, data){
  x = apply(sapply(X = c("NASH", "_QC", matrix), FUN = grepl, colnames(data)), 1, FUN = all)
  y = !grepl("RT_|score_", colnames(data)) 
  return(x & y)
}

# normalize QC 
normalizeQC = function(x){
  sum = apply(x, 2, function(x) sum(x, na.rm=TRUE))
  x2 = x
  for(j in 1:ncol(x2)) {
    x2[,j] = x2[,j]/sum[j]
  }
  x2
}

# Calculate CoV
CoV_QC = function(data){
  nor_d = normalizeQC(data)
  CoV = rep(NA, nrow(data))
  for(i in 1:nrow(data)){
    x = as.numeric(nor_d[i,])
    xx = x[!is.na(x)]
    if(length(xx) >= 3) {
      CoV[i] = sd(xx)/mean(xx)*100
    }
  }
  return(CoV)
}


Sample_area = function(matrix, data){
  x = apply(sapply(X = c("NASH",  matrix), FUN = grepl, colnames(data)), 1, FUN = all)
  y = !grepl("RT_|score_|_QC", colnames(data)) 
  return(x & y)
}



add_annot = function(data, tolerance){
  library(data.table)
  library(tibble)
  library(dplyr)
  
  a = data$RT
  a_ = sort(a, index.return = T)
  diff = abs(shift(a_$x) - a_$x)
  
  group = list()
  n = 1
  group[[n]] = a_$ix[1]
  
  if(length(diff) > 1){
    for (i in 2:length(diff)){
      if (diff[i] < tolerance){
        group[[n]] = c(a_$ix[i], unlist(group[n]))
      } else {
        n = n+1
        group[[n]] = a_$ix[i]
      }
    }
  }
  
  library(tibble)
  data = add_column(data, note = rep(NA, nrow(data)), .before = "name")
  x = unlist(lapply(group, function(x)length(x) > 1))
  
  if(any(x)){
    num = 0
    for (i in 1:length(group)){
      ind = unlist(group[[i]])
      if (length(ind) == 1){
        data$note[ind] = "remove"
      } else if (length(ind) > 1 & (sum(x) < 2)){
        data$note[ind] = "DiffAdduct"
      } else if (length(ind) > 1 & (sum(x) >= 2)){
        num = num + 1
        data$note[ind] = paste0("DiffAdduct", num)
      }
    }
  }
  mul_RT = is.na(data$note) &
    (duplicated(data$Compound_name) | duplicated(data$Compound_name, fromLast = T))
  data$note = ifelse(mul_RT, "Multi_RT", data$note)
  return(data)
}

Detection_Freq = function(data, matrix){
  x = apply(sapply(X = c("NASH", matrix, "RT"), FUN = grepl, colnames(data)), 1, FUN = all)
  y = !grepl("QC_", colnames(data))
  data = data[,x & y]
  
  n = ncol(data)
  freq = apply(data, 1, function(x)(length(which(x != 0))/n * 100))
  return(freq)
}



select_fragment = function(data){
  ind = unique(data$Var.2)
  nr = 1:nrow(data)
  nr_new = NULL
  
  # 1) remove the fragment ions with frequencies < 0.9 of the precursor ion's freq
  for (i in 1:length(ind)){
    cpd = which(data$Var.2 == ind[i])
    d = data[cpd,]
    ms1 = which(d$MS_level == "1")
    ms2 = which(d$MS_level == "2")
    x = which(d$detected_all[ms2] >= d$detected_all[ms1] * 0.9)
    nr_new = append(nr_new, c(cpd[ms1], cpd[ms2][x]))
  }
  data = data[nr_new,]
  
  # 2) filter fragment ions by criteria "correlation > 0.8" 
  tmp = data[,grep("score_", colnames(data))]
  Corr_median = apply(tmp, 1, function(x)median(x, na.rm = T))
  f = which(data$MS_level == "2")
  Corr_median[f] = sapply(Corr_median[f], function(x)ifelse(x > 0.8, x, NA))
  
  data = data[!is.na(Corr_median),]
  return(data)
}


quantifier = function(data){
  ind = unique(data$Var.2)
  nr = 1:nrow(data)
  nr_new = NULL
  
  for (i in 1:length(ind)){
    cpd = which(data$Var.2 == ind[i])
    d = data[cpd,]
    
    freq = grep("detected_liver|detected_plasma|detected_cm", colnames(d))
    
    # If there are equal values, 
    # "which.max" will take the first one (priority: liver > plasma > cm)
    m = strsplit(names(which.max(d[1, freq])), "_")[[1]][2]
    cov = apply(sapply(X = c("CoV",  m), FUN = grepl, colnames(d)), 1, all)
    nr_new = append(nr_new, cpd[which.min(d[,cov])])
  }
  data = data[nr_new,]
  return(data)
}


name_split = function(string, adduct){
  Adduct = unique(unlist(strsplit(names(table(adduct)), ", ")))
  ion = paste0(" ", paste0(Adduct, collapse = "| "))
  ion = gsub("+", "_", ion, fixed =T)
  n = strsplit(as.character(string), " --- |possibly an ISF --- ")    
  nn = unlist(n)
  
  for (i in 1:length(nn)) {  
    nn[i] = gsub("+", "_", nn[i], fixed =T)
    nn[i] = strsplit(nn[i], ion)[[1]][1]
    nn[i] = strsplit(nn[i], "; | cation| Cation")[[1]][1]  # "; " for annotations in MassBank
  }
  string_new = paste0(na.omit(unique(nn)), collapse = " --- ")
  return(string_new)
}

select_name = function(data){
  library(tibble)
  library(dplyr)
  library(stringr)
  
  name = sapply(data$name, function(x)name_split(x, data$adduct))
  f = rep(NA, nrow(data))
  
  for (i in 1:nrow(data)){
    name_in_lib = as.character(name[i])
    name_hmdb = data[i, "name (HMDB)"]
    abbrev = data[i, "abbrev (LMSD)"]
    if(!is.na(abbrev)){
      f[i] = abbrev
    } #else if (tolower(name_in_lib) %in% tolower(name_hmdb)) {
    #f[i] = name_hmdb
    #} #else if (!is.na(name_hmdb) & (nchar(name_in_lib) > nchar(name_hmdb))){
    #f[i] = name_hmdb
    #} 
    else {
      f[i] = name_in_lib
    }
  }
  
  data = add_column(data, Compound_name = f, .before = "name") %>% arrange(Compound_name)
  return(data)
}


Lnames_convert = function(string){
  x = gsub("\\(.*?\\)", "", string) 
  x = gsub("-SN\\d+", "", x)
  matched_colon = gregexpr(":", x, fixed = TRUE)
  n_colon = length(unlist(matched_colon))
  
  if (grepl("ST |FAHFA ", string)){
    x_new = x
  } else if ((n_colon == "1") & !(str_detect(x, "([0-9]+)O"))) {
    x_new = x
  } else{
    x2 = str_split(x, " ")
    species = x2[[1]][1]
    unit = str_split(x2[[1]][2], "_|;|O-|N-|/")
    
    ## remove unit without colon
    unit = lapply(unit, function(x)x[grep(":",x)])
    
    Cdb = str_extract_all(unlist(unit), "[^:]+")
    carbon = sum(as.numeric(lapply(Cdb, function(x)x[1])))
    db = sum(as.numeric(lapply(Cdb, function(x)x[2])))
    O = ifelse(str_detect(x, "O-"), "O-", "")
    N = ifelse(str_detect(x, "N-"), "N-", "")
    
    tmp = gsub(";O", ";1O", x2)
    #convert ";2O" into ";O2"
    if (str_detect(tmp, "([0-9]+)O")){
      tmp = unlist(str_extract_all(tmp, "([0-9]+)O"))
      n = sum(as.integer(gsub("O", "", tmp)))
      n = ifelse(n == 1, "", n)
      nO = paste0(";O", n)
      if (str_detect(string, "OH\\)")){
        nO = gsubfn("\\d+", function(x) as.numeric(x) + 1, nO)
      } 
    } else nO = ""
    x_new = paste0(species," ", O, N, carbon,":",db, nO)
  }
  
  # Add (FA xx) or (O-xx) chain
  if(str_detect(string, "\\(FA|\\(O-")){
    y = unlist(strsplit(string, "\\(|\\)"))
    y = y[grep("FA |O-", y)]
    
    if(length(y) > 0){
      y2 = strsplit(y, " |O-")[[1]][2]
      carbon2 = as.numeric(strsplit(y2, ":")[[1]][1])
      db2 = as.numeric(strsplit(y2, ":")[[1]][2])
      carbon = carbon + carbon2
      db = db + db2
      x_new = paste0(species," ", O, N, carbon,":",db, nO)
    }
  }
  return(x_new)
}

Simplify_lip = function(string){
  n = str_split(string, " --- ")    
  nn = unlist(n)
  
  for (i in 1:length(nn)){
    if (grepl(":", nn[i])){
      nn[i] = Lnames_convert(nn[i])
    }
  }
  nn = unique(nn)
  string_new = paste0(na.omit(nn), collapse = " --- ")
  return(string_new)
}



add_annot = function(data, tolerance){
  library(data.table)
  library(tibble)
  library(dplyr)
  
  a = data$RT
  a_ = sort(a, index.return = T)
  diff = abs(shift(a_$x) - a_$x)
  
  group = list()
  n = 1
  group[[n]] = a_$ix[1]
  
  if(length(diff) > 1){
    for (i in 2:length(diff)){
      if (diff[i] < tolerance){
        group[[n]] = c(a_$ix[i], unlist(group[n]))
      } else {
        n = n+1
        group[[n]] = a_$ix[i]
      }
    }
  }
  
  library(tibble)
  data = add_column(data, note = rep(NA, nrow(data)), .before = "name")
  x = unlist(lapply(group, function(x)length(x) > 1))
  
  if(any(x)){
    num = 0
    for (i in 1:length(group)){
      ind = unlist(group[[i]])
      if (length(ind) == 1){
        data$note[ind] = "remove"
      } else if (length(ind) > 1 & (sum(x) < 2)){
        data$note[ind] = "DiffAdduct"
      } else if (length(ind) > 1 & (sum(x) >= 2)){
        num = num + 1
        data$note[ind] = paste0("DiffAdduct", num)
      }
    }
  }
  mul_RT = is.na(data$note) &
    (duplicated(data$Compound_name) | duplicated(data$Compound_name, fromLast = T))
  data$note = ifelse(mul_RT, "Multi_RT", data$note)
  return(data)
}

remove_ms1_diffRT = function(data){
  id = unique(data$Compound_name)
  rid = NULL
  
  for (i in 1:length(id)){
    ind = which(data$Compound_name == id[i])
    if (length(ind > 1)){
      conf = data[ind, "Confidence_level"]
      multRT = data[ind, "note"]
      if (grepl("MSMS", conf) & grepl("MS1", conf) & grepl("Multi_RT", multRT)){
        rid = append(rid, ind[which(conf == "MS1")])
      }
    }
  }
  data = data[-rid,]
  return(data)
}


concat_adduct = function(data){
  library(stringr)
  
  concat = data.frame(matrix(NA, nrow = 1, ncol = ncol(data)))
  colnames(concat) = colnames(data)
  
  # Different operation for columns
  # 1. CoV"  1e7 for repeating CoV calculation
  col1 = which(colnames(data) %in% c("CoV_liver", "CoV_plasma", "CoV_cm"))
  # 2. detected_x/nonzero_rate, max
  col2 = c(grep("detected_", colnames(data)), grep("nonzero", colnames(data)))
  # 3. Samples: Sum up
  col3_l = which(grepl("Liver_", colnames(data)) & (str_detect(colnames(data),"[0-9]{3}") | grepl("QC", colnames(data))))
  col3_p = which(grepl("Plasma_", colnames(data)) & (str_detect(colnames(data),"[0-9]{3}") | grepl("QC", colnames(data))))
  col3_c = which(grepl("CM_", colnames(data)) & (str_detect(colnames(data),"[0-9]{3}") | grepl("QC", colnames(data))))
  # 4. RT
  col4 = grep("RT", colnames(data))
  # 5. information: paste with collapse", ", keep unique ones
  col5 = which(! 1:ncol(data) %in% c(col1, col2, col3_l, col3_p, col3_c, col4))
  
  concat[1, col1] = 1e7
  concat[1, col2] = apply(data[, col2, drop = FALSE], 2, max, na.rm = TRUE)
  concat[1, col4] = apply(data[, col4, drop = FALSE], 2, median, na.rm = TRUE)
  concat[1, col5] = apply(data[, col5, drop = FALSE], 2, function(x)ifelse(all(is.na(x)), NA, paste(unique(x[!is.na(x)]), collapse = ", ")))
  concat[1, "note"] = NA  # remove annotation
  
  r = which(data$CoV_liver < 30 & data$detected_liver > 30)
  concat[1, "CoV_liver"]  = ifelse(length(r) == 1, data$CoV_liver[r], 1e7) #CoV"  1e7 for repeating CoV calculation
  concat[1, col3_l] = apply(data[r, col3_l], 2, function(x)ifelse(length(r) >= 1, exp(mean(log(x))), NA))
  
  r = which(data$CoV_plasma < 30 & data$detected_plasma > 30)
  concat[1, "CoV_plasma"]  = ifelse(length(r) == 1, data$CoV_plasma[r], 1e7) 
  concat[1, col3_p] = apply(data[r, col3_p], 2, function(x)ifelse(length(r) >= 1, exp(mean(log(x))), NA))
  
  r = which(data$CoV_cm < 30 & data$detected_cm > 30)
  concat[1, "CoV_cm"]  = ifelse(length(r) == 1, data$CoV_cm[r], 1e7) 
  concat[1, col3_c] = apply(data[r, col3_c], 2, function(x)ifelse(length(r) >= 1, exp(mean(log(x))), NA))
  
  return(concat)
}




