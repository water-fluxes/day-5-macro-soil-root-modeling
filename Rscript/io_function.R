
message("loading libraries")
library(readxl)
source("./Rscript/getSUF.R") # MARSHAL
library(data.table)
library(Matrix)
library(ggplot2)
suppressWarnings(suppressMessages(library(dplyr)))
library(stringr)

message("loading function")

plot_rs<- function(all_roots){
  pl <- all_roots%>%
              ggplot()+
              geom_segment(aes(x = x1, xend = x2, y = z1, yend = z2, colour = age))+coord_fixed()+theme_classic()
  print(pl)
}


# rewrite the output file of CPlantBox
root_transform <- function(all_roots){
  all_roots <- as.data.table(all_roots)
  all_roots <- all_roots%>%
    transmute(node1ID = node1ID,
              node2ID = node2ID,
              branchID = branchID,
              x1 = x1, y1 = y1, z1 = z1, x2 = x2, y2= y2, z2 = z2,
              radius = radius,
              length = sqrt((all_roots$x2 - all_roots$x1)^2 + (all_roots$y2 - all_roots$y1)^2 + (all_roots$z2 - all_roots$z1)^2),
              time = time,
              type = type,
              age = age)%>%
    arrange(time)
  
  all_roots$node2ID <- 1:nrow(all_roots)
  all_roots$node1ID[all_roots$branchID == 1][1] <- 0
  all_roots$node1ID[all_roots$branchID == 1][-1] <- which(all_roots$branchID == 1)[-length(which(all_roots$branchID == 1))] # tap root ordination
  
  for(i in unique(all_roots$branchID)[-1]){
    all_roots$node1ID[all_roots$branchID == i][-1] <- which(all_roots$branchID == i)[-length(which(all_roots$branchID == 1))]
    if(all_roots$type[all_roots$branchID == i][1] %in% c(4,5)){ # connection with the collar
      all_roots$node1ID[all_roots$branchID == i][1] <- 0
    }
    if(all_roots$type[all_roots$branchID == i][1] %in% c(2,3)){ # connection with the parental root
      x1_child <- all_roots$x1[all_roots$branchID == i][1]
      y1_child <- all_roots$y1[all_roots$branchID == i][1]
      z1_child <- all_roots$z1[all_roots$branchID == i][1]
      
      tmp_time <- all_roots$time[all_roots$branchID == i][1]
      
      nearest <- all_roots%>%filter(branchID != i)%>%
        mutate(euc = sqrt((x1-x1_child)^2+ (y1 - y1_child)^2 + (z1 - z1_child)^2))
      nearest <- nearest[nearest$euc == min(nearest$euc), ]
      all_roots$node1ID[all_roots$branchID == i][1] <- nearest$node2ID[1] # oldest segments
    }
  }
  return(all_roots)
}


soil_initial <- function(soil_type = stype, field_capacity = 0.5){
  if(soil_type == "sandy loam"){
    soil <- data.frame(id=1:101,
                       z = sort(seq(-100,0,1), decreasing = T),
                       value = 1,
                       psi = rep(-73.5021, 101))
  }else if (soil_type == "loam") {
    soil <- data.frame(id=1:101,
                       z = sort(seq(-100,0,1), decreasing = T),
                       value = 1,
                       psi = rep(-132.869, 101))
    
  }else if (soil_type == "silty clay") {
    soil <- data.frame(id=1:101,
                       z = sort(seq(-100,0,1), decreasing = T),
                       value = 1,
                       psi = rep(-354.424, 101))
  }else{
    warning("This soil type is not referenced yet")
  }
  
  soil$psi <- soil$psi*(1/field_capacity)
  
  return(soil)
  
}


write.options.in <- function(project.path = "./Hydrus_1D_Couvreur/CouvreurV2", krs, kcomp){
  
  path <- paste0(project.path, "/Options.in")
  # path = "./Hydrus_1D_Couvreur/CouvreurV2/Options.IN"
  op = suppressWarnings(readLines(con = paste0(project.path,"/Options_1.IN"), n = -1L, encoding = "unknown"))
  tmp_krs <- op[10]
  op[10]<- paste(replace(str_split(tmp_krs, " ")[[1]], which(str_split(tmp_krs, " ")[[1]]== "0.0001"), krs), collapse = " ")
  tmp_kcomp <- op[13]
  op[13] <- paste(replace(str_split(tmp_kcomp, " ")[[1]], which(str_split(tmp_kcomp, " ")[[1]]== "0.0001"), kcomp), collapse = " ")
  
  write(op, path)
}


write.profile.dat <- function(project.path = "./Hydrus/CouvreurV2", SSF){
  
  out.file <- "PROFILE.DAT"
  prof = suppressWarnings(readLines(con = paste0(project.path,"/PROFILE_1.DAT"), n = -1L, encoding = "unknown"))
  h_ind = grep("h", prof)
  n <- as.numeric(str_split(prof[h_ind], " ")[[1]][3])
  if(nrow(SSF) != n){error("sum up SSf every two centimeter along the profile")}
  ile = data.table::fread(input = file.path(paste0(project.path,"/PROFILE_1.DAT")),
                          fill = TRUE, 
                          blank.lines.skip = FALSE, 
                          skip = 4)
  nopcol <- rep("na",(length(ile[h_ind+1,])-9))
  for(i in 1:length(nopcol)){
    nopcol[i] <- paste0("nop",i)
  }
  colnames(ile) <- c("id", "z", "h", "Mat", "Lay", "Beta", "Axz", "Bxz", "Dxz", nopcol)
  
  ile <- ile%>%select(-starts_with("nop"))
  ile <- ile[1:101, ]
  ile$h <- SSF$h
  ile$Beta <- SSF$suf
  data_fmt <- ile
  data_fmt = apply(data_fmt, MARGIN = 1, FUN = paste0, collapse = " ")
  
  prof_input1 = prof[1:h_ind]
  prof_input2 = data_fmt
  prof_input3 = prof[(1+h_ind+length(prof_input2)):length(prof)]
  
  prof_input_new = c(prof_input1, prof_input2, prof_input3)
  prof_in_file = file.path(project.path, out.file)
  write(prof_input_new, file = prof_in_file, append = F)
}

write.selector.in<- function(project.path, mesh = NA, soilparam = NA){
  
  out.file = "SELECTOR.IN"
  # default.filename = "ATMOSPH.IN"
  selector_data = readLines(con = paste0(project.path,"/SELECTOR_1.IN"), n = -1L, encoding = "unknown")
  
  if(file.exists(file.path(project.path, out.file))){
    file.remove(file.path(project.path, out.file))
  }
  if(!is.na(soilparam)){
    
    selector_data[27]<- paste(as.character(soilparam[,c("Q_r","Q_s","alpha","n","Ks","l")]),collapse = " ")
  }
  if(!is.na(mesh)){
    selector_data[30]
    mesh_1 <- selector_data[30]
    selector_data[30]<- paste(replace(str_split(mesh_1, " ")[[1]], which(str_split(mesh_1, " ")[[1]]== "1"), mesh), collapse = " ")
    selector_data[30]<- paste(replace(str_split(mesh_1, " ")[[1]], which(str_split(mesh_1, " ")[[1]]== "2"), mesh), collapse = " ")
    selector_data[30]<- paste(replace(str_split(mesh_1, " ")[[1]], which(str_split(mesh_1, " ")[[1]]== "3"), mesh), collapse = " ")
  }
  selector_in_file = file.path(project.path, out.file)
  write(selector_data, file = selector_in_file, append = F)
  
}

write.atmosph.in<- function(project.path, maxAL, deltaT, atm.bc.data, hCritS = 0, ...){
  
  out.file = "ATMOSPH.IN"
  # default.filename = "ATMOSPH.IN"
  atm_data = readLines(con = paste0(project.path,"ATMOSPH_1.IN"), n = -1L, encoding = "unknown")
  
  if(file.exists(file.path(project.path, out.file))){
    file.remove(file.path(project.path, out.file))
  }
  extinction_ind = grep("Extinction", atm_data)
  
  
  # write(atm_data, file = "ATMOSPH_IN.BAK", append = F)
  
  hcrits_ind = grep("hCritS", atm_data)
  atm_data[hcrits_ind + 1] = sprintf("%7.0f", hCritS)
  
  maxAL_ind = grep("MaxAL", atm_data)
  tAtm_ind = grep(" tAtm", atm_data)
  
  #tMax = maxAL*deltaT
  
  atm_data[(maxAL_ind + 1)] = sprintf("%7.0f", maxAL)
  end_line = atm_data[grep("end", atm_data)]
  
  # bc_data = atm_data[(tAtm_ind +1): (end_line - 1)]
  # data_ind = (tMax*(sim_ind-1) + 1):(sim_ind*tMax)
  
  # tAtm = seq(deltaT, tMax, by = deltaT)
  
  bc_data_vars = c("tAtm", "Prec", "rSoil", "rRoot", "hCritA", "rB",
                   "hB", "ht", "RootDepth")
  
  bc_data_new = atm.bc.data[1:maxAL, bc_data_vars]
  # bc_data_new = data.frame(tAtm = seq(deltaT, tMax, deltaT), bc_data_new, row.names = NULL)
  #  bc_data_new = bc_data_new[rep(seq_len(nrow(bc_data_new)), each = 4), ]
  #  bc_data_new$tAtm = seq(deltaT, tMax, by = deltaT)
  row.names(bc_data_new) = NULL
  
  tstep_decimals = get.decimalplaces(deltaT)
  
  fmt_vec = c("%12.4f", "%12.3f", "%12.4f", "%12.4f", "%12.0f", rep("%12.4f",8))
  fmt_vec[1] = sub(pattern = "0", replacement = tstep_decimals, fmt_vec[1])
  
  bc_data_fmt = bc_data_new
  
  for(a in 1:nrow(bc_data_fmt)) {
    bc_data_fmt[a, ] = sprintf(fmt = fmt_vec[1:ncol(bc_data_fmt)], bc_data_new[a, ])
  }
  bc_data_fmt = apply(bc_data_fmt, MARGIN = 1, FUN = paste, collapse = "")
  
  atm_input1 = atm_data[1:tAtm_ind]
  atm_input2 = bc_data_fmt
  atm_input3 = end_line
  
  atmosph_input_new = c(atm_input1, atm_input2, atm_input3)
  atmosph_in_file = file.path(project.path, out.file)
  write(atmosph_input_new, file = atmosph_in_file, append = F)
  
}

get.decimalplaces <- function (x) {
  if ((x%%1) != 0) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", 
                   fixed = TRUE)[[1]][[2]])
  }
  else {
    return(0)
  }
}


read.nod_inf <- function (project.path, out.file = "Nod_Inf.out", output = NULL, 
                          warn = FALSE, ...) 
{
  if (is.null(output) | missing(output)) {
    output = c("Head", "Moisture", "K", "C", "Flux", "Sink", 
               "Kappa", "v/KsTop", "Temp")
  }
  options(warn = -1)
  if (warn == TRUE) 
    options(warn = 0)
  nod_inf = data.table::fread(input = file.path(project.path, 
                                                out.file), fill = TRUE, blank.lines.skip = FALSE, skip = 9)
  time_lines = nod_inf[grepl("Time:", nod_inf[["Node"]]), 
                       ]
  times = c(0, as.numeric(time_lines$Depth))
  for (col in colnames(nod_inf)) set(nod_inf, j = col, value = as.numeric(nod_inf[[col]]))
  nod_inf = na.omit(nod_inf)
  nodes = sort(unique(nod_inf[["Node"]]))
  nod_inf[, `:=`(Time, rep(times, each = length(nodes)))]
  nod_split = split(nod_inf, f = nod_inf$Time)
  nrow_split = sapply(nod_split, nrow)
  extra_index = which(nrow_split > length(nodes))
  for (i in extra_index) {
    nod_split[[i]] = nod_split[[i]][1:length(nodes), ]
  }
  nod_inf = rbindlist(nod_split)
  output_names = intersect(output, colnames(nod_inf))
  output_names = c("Time", "Node", "Depth", output_names)
  nod_out = nod_inf[, .SD, .SDcols = output_names]
  options(warn = 0)
  return(nod_out)
}

add_hydraulics <- function(temp_roots, hydraulics){
  
  # Merge output of MARSHAL on specific root segment
  temp_roots$suf <- as.vector(hydraulics$suf)
  temp_roots$suf_eq <- as.vector(hydraulics$suf_eq)
  temp_roots$suf1 <- as.vector(hydraulics$suf1)
  temp_roots$kx <- as.vector(hydraulics$kx)
  temp_roots$kr <- as.vector(hydraulics$kr)
  temp_roots$jr <- as.vector(hydraulics$jr)
  temp_roots$jr_eq <- as.vector(hydraulics$jr_eq)
  temp_roots$jxl_eq <- as.vector(hydraulics$jxl_eq)
  temp_roots$psi_eq <- as.vector(hydraulics$psi_eq)
  temp_roots$psi <- as.vector(hydraulics$psi)
  temp_roots$jxl <- as.vector(hydraulics$jxl)
  temp_roots$psi_soil <- as.vector(hydraulics$psi_soil)
  return(temp_roots)
}



soil_map <- function(soil, Maxdepth = 100){
  
  soil_global <- soil
  VecAge <- c(0.8, unique(soil_global$value))
  water = NULL
  k = 1
  for(i in unique(soil_global$value)){
    table_soil <- soil_global[soil_global$value==i & soil_global$z >= -Maxdepth,]
    for(j in 1:nrow(table_soil)){
      # j=5
      tmp <- tibble(x = c(i,soil_global$value[soil_global$value>i][1], soil_global$value[soil_global$value>i][1], i),
                    y = c(table_soil$z[j], table_soil$z[j], table_soil$z[j+1], table_soil$z[j+1]), 
                    Psi = soil_global$psi[soil_global$value == i][j],
                    SSF = soil_global$SSF[soil_global$value == i][j],
                    Sink = soil_global$Sink[soil_global$value == i][j],
                    Flux = soil_global$Flux[soil_global$value == i][j],
                    moisture = soil_global$moisture[soil_global$value == i][j],
                    grid = k)
      
      k = k +1
      
      water <- rbind(water, tmp)
    }
    
  }
  
  return(water)
}



run_HydrusCouvreurMARSHAL<- function(all_roots, 
                          conductivities = "./www/conductivities.xlsx",
                          soil_type = "loam", 
                          rhizosphere = F){
  
  
  all_roots <- suppressWarnings(suppressMessages(root_transform(all_roots)))
  
  conductivities <- read_excel(conductivities)
  
  soil <- soil_initial(soil_type = soil_type, field_capacity = 1)
  
  soil_param <- read_excel("./www/Soil_type.xlsx")%>%
    filter(type == soil_type)%>%
    mutate(Ksat = Ks, lambda = l)
  
  tpots = -15000 # hPa
  SOIL = RWU <- NULL
  dt = 0.05 # 0.05 , 
  time_sequence <- seq(1, 30, dt)
  for(t in time_sequence){
    if(t %in% seq(2,30,2)){
      print(t)
    }
      
    temp_roots <- all_roots%>%
      filter(time <= t)
    
    # -----------------------------
    # Run MARSHAL
    # -----------------------------
    hydraulics <- getSUF(temp_roots, 
                         conductivities, 
                         soil, 
                         hetero = T, 
                         Psi_collar = tpots, soil_param)
    
    
    temp_roots <- add_hydraulics(temp_roots, hydraulics)
    if(rhizosphere){
      Tact <- hydraulics$tact_eq
      Tpot <- hydraulics$tpot_eq
      krs <- hydraulics$ksrs
      RLDWU <- temp_roots%>% # gather information by layer
        mutate(rz2 = round((z1+z2)/2))%>%
        dplyr::group_by(rz2)%>%
        dplyr::summarise(suf = sum(suf_eq),
                         ps = sum(psi),
                         jr = sum(jr_eq),
                         jx = sum(jxl_eq),
                         su = sum(suf))%>%
        ungroup()
      
    }else{
      Tact <- hydraulics$tact # cm3 d -1
      Tpot <- hydraulics$tpot
      krs <- hydraulics$krs # cm4 hPa-1 d-1
      RLDWU <- temp_roots%>% # gather information by layer
        mutate(rz2 = round((z1+z2)/2))%>%
        dplyr::group_by(rz2)%>%
        dplyr::summarise(suf = sum(suf1),
                         ps = sum(psi),
                         jr = sum(jr),
                         jx = sum(jxl),
                         su = sum(suf))%>%
        ungroup()
    }

    Beta <- rep(0, 101)
    
    RWU_tmp <- temp_roots%>% # gather information by layer
      mutate(rz2 = round((z1+z2)/2))%>%
      dplyr::group_by(rz2)%>%
      dplyr::summarise(suf = sum(suf1),
                       suf_eq = sum(suf_eq),
                       ps = sum(psi),
                       jr = sum(jr),
                       jx = sum(jxl),
                       su = sum(suf),
                       jr_eq = sum(jr_eq),
                       jx_eq = sum(jxl_eq))%>%
      ungroup()
    RWU <- rbind(RWU, RWU_tmp)
    
    Q_dou = sum(temp_roots$jr)
    #message(paste0("the root system is suppose to take up: ",Q_dou,"[cm3 d-1]"))
    
    #include the Standard upatke fraction
    Beta[which(soil$z %in% RLDWU$rz2 )] <- rev(RLDWU$suf) # from above to below
    SSF <- data.frame(suf = Beta, h = soil$psi)
    SUF <- rev(RLDWU$suf[RLDWU$rz2 >= -100 & RLDWU$rz2 <= 0])
    Jr <- rev(RLDWU$jr[RLDWU$rz2 >= -100 & RLDWU$rz2 <= 0])
    # overwirte the profile boundary condition
    write.profile.dat(project.path = "./CouvreurV2", SSF)
    # print("profile.dat correctly overwrited")
    # Calculate Kcomp
    Hsr <- soil$psi[which(soil$z %in% RLDWU$rz2 )]
    kcomp <- krs[1]
    if(length(unique(Hsr))> 1){
      Hseq <- t(Hsr) %*% t(t(SUF))
      kcomp = (Jr -  Q_dou*SUF) %*% ((Hsr-Hseq[1])*SUF)^(-1)
    }else{
      Hseq <- unique(Hsr)
    }
    kcomp <- kcomp[1]
    write.options.in(project.path = "./CouvreurV2", krs/75/15, kcomp/75/15)
    
    
    # Good idea to check if Tpot works similarly
    atm_bc_data <- data.frame(tAtm = 1, Prec = 0, rSoil = 0, 
                              rRoot = Tact/75/15, hCritA = 15000, rB = 0, hB = 0, ht = 0, 
                              RootDepth = 0)
    # overwrite the atmposheric boundary condition of hydrus.
    write.atmosph.in("./CouvreurV2/",
                     maxAL = 1,
                     deltaT = 1,
                     atm_bc_data,
                     hCritS = 15000,
                     input.pet = F)
    
    system("./H1D_calc.exe", show.output.on.console = F)
    
    hydrus <- read.nod_inf(project.path = "./CouvreurV2", 
                           out.file = paste0("Nod_Inf.out"))  
    soil <- hydrus%>%
      filter(Time == dt)%>%
      transmute(id = Node,
                z = Depth,
                value = t,
                psi = Head,
                moisture = Moisture,
                SSF = SSF$suf,
                Sink = Sink, Flux = Flux)
    SOIL <- rbind(SOIL, soil)
    
    Tact_hydrus <- sum(hydrus$Sink)*75*15
    
  }
  return(data <- list(SOIL, RWU))
}

message("...done")


plot_soil <- function(SOIL, type = "SSF"){
  
  water <- soil_map(SOIL, Maxdepth = 50)
  
  if(type == "SSF"){
    pl <- water %>%
      ggplot()+
      geom_polygon(aes(x, y, group = grid, fill = SSF))+
      scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = 0.25)+
      theme_classic()
  }
  if(type == "Sink"){
    pl <- water %>%
      ggplot()+
      geom_polygon(aes(x, y, group = grid, fill = Sink*75*15))+
      scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = 0.3)+
      theme_classic()
  }
  if(type == "Flux"){
    pl <- water %>%
      ggplot()+
      geom_polygon(aes(x, y, group = grid, fill = Flux*75*15))+
      scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = 0)+
      theme_classic()
  }
  if(type == "Psi"){
    pl <- water %>%
      ggplot()+
      geom_polygon(aes(x, y, group = grid, fill = Psi))+
      scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = -1800)+
      theme_classic()
  }
  
  print(pl)
}

read.tlevel.out <- function (project.path, out.file = "T_Level.out", output = NULL, 
          warn = FALSE, ...) 
{
  if (is.null(output) | missing(output)) {
    output = output = c("rTop", "rRoot", "vTop", "vRoot", 
                        "vBot", "sum(rTop)", "sum(rRoot)", "sum(vTop)", 
                        "sum(vRoot)", "sum(vBot)", "hTop", "hRoot", "hBot", 
                        "RunOff", "sum(Runoff)", "Volume", "sum(Infil)", 
                        "sum(Evap)", "TLevel", "Cum(WTrans)", "SnowLayer")
  }
  options(warn = -1)
  tlevel_out = data.table::fread(input = file.path(project.path, 
                                                   out.file), fill = TRUE, blank.lines.skip = T, skip = 6, 
                                 header = T)
  tlevel_out = apply(tlevel_out, MARGIN = 2, FUN = as.numeric)
  tlevel_out = na.omit(tlevel_out)
  tlevel_out = data.frame(tlevel_out, check.names = FALSE, 
                          row.names = NULL)
  tstart_ind = which(tlevel_out$TLevel == 1)
  sum_cols_ind = grep("sum", names(tlevel_out))
  sum_col_names = names(tlevel_out)[sum_cols_ind]
  if (length(tstart_ind) == 1) {
    tlevel_out = tlevel_out
  }
  else {
    for (i in 2:length(tstart_ind)) {
      run1_totals = tlevel_out[(tstart_ind[i] - 1), sum_cols_ind]
      if (i == length(tstart_ind)) {
        run_i_ind = tstart_ind[i]:nrow(tlevel_out)
      }
      else {
        run_i_ind = tstart_ind[i]:(tstart_ind[i + 1] - 
                                     1)
      }
      tout_j = tlevel_out[run_i_ind, ]
      for (j in sum_col_names) {
        tlevel_out[run_i_ind, j] = tout_j[, j] + run1_totals[[j]]
      }
    }
  }
  options(warn = 0)
  return(tlevel_out)
}

