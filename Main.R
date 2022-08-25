


library(tidyverse)
setwd("C:/Users/heymansad/Documents/GitHub/HydrusCouMa/")
source("./Rscript/io_function.R") # function for this script to works
source("./Rscript/getSUF.R") # MARSHAL


data <- data.table::fread("C:/Users/heymansad/Documents/GitHub/HydrusCouMa/www/ET0badlauchstaedt_2019", sep = " ")
rain <- data.table::fread("C:/Users/heymansad/Documents/GitHub/HydrusCouMa/www/rainbadlauchstaedt_2019", sep = " ")




data %>% ggplot(aes(V1,V2))+geom_line()+xlim(65,75)

all_root <- data.table::fread("C:/Users/heymansad/Documents/GitHub/HydrusCouMa/www/63_rootsystem.txt")%>%
  mutate(length = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2))

all_root = root_transform(all_root)

plot_rs(all_root)

conductivities <- read_excel("./www/conductivities.xlsx")

soil <- soil_initial(soil_type = "loam", field_capacity = 0.1)

soil_param <- read_excel("./www/Soil_type.xlsx")%>%
  filter(type == "loam")%>%
  mutate(Ksat = Ks, lambda = l)

tpots = -15000
OUT = SOIL <- NULL
dt = 1
time_sequence <- seq(20, 60, dt)
for(t in time_sequence){
  tt = t-min(time_sequence)
  print(t)
  temp_roots <- all_root%>%
    filter(time <= t)
  
  # -----------------------------
  # Run MARSHAL
  # -----------------------------
  hydraulics <- getSUF(temp_roots, 
                       conductivities, 
                       soil, 
                       hetero =F, 
                       Psi_collar = tpots, soil_param)
  
  
  temp_roots <- add_hydraulics(temp_roots, hydraulics)
  krs <- hydraulics$krs # cm4 hPa-1 d-1
  kcomp <- krs
  
  #############################
  # # Calculate Kcomp
  # Hsr <- soil$psi[which(soil$z %in% temp_roots$rz2 )]
  # kcomp <- krs[1]
  # if(length(unique(Hsr))> 1){
  #   Hseq <- t(Hsr) %*% t(t(SUF))
  #   kcomp = (Jr -  Q_dou*SUF) %*% ((Hsr-Hseq[1])*SUF)^(-1)
  # }else{
  #   Hseq <- unique(Hsr)
  # }
  # kcomp <- kcomp
  ##############################
  
  Beta <- rep(0, 101)
    RLDWU <- temp_roots%>% # gather information by layer
    mutate(rz2 = round((z1/2+z2/2)/2)*2)%>%
    dplyr::group_by(rz2)%>%
    dplyr::summarise(suf = sum(suf1),
                     ps = sum(psi),
                     jr = sum(jr),
                     jx = sum(jxl),
                     su = sum(suf),
                     jr_eq = sum(jr),
                     jx_eq = sum(jxl))%>%
    ungroup()
  Beta[which(soil$z %in% RLDWU$rz2 )] <- rev(RLDWU$suf) # from above to below
  SSF <- data.frame(suf = Beta, h = soil$psi)

  # overwirte the profile boundary condition
  write.profile.dat(project.path = "./Day5", SSF)
  # message("profile.dat is correctly written")


  write.options.in(project.path = "./Day5", krs/75/15, kcomp/75/15)
  # message("options.in is correctly written")
  
  
  atm_bc_data <- data.frame(tAtm = round(seq(1/24,1,1/24),4), Prec = rain$V1[round(tt)+65]/10, rSoil = 0, 
                            rRoot = data$V2[(2+tt*24+65*24):(1+(tt+1)*24+65*24)], hCritA = 15000, rB = 0, hB = 0, ht = 0, 
                            RootDepth = 0)
  # overwrite the atmposheric boundary condition of hydrus.
  write.atmosph.in("./Day5/",
                   maxAL = 24,
                   deltaT = 1,
                   atm_bc_data,
                   hCritS = 15000,
                   input.pet = F)
  # message("atmosph.in is correctly written")
  
  system("./H1D_calc.exe", show.output.on.console = F)
  
  hydrus <- read.nod_inf(project.path = "./Day5", 
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
  SOIL <- rbind(SOIL, hydrus%>% 
                  mutate(SSF = rep(SSF$suf,25), krs = krs,id = Node, z = Depth, value = Time, psi = Head, moisture = Moisture)%>%
                  filter(Time != 0)%>%
                  mutate(Time = sort(rep(data$V1[(2+tt*24+65*24):(1+(tt+1)*24+65*24)],101))))
  
  out_data = read.tlevel.out(project.path = "./Day5", out.file = paste0("T_Level.OUT"))%>%
    mutate(Time = Time + tt+65)
  OUT <- rbind(OUT, out_data)
  
  # Tact_hydrus <- sum(hydrus$Sink)*75*15
}



SOIL%>%
  ggplot(aes(z, psi))+
  geom_line(aes(colour = Time, group = factor(Time)), alpha = 0.1)+
  coord_flip()+
  viridis::scale_colour_viridis()+
  xlim(-60,0)

SOIL%>%
  ggplot(aes(z, SSF))+
  geom_line(aes(colour = Time, group = factor(Time)), alpha = 0.1)+
  coord_flip()+
  viridis::scale_colour_viridis()+
  xlim(-60,0)

rain$day = 1:nrow(rain)

SoilFlux<- SOIL %>%
  dplyr::group_by(Time)%>%
  dplyr::summarise(Soil_Water_balance = -sum(Flux))%>%
  ungroup()

ggplot()+
  geom_line(aes(Time, Soil_Water_balance), data = SoilFlux)+
  geom_segment(aes(x = day, xend = day, y = 0,yend = V1/10), size = 3, alpha = 0.2, data = rain, colour = "blue")+
  xlim(65,105)

OUT %>%
  ggplot()+
  geom_line(aes(Time, rRoot))+
  geom_line(aes(Time, vRoot), colour = 'red')





