

library(readxl)
setwd("~/HydraMaize/22-04 Hydrus_Cou/")
source("./Rscript/io_function.R") # function for this script to works
source("./Rscript/getSUF.R") # MARSHAL
library(data.table)
library(Matrix)
library(tidyverse)


all_roots <- fread(paste0("./www/35_rootsystem.txt"), header = T)
all_roots <- fread(paste0("./www/30_rhizo.txt"), header = T)
all_roots <- fread(paste0("./www/30_core.txt"), header = T)

all_roots%>%
  ggplot()+
  geom_segment(aes(x = x1, xend = x2, y = z1, yend = z2))+coord_fixed()+theme_classic()

all_roots <- root_transform(all_roots)

conductivities <- read_excel("./www/conductivities.xlsx")

soil <- soil_initial(soil_type = "loam", field_capacity = 0.1)

soil_param <- read_excel("./www/Soil_type.xlsx")%>%
  filter(type == "loam")%>%
  mutate(Ksat = Ks, lambda = l)

tpots = -15000 # hPa
SOIL <- NULL
dt = 0.1 # 0.05 , 
time_sequence <- seq(1, 30, dt)
for(t in time_sequence){
  print(t)
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
  Tact <- hydraulics$tact # cm3 d -1 
  Tpot <- hydraulics$tpot
  krs <- hydraulics$krs # cm4 hPa-1 d-1
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
  Q_dou = sum(temp_roots$jr)
  message(paste0("the root system is suppose to take up: ",Q_dou,"[cm3 d-1]"))
  
    #include the Standard upatke fraction
  Beta[which(soil$z %in% RLDWU$rz2 )] <- rev(RLDWU$suf) # from above to below
  SSF <- data.frame(suf = Beta, h = soil$psi)
  SUF <- rev(RLDWU$suf[RLDWU$rz2 >= -200 & RLDWU$rz2 <= 0])
  Jr <- rev(RLDWU$jr[RLDWU$rz2 >= -200 & RLDWU$rz2 <= 0])
   # overwirte the profile boundary condition
  write.profile.dat(project.path = "./CouvreurV2", SSF)
  # print("profile.dat correctly overwrited")
  # Calculate Kcomp
  Hsr <- soil$psi[which(soil$z %in% RLDWU$rz2 )]
  kcomp <- krs[1]
  if(length(unique(Hsr))> 1){
    Hseq <- t(Hsr) %*% t(t(SUF))
    #kcomp <- t(Hsr-Hseq[1]) %*% (Q_dou/SUF - Tact) / (t(Hsr-Hseq[1]) %*% (Hsr-Hseq[1]))[1]
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


water <- soil_map(SOIL, Maxdepth = 50)

water %>%
  ggplot()+
  geom_polygon(aes(x, y, group = grid, fill = SSF))+
  scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = 0.25)+
  theme_classic()

water %>%
  ggplot()+
  geom_polygon(aes(x, y, group = grid, fill = Sink*75*15))+
  scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = 0.3)+
  theme_classic()

water %>%
  ggplot()+
  geom_polygon(aes(x, y, group = grid, fill = Flux*75*15))+
  scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = 0)+
  theme_classic()

water %>%
  ggplot()+
  geom_polygon(aes(x, y, group = grid, fill = Psi))+
  scale_fill_gradient2(low = "firebrick1", mid = "white", high = "dodgerblue2", midpoint = -1800)+
  theme_classic()

SOIL %>%
  dplyr::group_by(value)%>%
  dplyr::summarise(Sink = sum(Sink))%>%
  ungroup()%>%
  ggplot()+
  geom_line(aes(x = value, y = Sink))

