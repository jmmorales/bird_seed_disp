#######################################################
###ARANGE DATA########################
require("mvtnorm")
require("tmvtnorm")
library("bbmle")
library("vegan")
library("emdbook")
library("lattice")
require("VGAM")
require("R.utils")

load("Matrices_Llao.RData")
ncol(bird_ab)
#Comprobamos que los tags y el nombre completo
#tanto de los pájaros como de las plantas siguen el mismo orden

View(cbind(as.character(plant_ids), as.character(plants_name)))
View(cbind(as.character(bsp), as.character(frug_name)))
#bird_ab = bird_ab[,-which(colnames(bird_ab) == "DIXPIPR")]
##La elimino porque no había filogenia

rowSums(bird_ab)#

bsp_community = sort(as.character(phy_Llao$X))
spf_community = seq(from =1, to = length(bsp_community), by = 1)
MT1 = read.csv("traits.csv")
View(cbind(bsp_community, MT1))

#reordenamos MT según bsp_community
MT = NULL
for(i in 1:length(bsp_community))
{
  tmp = MT1[which(MT1$frug_name == bsp_community[i]),]
  MT = rbind(MT, tmp)
}

colnames(MT)
TT = cbind(rep(1, nrow(MT)), MT)
TT = TT[,-2]
View(TT)
colnames(TT)[1] = "Intercept"
TT = as.matrix(TT)

###Rorganize phylogeny (following alphabetic order)
CC1 <- as.matrix(phy_Llao)
CC1 = CC1[,-1]
CC1 = apply(CC1,2,as.numeric)

CC =  matrix(0,nrow=length(bsp_community), ncol=length(bsp_community)) 
for(i in 1:length(bsp_community))
{
  tmp = which(colnames(CC1) == bsp_community[i])
  for(j in 1:length(bsp_community))
  {
    tmp2 = which(colnames(CC1) == bsp_community[j])
    CC[i, j] = CC1[tmp, tmp2]
  }
  
}

#View(CC)
colnames(CC) = bsp_community

#Estimate GPT following Herrera 1984
write.csv(bird_data, "Llao_traits_raw.csv", row.names = F)#convierto el tamaño en gramos
traits_raw = read.csv("Llao_traits_raw.csv")
mean_GPT=numeric(nrow(TT))
size_gr = traits_raw$Size*1000
scale_GPT = numeric(nrow(TT))
shape = 1.59 #ejemplo de Turdus de Morales 2013

for(i in 1:nrow(TT))
{
  tmp = TT[i,]
  if(TT[i,4] == 1)
  {tmp2 =  26.173+0.481*size_gr[i] 
  mean_GPT[i] = tmp2
  scale_GPT[i] = tmp2/1.59
  }else{
    tmp2 = 62.507+0.957*size_gr[i] 
    mean_GPT[i] = tmp2
    scale_GPT[i] = tmp2/1.59 
  }
}


#i =2
#mean(rgamma(1000, scale = scale_GPT[i], shape = 1.59))
#mean_GPT[i]


#Define the scale of landscape matrices (use the same scale as MCMC)

Fc = Fc1
nF = nF1 #en la simulación le aplica el logaritmo más 1
Distance = Dist / 10 #cuando se ajustó estaba expresado en dm
distos = Distance##matriz distancias entre celdas precalculadas
distos= rbind(distos, numeric(ncol(Distance))*NA)
distos = cbind(distos, numeric(nrow(distos))*NA)

map_cell = read.csv("map_Ncell.csv")
Cells1 = NULL
for(i in 1:nrow(map_cell))
{
  for(j in 1:ncol(map_cell))
  {
    out = c(map_cell[i,j], i, j)
    Cells1 = rbind(Cells1,out)
  }
}
View(Cells1)
Cells1 = data.frame(Cells1)
colnames(Cells1) = c("Ncell", "y_coord", "x_coord")
xp = Cells1$x_coord
yp = Cells1$y_coord

#plot_v
plot_v = c()
for(i in 1:nrow(map))
{
  plot_v = c(plot_v, map[i,])
}

Ab_move = Ab_move1##abundancia relativa de las especies de las celdas en el paisaje
#Fl = Fl44#abundancia de las especies en el paisaje
Fc = Fc1#abundancia relativa de las especies en las celdas
B = B1#distancia al borde de las celdas
rows = max(Cells1$y_coord)
cols = max(Cells1$x_coord)
ncells = rows*cols
Forb_int = Fi[,-1]
View(Forb_int)
Forb_int = Forb_int[,-ncol(Forb_int)]
ncol(Forb_int) == length(plant_ids)

cover_v = c()
for(i in 1:nrow(cover))
{
  cover_v = c(cover_v, cover[i,])
}
cover = cover_v

open = cover
open[which(cover==0)]=1
open[which(cover==1)]=0
View(cbind(cover,open))


###pids empezando en aquellas celdas donde se ha observado la specie
pids = seq(from = 1, to =ncells)
seqs = seq[which(!is.na(seq$Ncell)),]

#Celdas de las que va a samplear los movimientos iniciales

initC = vector('list', length(bsp_community))

for(i in 1:length(bsp_community))
{
  tmp = seqs[which(seqs$frug_name == bsp_community[i]),]
  tmp2 = unique(tmp$Ncell)
  tmp2 = tmp2[which(tmp2 != (ncells+1))]
  initC [[i]] = tmp2
}

# Load MCMC fit and organize ----------------------------------------------


load("MCMCfit_Llao_model_input.RData")
ns = length(bsp_community)
activities = c("move", "choose", "consume", "perch")
params = vector("list", 4)##una lista de THETAS para cada actividad
params[[1]] = array(NA, c(ns, 7)) ##ad, b_c, afs, bf, bs, a0, b0
params[[2]] = array(NA, c(ns, 2))#a_sc, b_sc
params[[3]] = array(NA, c(ns, 2))#lambda, p0
params[[4]] = array(NA, c(ns, 2))#shape, rate
nch = 3 ##número de cadenas que corrieron en el MCMC
nreps = 100
iters = vector("list", 4)
store_length= c(ncol(THE_M[1,,]), ncol(THE_CH[1,,]), ncol(THE_CNS[1,,]), ncol(THE_PER[1,,]))
for(i in 1:4)#esto es simplemente porque para choice hubo #mientras que hubo 10000 en los otros20000 iteraciones
{
 iters[[i]]  = sample(seq((store_length[i]-1)), size =nreps, replace = T)
 }
#(store_length[i]-1) hago -1 para que en cso de fallo (de rmvnorm) como suma 1 nunca se quede fuera de rango
##LIKELIHOOD FUNCTIONS

source("like_mov_May18.R")

save.image("Llao_data_input_sim.RData")


####SIMULATION
rm(list =ls())
require("mvtnorm")
require("tmvtnorm")
library("bbmle")
library("vegan")
library("emdbook")
library("lattice")
require("VGAM")
require("R.utils")

load("Llao_data_input_sim.RData")

Ntracks = 1000
ntracks = ceiling(Ntracks/ nrow(Fl))



#"n" listas (n= nreps) cada lista guarda información sobre los movimientos (ntracks*mov.por tracks)
#Consumption and dispersal outputs
KK <- vector('list', nreps)#Guarda las celdas en las que se defecan semillas
ZZ <- vector('list', nreps)#Guarda las celdas entre las que se producen los movimientos
IDC <- vector('list', nreps)#Guarda la identidad de las semillas que ha comido
NC <- vector("list", nreps)#Guarda el número de frutos que ha comido
BB <- vector('list', nreps)#Guarda la especie de pájaro
WEEK <- vector('list', nreps)#Guarda la info de las semanas
Nos <- c() # will hold number of defecations out of plot length == rep*weeks
Tr = vector('list', nreps)#Guarda la info de los tracks
#Movement outputs
ZZ_mov <- vector('list', nreps) #Guarda las celdas de los movimientos (indep.) de que haya consumo o no
RR_mov <- vector('list', nreps)#Distancias de vuelo (indep del consumo)
BB_mov <- vector('list', nreps)#Guarda los identificadores de los pájaros para todos los movimientos
WEEK_mov <- vector('list', nreps)#Guarda la semana que está siendo simulada
Tr_mov <-  vector('list', nreps)
#Changes in diversity, equitativity...
Div<-  vector('list', nreps)
Ri<-  vector('list', nreps)
J<-  vector('list', nreps)
#Summaries of repetition
Mov <- vector('list', nreps)
Inter <- vector('list', nreps)
Seedrain <- vector('list', nreps)


  for (j in 1:nreps){
    print(paste("rep", " ", j, sep = ""))
    print(paste(j/nreps, " ", "%", sep = ""))
    ##Create vectors to hold values of week
    K_w = c()
    Zm_w = c()
    idf_w = c()
    consumed_w = c()
    B_w = c()
    No_w = c()
    Week_w = c()
    Track_w = c()
    #los que tienen un s guarda las celdas independientemente de que haya 
    #consumo de frutos o no
    Zs_mov_w = c()
    BB_mov_w = c()
    R_mov_w = c()
    Track_mov_w = c()
    Week_mov_w = c()
    #changes in diversity
    deltaH_w = c()
    deltaJ_w = c()
    deltaS_w = c()
    Rain = matrix(0, ncells,length(plant_ids))#va a ir guardando el seed rain
   #Avoid depletion among replicates
     nFj = nF #evitar que en la siguiente repetición haya un depletion de los frutos
    #que se han consumido en la anterior
    #:nrow(Fl)
####################################################
#SAMPLE SPECIES SPECIFIC RESPONSES FROM POSTERIORS
# With direct and indirect (derived from traits) information

    ###Movement##########################
    tmp = which(colnames(TT) %in% move_traits)
    TT_M = TT[, tmp]
    M_M = TT_M %*% ZT_M[,,iters[[1]][j]]
    nt = nrow(ZT_M[,,iters[[1]][j]])
    np = ncol(ZT_M[,,iters[[1]][j]])
    M_Mv = as.vector(M_M)#vectorizo agrupando por parametros
    #parametro 1 (todas las spp), param2 (todas las spp) etc.
    rho = RH_M[iters[[1]][j]]
    C2_M = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
    V_M = kronecker(SI_M[,,iters[[1]][j]], C2_M)
    
    res <- try(withTimeout({
      THE_M_i_v = rmvnorm(n =1, mean = M_Mv, sigma = V_M);
    }, timeout=60, onTimeout = "error"))
    
    if(inherits(res, "try-error"))
    {
      #error handling code, go to another parameterset (sum 1 to the row sampled)
      tmp = which(colnames(TT) %in% move_traits)
      TT_M = TT[, tmp]
      M_M = TT_M %*% ZT_M[,,(iters[[1]][j]+1)]
      nt = nrow(ZT_M[,,(iters[[1]][j]+1)])
      np = ncol(ZT_M[,,(iters[[1]][j]+1)])
      M_Mv = as.vector(M_M)#vectorizo agrupando por parametros
      #parametro 1 (todas las spp), param2 (todas las spp) etc.
      rho = RH_M[(iters[[1]][j]+1)]
      C2_M = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
      V_M = kronecker(SI_M[,,(iters[[1]][j]+1)], C2_M)
      THE_M_i_v = rmvnorm(n =1, mean = M_Mv, sigma = V_M)
      THE_M_i = matrix(THE_M_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% move_spp)
        {
          tmp = which(move_spp == bsp_community[s])
          params[[1]][s,] = THE_M[tmp,, (iters[[1]][j]+1)]
        }else{
          params[[1]][s,] = THE_M_i[s,]
        }
      }
      
     }else{
      THE_M_i = matrix(THE_M_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% move_spp)
        {
          tmp = which(move_spp == bsp_community[s])
          params[[1]][s,] = THE_M[tmp,, iters[[1]][j]]
        }else{
          params[[1]][s,] = THE_M_i[s,]
        }
      }
      }
    
    ###Fruit choice###############################
    tmp = which(colnames(TT) %in% choice_traits)
    TT_CH = TT[, tmp]
    M_CH = TT_CH %*% ZT_CH[,,iters[[2]][j]]
    nt = nrow(ZT_CH[,,iters[[2]][j]])
    np = ncol(ZT_CH[,,iters[[2]][j]])
    M_CHv = as.vector(M_CH)#vectorizo agrupando por parametros
    #parametro 1 (todas las spp), param2 (todas las spp) etc.
    rho = RH_CH[iters[[2]][j]]
    C2_CH = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
    V_CH = kronecker(SI_CH[,,iters[[2]][j]], C2_CH)
    
    res <- try(withTimeout({
      THE_CH_i_v = rmvnorm(n =1, mean = M_CHv, sigma = V_CH);
    }, timeout=60, onTimeout = "error"))
    
    if(inherits(res, "try-error"))
    {
      #error handling code, go to another parameterset (sum 1 to the row sampled)
      tmp = which(colnames(TT) %in% choice_traits)
      TT_CH = TT[, tmp]
      M_CH = TT_CH %*% ZT_CH[,,(iters[[2]][j]+1)]
      nt = nrow(ZT_CH[,,(iters[[2]][j]+1)])
      np = ncol(ZT_CH[,,(iters[[2]][j]+1)])
      M_CHv = as.vector(M_CH)#vectorizo agrupando por parametros
      #parametro 1 (todas las spp), param2 (todas las spp) etc.
      rho = RH_CH[(iters[[2]][j]+1)]
      C2_CH = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
      V_CH = kronecker(SI_CH[,,(iters[[2]][j]+1)], C2_CH)
      THE_CH_i_v = rmvnorm(n =1, mean = M_CHv, sigma = V_CH)
      THE_CH_i = matrix(THE_CH_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% choice_spp)
        {
          tmp = which(choice_spp == bsp_community[s])
          params[[2]][s,] = THE_CH[tmp,, (iters[[2]][j]+1)]
        }else{
          params[[2]][s,] = THE_CH_i[s,]
        }
      }
      
    }else{
      THE_CH_i = matrix(THE_CH_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% choice_spp)
        {
          tmp = which(choice_spp == bsp_community[s])
          params[[2]][s,] = THE_CH[tmp,, iters[[2]][j]]
        }else{
          params[[2]][s,] = THE_CH_i[s,]
        }
      }
   }
    
    
    ###Fruit consumption####################################
    tmp = which(colnames(TT) %in% cons_traits)
    TT_CNS = TT[, tmp]
    M_CNS = TT_CNS %*% ZT_CNS[,,iters[[3]][j]]
    nt = nrow(ZT_CNS[,,iters[[3]][j]])
    np = ncol(ZT_CNS[,,iters[[3]][j]])
    M_CNSv = as.vector(M_CNS)#vectorizo agrupando por parametros (columnas)
    #parametro 1 (todas las spp), param2 (todas las spp) etc.
    rho = RH_CNS[iters[[3]][j]]
    C2_CNS = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
    V_CNS = kronecker(SI_CNS[,,iters[[3]][j]], C2_CNS)
    
    res <- try(withTimeout({
      THE_CNS_i_v = rtmvt(n =1, mean = M_CNSv, sigma = V_CNS,
                          lower = rep(0, ns*2),
                          upper = c(rep(Inf, ns), rep(1, ns)));
    }, timeout=60, onTimeout = "error"))
    
    if(inherits(res, "try-error"))
    {
      #error handling code, go to another parameterset (sum 1 to the row sampled)
      tmp = which(colnames(TT) %in% cons_traits)
      TT_CNS = TT[, tmp]
      M_CNS = TT_CNS %*% ZT_CNS[,,(iters[[3]][j]+1)]
      nt = nrow(ZT_CNS[,,(iters[[3]][j]+1)])
      np = ncol(ZT_CNS[,, (iters[[3]][j]+1)])
      M_CNSv = as.vector(M_CNS)#vectorizo agrupando por parametros (columnas)
      #parametro 1 (todas las spp), param2 (todas las spp) etc.
      rho = RH_CNS[(iters[[3]][j]+1)]
      C2_CNS = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
      V_CNS = kronecker(SI_CNS[,,(iters[[3]][j]+1)], C2_CNS)
      THE_CNS_i_v = rtmvt(n =1, mean = M_CNSv, sigma = V_CNS,
                          lower = rep(0, ns*2),
                          upper = c(rep(Inf, ns), rep(1, ns)))
      THE_CNS_i = matrix(THE_CNS_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% cons_spp)
        {
          tmp = which(cons_spp == bsp_community[s])
          params[[3]][s,] = THE_CNS[tmp,, (iters[[3]][j]+1)]
        }else{
          params[[3]][s,] = THE_CNS_i[s,]
        }
      }
      
    }else{
      THE_CNS_i = matrix(THE_CNS_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% cons_spp)
        {
          tmp = which(cons_spp == bsp_community[s])
          params[[3]][s,] = THE_CNS[tmp,, iters[[3]][j]]
        }else{
          params[[3]][s,] = THE_CNS_i[s,]
        }
      }
    }
    
    
    ###Perching times#################
    tmp = which(colnames(TT) %in% perch_traits)
    TT_PER= TT[, tmp]
    M_PER = TT_PER %*% ZT_PER[,,iters[[4]][j]]
    nt = nrow(ZT_PER[,,iters[[4]][j]])
    np = ncol(ZT_PER[,,iters[[4]][j]])
    M_PERv = as.vector(M_PER)#vectorizo agrupando por parametros (columnas)
    #parametro 1 (todas las spp), param2 (todas las spp) etc.
    rho = RH_PER[iters[[4]][j]]
    C2_PER = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
    V_PER = kronecker(SI_PER[,,iters[[4]][j]], C2_PER)
    res <- try(withTimeout({
      THE_PER_i_v = rtmvt(n =1, mean = M_PERv, sigma = V_PER,
                          lower = rep(0, ns*2),
                          upper = rep(Inf, ns*2));
    }, timeout=60, onTimeout = "error"))
    
    if(inherits(res, "try-error"))
    {
      tmp = which(colnames(TT) %in% perch_traits)
      TT_PER= TT[, tmp]
      M_PER = TT_PER %*% ZT_PER[,,(iters[[4]][j]+1)]
      nt = nrow(ZT_PER[,,(iters[[4]][j]+1)])
      np = ncol(ZT_PER[,,(iters[[4]][j]+1)])
      M_PERv = as.vector(M_PER)#vectorizo agrupando por parametros (columnas)
      #parametro 1 (todas las spp), param2 (todas las spp) etc.
      rho = RH_PER[(iters[[4]][j]+1)]
      C2_PER = (rho*CC)+ ((1-rho)*diag(nrow(CC)))
      V_PER = kronecker(SI_PER[,,(iters[[4]][j]+1)], C2_PER)
      THE_PER_i_v = rtmvt(n =1, mean = M_PERv, sigma = V_PER,
                          lower = rep(0, ns*2),
                          upper = rep(Inf, ns*2))
      THE_PER_i = matrix(THE_PER_i_v, nrow = ns, ncol = np)
      #Create theta with direct and indirect data
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% perch_spp)
        {
          tmp = which(perch_spp == bsp_community[s])
          params[[4]][s,] = THE_PER[tmp,, (iters[[4]][j]+1)]
        }else{
          params[[4]][s,] = THE_PER_i[s,]
        }
      }
    }else{
      THE_PER_i = matrix(THE_PER_i_v, nrow = ns, ncol = np)
      for(s in 1:ns)
      {
        if(as.character(bsp_community[s]) %in% perch_spp)
        {
          tmp = which(perch_spp == bsp_community[s])
          params[[4]][s,] = THE_PER[tmp,, iters[[4]][j]]
        }else{
          params[[4]][s,] = THE_PER_i[s,]
        }
      } 
    }
############################################################################### 
    # Begin simulation  -------------------------------------------------------

        for (week in 1:nrow(Fl))##Week loop
    {
    bab = bird_ab[week,]
    spp <- sample(spf_community,size=ntracks, prob=bab, replace=T)  # bird species sampled from relative abundance values   ###DUDA CON ESE 5
    Bird =NA
    Zs = NA
    Track = NA
    # Bird movement -----------------------------------------------------------
    for (i in 1:ntracks){
      nid = spp[i]#species identification
      iC_spp = initC [[nid]]
      if(length(iC_spp) > 0)##Some observations began out of plots 
      {
        z = sample(iC_spp,1) #choose the cell where track begins
      }else{
        z = sample(pids,1) #choose randomly a cell in which the species was observed
      }
      
      nz = z
      count = 1
      Z = z
    while(nz != ncol(Ab_move) + 1){ ##Simulates until it gets out of plot
        v = plogis(params[[1]][nid,6] +params[[1]][nid,7] * B[z]/100)#B distance-to-edge vector
        pm <- like_mov(distance=Distance[z,], covers=open, pl = plot_v,fruits=log(nFj[week,]+1), ab_mv = Ab_move[week,] ,
                        a_d =params[[1]][nid,1], b_c=params[[1]][nid,2] , 
                        a_fs=params[[1]][nid,3],b_f=params[[1]][nid,4],
                         b_s=params[[1]][nid,5]  
                       )#Movement function
        nz = sample(pids,size=1,prob=pm) #Sample next cell
       if(runif(1)<v){
          nz =(ncol(Ab_move) + 1)
               }
        z = nz
        count = count+1 #Count sirve para ir contanto el número de movimientos que hace antes de marcharse fuera
        Z[count] = z  #Guarda en el vector Z a qué celda se ha movido

      }
      Bird = c(Bird, rep(nid,count)) #identify species with steps
      Zs = c(Zs, Z)#tracks the steps of speceis
      Track = c(Track, rep(i, count))#identify of track
    } 
  
    Bird <- Bird[-1] ##Remove first data (randomly sampled)
    Zs <- Zs[-1]
    Track = Track[-1]
    
    
# Quantify distance of steps    
   R <- numeric((length(Zs)-1))*NA #
    for(i in 1:length(R)){
      if(Track[i]==Track[i+1]) {R[i] <- distos[Zs[i],Zs[i+1]] }
    }
    
    
#########End movement simulation ###########################
  
   
##############Fruit consumption####################
    consumed = rep(NA, length(Zs))  
    for(i in 1:length(Zs))
    {
      if(Zs[i] == ncol(Ab_move)+1)
      {consumed[i] = NA #out of plot
      }else{
        consumed[i] = rzipois(1, params[[3]][Bird[i],1] , params[[3]][Bird[i],2])
        }
      }
   
    idf <- numeric(length(Zs))  # Identity of fruits consumed
   #idf ==0 no consumption y NA out-of-plot movement
     fsp = 1:ncol(Fc[,,1])   # indicator for plant species

     for(i in 1:length(Zs)){
      if(Zs[i] == (ncol(Ab_move)+1)){
      idf[i]  = NA #p
      }else{
      tm <- nFj[week, Zs[i]] #  number of fruits (just to see if there are fruits in the cell)
      tm2 <-  length(Fc[Zs[i],,week])-length(which(is.na(Fc[Zs[i],,week]))) 
    #if the bird chose to consume, cell holds fruits y and not NA (some cells had information missing in certain weeks)
      if(consumed[i] > 0 & tm > 0 & is.finite(tm) & tm2 > 0){  #
         ab = as.numeric(Fc[Zs[i],,week])
         Ac <- plogis(params[[2]][Bird[i],1]+ params[[2]][Bird[i],2]*ab)*Forb_int[Bird[i],]#Fruit selection function
         loss_spp = which(ab == 0)
         Ac[loss_spp] = 0 #Skip a problem with dataset(unidentified fruits)
         if(sum(Ac) == 0)#
         {
           consumed[i] =0
           idf[i]=0
         }else{
           k = Ac /sum(Ac)
           idf[i] <- sample(fsp,1,prob=k) #Sample the fruits consumed 
         }
         # 
       ##Assess real number of fruits consumed (according to cell availability)
        consumed[i] <- min(consumed[i], ceiling(Fc[Zs[i],idf[i],week]*nFj[week,Zs[i]]))
        nFj[week, Zs[i] ] = nFj[week, Zs[i] ] - consumed[i] #update number of fruits in cell
        }else{
        idf[i] = 0
        consumed[i] = 0}#End of ifelse fruit consumtion
    }#end if (Zs[1] != 441)
  }#end fruit consumption iteration


 #setdiff(which(consumed == 0), which(idf == 0)) # Check ok
     
#Save vectors for movement comparisons

WW = rep(week, length(idf))


# Simulate seed defecation ------------------------------------------------


    flightspeed = 360 # meters per minute (see refs in Morales & Carlo 2006)
    TFl <- (R*10)/flightspeed # time flying (minutes). 10 m if flight R = 0 (same cell, cell size)
   TFl[TFl==0] <- 10/flightspeed #Same cell
   TFl[is.na(TFl)] <-0  #Out of plot
   TFl <- c(TFl,0) #Add one 0 in order to have the same length between time flying and consumed fruits 
#remember that time flying since they are jumps between cells the have ncells used -1
#Perching time and GPT
     PT = rep(NA, length(Bird))
    for(i in 1:length(Bird))
    {
       if(Zs[i] != (ncol(Ab_move)+1))
       {tmp = rgamma(1, shape = params[[4]][Bird[i],1], rate = params[[4]][Bird[i],2])
        PT[i] = tmp/60 #transform into minutes
       } 
       }
    PT[is.na(PT)] <-0 #out-of-plot
     
  #  View(cbind(Track, Bird, Zs, TFl, PT, consumed, idf))#Check ok

#Time in Gut
    TinG <- rgamma(length(Zs), scale=scale_GPT[Bird[i]],shape=1.59)  # time for seeds in gut (in minutes)
    ###Time in Gut same for all seeds consumed at the same time (not individualized per seed)
  
    Zc = numeric(length(Zs)) # will hold location of defecated seeds
    No = numeric(length(Zs)) # will hold seeds out of plot
    K =  numeric(length(Zs))*NA # will hold dispersal distance
    Zm = matrix(NA,length(Zs),5) #guarda las coordenadas del movimiento (celda inicial y celda final), spp planta, spp pájaro, frutos consumidos
    sem <- consumed*idf/idf  # actual number of seeds consumed (to remove those with idf = 0, no fruits consumed) 
    #setdiff(which(sem == 0), which(idf == 0)) Ok
    sem[is.na(sem)] <- 0    # set divided y zero equal to zero
    
##Pooping
    #(A) same cell as consumption
    #(B) Different cell: while flying/perching
    #(C) Out of plot
  
    tmp = c()
    tmp2 = c()
   # tmp3 = c()
    for(i in 1:length(Zc)){
      if(sem[i]>0){
            idt <- which(Track==Track[i] & Zs != 0)  # find which observations correspond to the current track
             TTime <- PT[i:max(idt)] + TFl[i:max(idt)]  # add perching time and flying time throughout the track
        Tcum <- cumsum(TTime) #Track clock "cumulative time of movements within tracks"
       ##Has it pooped seeds throughout the track
         if(TinG[i]<max(Tcum)){ ## 
       #   tmp = c(tmp,1)
          idc <- which(log(Tcum-TinG[i])== min(log(Tcum-TinG[i]), na.rm=T)) 
          idc = idc[[1]] ##avoid problems with out-of-plot data (length(idc = 2))
        #  if(length(idc)> 1)
        #  {tmp3 = c(tmp3, i)}
         #which is the closest largest value to TinG (idc = 1, same cell, idc> 1 different cell)
          Zc[i] <- Zs[(i-1)+(idc)] 
        #which cell should it poop (i-1 because it TFl represent jumps)
          K[i] <- distos[Zs[i],Zc[i]] #calcula la distancia de dispersion
          Rain[Zc[i], idf[i]] = Rain[Zc[i], idf[i]] + consumed[i]
          if(is.na(K[i]))#pooping out-of-plot (end of track not pooped yet)
          {tmp2 = c(tmp2, 1)}
          Zm[i,] <- c(Zs[i], Zc[i], idf[i],Bird[i], consumed[i])
##poops in another cell that is not the same as consumption 
            if(idc>1){
            #if it poops while flying
            if( (TinG[i]-Tcum[idc-1]) > PT[(i-1)+(idc)]){ # 
              delta <- TinG[i]-Tcum[idc-1]-PT[(i-1)+(idc)] # time the seed is in flight
              speedx <- ( xp[Zs[(i)+(idc)]] - xp[Zs[(i-1)+(idc)]] )/TFl[(i-1)+(idc)] #coordinate x for deposition cell
              speedy <- ( yp[Zs[(i)+(idc)]] - yp[Zs[(i-1)+(idc)]] )/TFl[(i-1)+(idc)] #coordinate y for deposiition cell
              xn <- round( xp[Zs[(i-1)+(idc)]] + delta*speedx ) #p
              yn <- round( yp[Zs[(i-1)+(idc)]] + delta*speedy ) 
              Zc[i] <- which(xp==xn & yp==yn)#identify target cell
              K[i] <- distos[Zs[i],Zc[i]] #seed dispersal distance
              Rain[Zc[i], idf[i]] = Rain[Zc[i], idf[i]] + consumed[i]#Seed rain matrix
              Zm[i,] <- c(Zs[i], Zc[i], idf[i],Bird[i], consumed[i])#Information of movement coordinates
            }
          }#termina el idc > 1
        
        }else{
          No[i] <- 1
          Zm[i,] <- c(Zs[i],NA, idf[i],Bird[i], consumed[i])
          }
        } #end of if de sem>0 section (i.e., seeds consumed)
    }#end of Zc loop (Zc had info of tracks)
    
    #Checks
    #length(which(sem>0))# las que han comido y no han sido descartadas por ser los primeros tracks
   #  sum(No) #cagados afuera del track
   #  sum(tmp2)#cagados en el borde del plot Z = 401
   #  length(which(!is.na(K))) #los que han cagado dentro del plot
  
    
##Safe weeky outputs
 ##Consumption and dispersal outputs
    K_w = c(K_w, K[which(sem>0)] )#Seed dispersal distances
    Zm_w = c( Zm_w,Zm[which(sem>0),])#Coodinates of seed movement
    idf_w = c(idf_w, idf[which(sem>0)]) # Identity of fruit spececies
    consumed_w = c(consumed_w, consumed[which(sem>0)])#Number of fruits consumed
    B_w = c(B_w, Bird[which(sem>0)])#Bird identifier
    No_w = c(No_w, sum(No))#Number of seeds out of plot
    Week_w = c(Week_w, WW[which(sem>0)])#week identifier
    Track_s = Track[which(sem>0)]
    if(week > 1)
    {Track_s = Track_s + max(Track_s) 
    }else
    {
      Track_s = Track_s
    }
    Track_w = c(Track_w, Track_s)
    
##Save information of bird movements (with and without consumption)
    Zs_mov_w = c(Zs_mov_w, Zs[which(Zs != 0)]) #
    R_mov_w = c(R_mov_w,R[which(Zs != 0)])
    BB_mov_w = c(BB_mov_w, Bird[which(Bird != 0)])
    Week_mov_w = c(Week_mov_w, WW[which(!is.na(WW))]) #Identificador de las semanas
   
 #Renumber tracks to avoid repetitions
    Track_m = Track[!is.na(Track)]
    if(week > 1)
    {Track_m = Track_m + max(Track_m) 
    }else
    {
      Track_m = Track_m
    }
    Track_mov_w = c(Track_mov_w, Track_m)
 
  
   ##Ourputs related to seed diversity
    out = cbind(idf[which(sem>0)], consumed[which(sem>0)],Bird[which(sem>0)])
    out = data.frame(out)
    colnames(out) = c("species", "N", "birds")
    fr_partes = colnames(fr_pico)[-1]
    
    disp = rep(0, length(plant_ids))
    for (k in 1:length(plant_ids))
    {
      tmp = out[which(out$species == k),]
      if(nrow(tmp) > 0)
      {
        if(plants_name[k]%in% fr_partes)
        {
          tmp_disp = 0
          for(ss in 1:length(bsp))
          {
            tmp3 = tmp[which(tmp$birds == ss),]
            if(nrow(tmp3))
            {
              tmp_nseeds = fr_pico[ss,which(colnames(fr_pico) == plants_name[k])]
              tmp_disp = tmp_disp + (sum(tmp3$N)*tmp_nseeds)
            }
          }
          disp[k] = tmp_disp
        }else{
          #No se comen por partes así que el fruto entero * semillas
          disp[k] = sum(tmp$N)*ceiling(Nseeds[k])
        }
        
      }
      
    }
    
    H0 = av_div[week]
    H1 = diversity(disp, "shannon")
    expH0 = exp(H0)
    expH1 = exp(H1)
    
    S0<- av_Ri[week]
    J0 <- H0/log(S0)
    S1<- specnumber(disp)
    J1 <- H1/log(S1)
    
    
    delta_H = (expH1-expH0)/expH0
    delta_S  = S1 - S0
    delta_J = (J1-J0)/J0
    
    deltaH_w = c(deltaH_w, delta_H)
    deltaS_w= c(deltaS_w, delta_S)
    deltaJ_w = c(deltaJ_w, delta_J)
    
    #===============================================================================
    
  }##End of week loop
    #Guarda la info of repetition (j)
  #Fruit consumption & dispersal
    KK[[j]] <- K_w #dispersal distances
    ZZ[[j]] <- Zm_w#Movement and consumption coordinates. KEY output Interactions
    IDC[[j]] <- idf_w # Identity of fruit spp
    NC[[j]] <- consumed_w #Number of fruits consumed
    BB[[j]] <- B_w#Identity of bird spp
    WEEK[[j]] <- Week_w #week
    Tr[[j]] = Track_w
    Nos <- c(Nos,No_w)#seeds out of plot
  #Movement
    ZZ_mov[[j]] <- Zs_mov_w
    RR_mov[[j]] <- R_mov_w 
    BB_mov[[j]] <- BB_mov_w
    WEEK_mov[[j]] <- Week_mov_w
    Tr_mov[[j]] <- Track_mov_w
#Changes in diversity, equitativity...
    Div[[j]] <- deltaH_w
    Ri[[j]] <- deltaS_w
    J[[j]] <- deltaJ_w
###summaries
    Mov[[j]] <- cbind(Week_mov_w, Track_mov_w, BB_mov_w, Zs_mov_w, R_mov_w)
    Inter[[j]] = cbind(Week_w, B_w, idf_w, consumed_w)
    Seedrain[[j]] = Rain
    
 ###SAVE DATA EVERY 10 REPS
    
  if(j%%10 == 0)
  {save.image("Llao_simulation.RData")}
    
  }##FIN DEL BUCLE nreps
save.image("Llao_simulation_all_community.RData")