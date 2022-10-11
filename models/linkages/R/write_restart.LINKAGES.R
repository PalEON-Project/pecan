#-------------------------------------------------------------------------------
# Copyright (c) 2012 University of Illinois, NCSA.
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

##' @title write_restart.LINKAGES
##' @name  write_restart.LINKAGES
##' @author Ann Raiho \email{araiho@@nd.edu}
##' 
##' @param outdir      output directory
##' @param runid       run ID
##' @param time        year that is being read
##' @param settings    PEcAn settings object
##' @param new.state    analysis vector
##' @param RENAME      flag to either rename output file or not
##' @param variables
##' @param sample_parameters
##' @param trait.values
##' 
##' @description Write restart files for LINKAGES
##' 
##' @return NONE
##' @export
##' 

# outdir, runid, time, settings, new.state, variables, sample_parameters = FALSE, trait.values =
# NA,met=NULL,RENAME = TRUE

write_restart.LINKAGES <- function(outdir, runid, start.time, stop.time,
                                   settings, new.state, 
                                   RENAME = TRUE, new.params, inputs) {
  
  ### TO DO : needs to be vectorized to improve SDA speed for runs that are longer than 50 years
  ### TO DO : distance matrix needs fixing
  
  ### Removing negative numbers because biomass can't be negative ###
  new.state[new.state < 0] <- 0
  
  names.keep <- names(new.state)
  new.state <- as.matrix(new.state)
  names(new.state) <- names.keep
  
  new.state.save <- new.state
  
  if(any(grep('Fcomp',names.keep))){
    new.state <- new.state.save[grep("Fcomp", names(new.state.save))]
    new.state.other <- new.state.save[grep("Fcomp", names(new.state.save), invert = TRUE)]
  }
  
  if(any(grep('AGB.pft',names.keep))){
    new.state <- new.state.save[grep("AGB.pft", names(new.state.save))]
    new.state.other <- new.state.save[grep("AGB.pft", names(new.state.save), invert = TRUE)]
  }
  
  variables <- names(new.state)
  ### Going to need to change this... ### Get some expert opinion
  N <- length(new.state)
  #distance.matrix <- matrix(1, N, N)
  #for (i in seq_len(N)) {
  #  distance.matrix[i, ] <- sample(c(seq(0, N-1, 1)), size = N)
  #  if(which(distance.matrix[i,]==0)!=i){
  #    distance.matrix[i,which(distance.matrix[i,]==0)] <- distance.matrix[i,i]
  #    distance.matrix[i,i] <- 0
  #  } 
  #}
  #diag(distance.matrix) <- 0
  
  if(FALSE){
    distance.matrix <- rbind(c(0, 1, 4, 3, 2, 6, 5, 8, 7, 9, 10, 11, 12, 13, 14), 
                             c(5, 0, 3, 4, 8, 1, 2, 7, 6, 9, 10, 11, 12, 13, 14), 
                             c(5, 3, 0, 1, 8, 4, 2, 7, 6, 9, 10, 11, 12, 13, 14), 
                             c(6, 2, 1, 0, 8, 4, 3, 7, 5, 9, 10, 11, 12, 13, 14), 
                             c(2, 7, 5, 4, 0, 8, 6, 1, 3, 9, 10, 11, 12, 13, 14), 
                             c(6, 1, 3, 4, 8, 0, 2, 7, 5, 9, 10, 11, 12, 13, 14), 
                             c(5, 3, 1, 2, 8, 6, 0, 7, 4, 9, 10, 11, 12, 13, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 0, 2, 9, 10, 11, 12, 13, 14), 
                             c(1, 5, 3, 2, 7, 6, 4, 8, 0, 9, 10, 11, 12, 13, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 9, 2, 0, 10, 11, 12, 13, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 10, 2, 9, 0, 11, 12, 13, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 11, 2, 9, 10, 0, 12, 13, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 12, 2, 9, 10, 11, 0, 13, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 13, 2, 9, 10, 11, 12, 0, 14), 
                             c(3, 6, 4, 5, 1, 7, 8, 14, 2, 9, 10, 11, 12, 13, 0))
    
  }
  #distance.matrix <- rbind(c(0,3,1,2), c(3,0,2,1), c(1,2,0,3), c(2,1,3,0))
  
  ## HACK
  #spp.params.default <- read.csv(system.file("spp_matrix.csv", package = "linkages"))  #default spp.params
  spp.params.default <- read.csv("/data/dbfiles/spp_matrix.csv")  # default spp.params
  nspec <- length(settings$pfts)
  spp.params.save <- numeric(nspec)
  for (i in seq_len(nspec)) {
    spp.params.save[i] <- which(spp.params.default[, 1] %in% settings$pfts[i]$pft$name)
  }
  
  spp.params <- spp.params.default[spp.params.save, ]
  biomass_spp_params <- function(new.params, default.params, pft) {
    if ("SLTA" %in% names(new.params[[as.character(pft)]])) {
      slta <- new.params[[as.character(pft)]]$SLTA
    } else {
      slta <- default.params[default.params$Spp_Name == pft, ]$SLTA
    }
    if ("SLTB" %in% names(new.params[[as.character(pft)]])) {
      sltb <- new.params[[as.character(pft)]]$SLTB
    } else {
      sltb <- default.params[default.params$Spp_Name == pft, ]$SLTB
    }
    if ("SLA" %in% names(new.params[[as.character(pft)]])) {
      sla_use <- (1/new.params[[as.character(pft)]]$SLA)*1000
      sla_use[sla_use > 5000] <- rnorm(1,4000,100)
      fwt <- sla_use  #(1 / new.params[[as.character(pft)]]$SLA) * 10000
    } else {
      fwt <- default.params[default.params$Spp_Name == pft, ]$FWT
    }
    if ("FRT" %in% names(new.params[[as.character(pft)]])) {
      frt <- new.params[[as.character(pft)]]$FRT
    } else {
      frt <- default.params[default.params$Spp_Name == pft, ]$FRT
    }
    if(pft == 'ACRU'){b_0 = -2.047; b_1 = 2.3852} # acer rubrum, red maple
    if(pft == 'ACSA3'){b_0 = -1.8011; b_1 = 2.3852} # acer saccharum, sugar maple
    if(pft == 'BEAL2'){b_0 = -1.8096; b_1 = 2.348} # betula alleghaniensis, yellow birch
    if(pft == 'BELE'){b_0 = -2.2652; b_1 = 2.5349} # betula lenta, sweet birch
    if(pft == 'FAGR'){b_0 = -2.0705; b_1 = 2.441} # fagus grandifolia, american beech
    if(pft == 'PIRU'){b_0 = -2.1364; b_1 = 2.3233} # picea rubens, red spruce
    if(pft == 'PIST'){b_0 = -2.6177; b_1 = 2.4638} # pinus strobus, white pine
    if(pft == 'QUAL'){b_0 = -2.0705; b_1 = 2.441} # quercus alba, white oak
    if(pft == 'QUMO'){b_0 = -2.0705; b_1 = 2.441} # quercus montana, chestnut oak
    if(pft == 'QURU'){b_0 = -2.0705; b_1 = 2.441} # quercus rubra, red oak
    if(pft == 'QUVE'){b_0 = -2.0705; b_1 = 2.441} # quercus velutina, black oak
    if(pft == 'THOC2'){b_0 = -1.9615; b_1 = 2.1063} # thuja occidentalis, northern white-cedar
    if(pft == 'TSCA'){b_0 = -2.348; b_1 = 2.3876} # tsuga canadensis, eastern hemlock
  # acer rubrum, red maple
  if(pft == 'ACRU'){bwood1 = 0.1014; bwood2 = 2.3448; bbark1 = 0.0291; bbark2 = 2.0893;
  bbranches1 = 0.0175; bbranches2 = 2.4846; bfoliage1 = 0.0515; bfoliage2 = 1.5198}
  # acer saccharum, sugar maple
  if(pft == 'ACSA3'){bwood1 = 0.1315; bwood2 = 2.3129; bbark1 = 0.0631; bbark2 = 1.6241;
  bbranches1 = 0.033; bbranches2 = 2.3741; bfoliage1 = 0.0393; bfoliage2 = 1.693}
  # betula alleghaniensis, yellow  birch
  if(pft == 'BEAL2'){bwood1 = 0.1932; bwood2 = 2.1569; bbark1 = 0.0192; bbark2 = 2.2475;
  bbranches1 = 0.0305; bbranches2 = 2.4044; bfoliage1 = 0.1119; bfoliage2 = 1.3973}
  # fagus grandifolia, american beech
  if(pft == 'FAGR'){bwood1 = 0.1478; bwood2 = 2.2986; bbark1 = 0.012; bbark2 = 2.2388;
  bbranches1 = 0.037; bbranches2 = 2.368; bfoliage1 = 0.0376; bfoliage2 = 1.6164}
  # picea rubens, red spruce
  if(pft == 'PIRU'){bwood1 = 0.0989; bwood2 = 2.2814; bbark1 = 0.022; bbark2 = 2.0908;
  bbranches1 = 0.0005; bbranches2 = 3.275; bfoliage1 = 0.0066; bfoliage2 = 2.4213}
  # pinus strobus, white pine
  if(pft == 'PIST'){bwood1 = 0.0997; bwood2 = 2.2709; bbark1 = 0.0192; bbark2 = 2.2038;
  bbranches1 = 0.0056; bbranches2 = 2.6011; bfoliage1 = 0.0284; bfoliage2 = 1.9375}
  # quercus alba, white oak
  if(pft == 'QUAL'){bwood1 = 0.0762; bwood2 = 2.3335; bbark1 = 0.0338; bbark2 = 1.9845;
  bbranches1 = 0.0113; bbranches2 = 2.6211; bfoliage1 = 0.0188; bfoliage2 = 1.7881}
  # quercus rubra, red oak
  if(pft == 'QURU'){bwood1 = 0.1754; bwood2 = 2.1616; bbark1 = 0.0381; bbark2 = 2.0991;
  bbranches1 = 0.0085; bbranches2 = 2.779; bfoliage1 = 0.0373; bfoliage2 = 1.674}
  # tsuga canadensis, eastern hemlock
  if(pft == 'TSCA'){bwood1 = 0.0619; bwood2 = 2.3821; bbark1 = 0.0139; bbark2 = 2.3282;
  bbranches1 = 0.0217; bbranches2 = 2.2653; bfoliage1 = 0.0776; bfoliage2 = 1.6995}
  # betula lenta, sweet birch (used values for hop-hornbeam, classification from chojnacky)
  if(pft == 'BELE'){bwood1 = 0.1929; bwood2 = 1.9672; bbark1 = 0.0671; bbark2 = 1.5911;
  bbraches1 = 0.0278; bbranches2 = 2.1336; bfoliage1 = 0.0293; bfoliage2 = 1.9502}
  # quercus montana, chestnut oak (used values for white oak, classification from chojnacky)
  if(pft == 'QUMO'){bwood1 = 0.0762; bwood2 = 2.3335; bbark1 = 0.0338; bbark2 = 1.9845;
  bbranches1 = 0.0113; bbranches2 = 2.6211; bfoliage1 = 0.0188; bfoliage2 = 1.7881}
  # quercus velutina, black oak (used values for white oak, classification from chojnacky)
  if(pft == 'QUVE'){bwood1 = 0.0762; bwood2 = 2.3335; bbark1 = 0.0338; bbark2 = 1.9845;
  bbranches1 = 0.0113; bbranches2 = 2.6211; bfoliage1 = 0.0188; bfoliage2 = 1.7881}
  # thuja occidentalis, northern white-cedar
  if(pft == 'THOC2'){bwood1 = 0.0654; bwood2 = 2.2121; bbark1 = 0.0114; bbark2 = 2.1432;
  bbranches1 = 0.0335; bbranches2 = 1.9367; bfoliage1 = 0.0499; bfoliage2 = 1.7278}
    
  return(list(slta = slta, sltb = sltb, fwt = fwt, frt = frt, 
              b_0 = b_0, b_1 = b_1,
              bwood1 = bwood1, bwood2 = bwood2, bbark1 = bbark1, bbark2 = bbark2,
              bbranches1 = bbranches1, bbranches2 = bbranches2, bfoliage1 = bfoliage1,
              bfoliage2 = bfoliage2))
  } # biomass_spp_params
  
  biomass_function <- function(dbh, spp.biomass.params) {
    # Default
    # kg/tree
    if(bio_method == 'default'){0.1193 * dbh ^ 2.393 + 
      ((spp.biomass.params$slta + spp.biomass.params$sltb * dbh) / 2) ^ 2 * 
      3.14 * spp.biomass.params$fwt * spp.biomass.params$frt * 0.001}
    # Chojnacky equations
    if(bio_method == 'chojnacky'){(exp(spp.biomass.params$b_0)*(dbh^spp.biomass.params$b_1))}
    # Lambert equations
    if(bio_method == 'lambert'){
      spp.biomass.parameters$bwood1*(dbh^spp.biomass.parameters$bwood2) +
      spp.biomass.parameters$bbark1*(dbh^spp.biomass.parameters$bbark2) +
      spp.biomass.parameters$bfoliage1*(dbh^spp.biomass.parameters$bfoliage2) +
      spp.biomass.parameters$bbranches1*(dbh^spp.biomass.parameters$bbranches2)}
  } # biomass_function
  
  merit <- function(dbh, b_obs, spp.biomass.params) {
    (b_obs - biomass_function(dbh, spp.biomass.params)) ^ 2
  } # merit

  ## distance matrix calculation :: identify ranking for cloning species by identifying which species are closer together in parameter space
  # gather all species parameters into a dataframe 
  all.params = spp.params
  
  for (pft in spp.params$Spp_Name){
      
    # get information for specific PFT 
    pft.ind = which(all.params$Spp_Name == pft)
    available.vals = names(new.params[[as.character(pft)]])
      
    # overwrite default parameters where prior-drawn parameter are available
      
    # MK: if original is over 5000, this should probably be the same value drawn in the other part of the script
    if ("SLA" %in% available.vals) {
      sla_use <- (1/new.params[[as.character(pft)]]$SLA)*1000
      sla_use[sla_use>5000] <- rnorm(1,4000,100)
      all.params$FWT[pft.ind] <- sla_use
    }
    if ("HTMAX" %in% available.vals & "DBHMAX" %in% available.vals) {
      all.params$B2[pft.ind] <- 2 * (((new.params[[as.character(pft)]]$HTMAX * 100) - 137) / 
                                       (new.params[[as.character(pft)]]$DBHMAX * 100))
      all.params$B3[pft.ind] <- (new.params[[as.character(pft)]]$HTMAX * 100 - 137) / (new.params[[as.character(pft)]]$DBHMAX * 100)^2
    }
    if ("root2shoot" %in% available.vals) {
      all.params$RTST[pft.ind] <- new.params[[as.character(pft)]]$root2shoot
    }
    if ("DMAX" %in% available.vals) {
      all.params$DMAX[pft.ind] <- new.params[[as.character(pft)]]$DMAX
    }
    if ("DMIN" %in% available.vals) {
      all.params$DMIN[pft.ind] <- new.params[[as.character(pft)]]$DMIN
    }
    if ("AGEMX" %in% available.vals) {
      all.params$AGEMX[pft.ind] <- new.params[[as.character(pft)]]$AGEMX
    }
    if ("Gmax" %in% available.vals) {
      all.params$G[pft.ind] <- new.params[[as.character(pft)]]$Gmax
    }
    if ("SPRTND" %in% available.vals) {
      all.params$SPRTND[pft.ind] <- new.params[[as.character(pft)]]$SPRTND
    }
    if ("SPRTMN" %in% available.vals) {
      all.params$SPRTMN[pft.ind] <- new.params[[as.character(pft)]]$SPRTMN
    }
    if ("SPRTMX" %in% available.vals) {
      all.params$SPRTMX[pft.ind] <- new.params[[as.character(pft)]]$SPRTMX
    }
    if ("MPLANT" %in% available.vals) {
      all.params$MPLANT[pft.ind] <- new.params[[as.character(pft)]]$MPLANT
    }
    if ("D3" %in% available.vals) {
      all.params$D3[pft.ind] <- new.params[[as.character(pft)]]$D3
    }
    if ("FROST" %in% available.vals) {
      all.params$FROST[pft.ind] <- new.params[[as.character(pft)]]$FROST
    }
    if ("CM1" %in% available.vals) {
      all.params$CM1[pft.ind] <- new.params[[as.character(pft)]]$CM1
    }
    if ("CM2" %in% available.vals) {
      all.params$CM2[pft.ind] <- new.params[[as.character(pft)]]$CM2
    }
    if ("CM3" %in% available.vals) {
      all.params$CM3[pft.ind] <- new.params[[as.character(pft)]]$CM3
    }
    if ("CM4" %in% available.vals) {
      all.params$CM4[pft.ind] <- new.params[[as.character(pft)]]$CM4
    }
    if ("CM5" %in% available.vals) {
      all.params$CM5[pft.ind] <- new.params[[as.character(pft)]]$CM5
    }
    if ("SLTA" %in% available.vals) {
      all.params$SLTA[pft.ind] <- new.params[[as.character(pft)]]$SLTA
    }
    if ("SLTB" %in% available.vals) {
      all.params$SLTB[pft.ind] <- new.params[[as.character(pft)]]$SLTB
    }
    if ("FRT" %in% available.vals) {
      all.params$FRT[pft.ind] <- new.params[[as.character(pft)]]$FRT
    }
  }
    
  # remove all parameters not to be used in distance calculation 
  # TL is a categorial variable so doesn't make sense to consider for distance 
  remove.ids <- which(names(all.params) %in% c('Spp_Name', 'Spp_Number', 'TL'))
  all.params <- all.params[,-c(remove.ids)]
  rownames(all.params) <- spp.params$Spp_Name
  
  # calculate distance matrix of standardized parameter matrix 
  distances <- as.matrix(dist(scale(all.params), method = 'euclidean',upper = TRUE, diag=TRUE))
  distance.matrix <- distances
  
  # place similarity rankings in rows of distance.matrix for each species
  for (i in 1:nrow(distances)){
    ord <- sort(distances[i,], index.return=TRUE)$ix
    distance.matrix[i,] = sapply(c(1:ncol(distances)), function(col.ind){which(ord == col.ind)-1})
  }
  

  ## HACK
  
  # skip ensemble member if no file availible
  outfile <- file.path(outdir, runid, "linkages.out.Rdata")
  if (!file.exists(outfile)) {
    outfile <- file.path(outdir, runid, paste0(start.time, "linkages.out.Rdata"))
    if (!file.exists(outfile)) {
      logger.severe(paste0("missing outfile ens #", runid))
    }
  }
  print(paste0("runid = ", runid))
  
  # load output
  load(outfile)
  
  ntrees <- ntrees.kill[, ncol(ntrees.kill), 1]  # number of trees
  
  if(sum(ntrees)==0) {
    #reloads spin up if theres nothing in the output file
    print('No survivors. Reusing spinup.')
    load(file.path(outdir, runid,list.files(file.path(outdir, runid))[grep(list.files(file.path(outdir, runid)),pattern='linkages')][1]))
    ntrees <- ntrees.kill[, ncol(ntrees.kill), 1]  # number of trees
    
  }
  
  nspec  <- length(settings$pfts)
  ncohrt <- ncohrt
  tyl    <- tyl
  C.mat  <- C.mat
  
  nogro  <- as.vector(nogro.save[, ncol(nogro.save), 1])  ## no growth indicator
  ksprt  <- matrix(0, 1, nspec)  ## kill sprout indicator ## LOOK INTO THIS
  iage   <- as.vector(iage.save[, ncol(iage.save), 1])  # individual age
  
  dbh    <- as.vector(dbh.save[, ncol(dbh.save), 1])
  
  n.index <- c(rep(1, ntrees[1]))
  for (i in 2:length(settings$pfts)) {
    n.index <- c(n.index, rep(i, ntrees[i]))
  }
  
  #if(max(dbh) < 20){ # if all trees are small than large trees are 95th percentile otherwise trees bigger than 5 cm
    #large.trees <- which(dbh >= (max(dbh) / 1.05))
  #}else{
    #large.trees <- which(dbh >= 20)
  #}
  
  #large.trees <- which(dbh > 0)

  # following current assumptions for biomass data processing, only biomass for trees with DBH greater than 5 cm is considered in the data 
  large.trees <- which(dbh > 5)
  if (length(large.trees) == 0) large.trees <- which(dbh >= quantile(dbh, 0.95))
  
  for (s in seq_along(settings$pfts)) {
    ntrees[s] <- length(which(n.index[large.trees] == s))
  }
  
  n.index <- n.index[large.trees]
  
  dbh <- dbh[large.trees]
  iage <- iage[large.trees]
  nogro <- nogro[large.trees]
  
  new.ntrees <- numeric(length(settings$pfts))
  
  print(paste0("ntrees (large trees) =", ntrees))  #these are the large trees
  
  ##### This takes the average individual biomass of each species from the model and computes how many
  ##### individuals you should keep to match the biomass estimated from the data.  Still have to correct
  ##### for the total species biomass in the next step.
  
  ind.biomass <- numeric(sum(ntrees))
  
  # calculate biomass of each individual
  for (j in seq_len(sum(ntrees))) {
    # slta <- spp.params$SLTA[n.index[j]] sltb <- spp.params$SLTB[n.index[j]] fwt <-
    # spp.params$FWT[n.index[j]] frt <- spp.params$FRT[n.index[j]]
    pft <- spp.params$Spp_Name[n.index[j]]
    spp.biomass.params <- biomass_spp_params(new.params = new.params, 
                                             default.params = spp.params.default, 
                                             pft = pft)
    ind.biomass[j] <- biomass_function(dbh[j], spp.biomass.params) * (1 / 833) * 0.48  # changing units to be kgC/m^2
  }
  
  data2 <- data.frame(ind.biomass = ind.biomass,
                      n.index = n.index)
  mean.biomass.spp <- aggregate(ind.biomass ~ n.index, mean, data = data2)   # calculate mean individual biomass for each species
  #browser()
  # calculate number of individuals needed to match new.state
  for (s in seq_along(settings$pfts)) {
    
    if (ntrees[s] > 0) {
      fix_adjust <- new.state[s]/mean.biomass.spp[mean.biomass.spp[, 1] == s, 2]  # number of individuals needed to agree with new.state      
    } else {
      for (r in 1:(length(settings$pfts) - 1)) {
        s.select <- which(distance.matrix[s, ] == r)  # select a new spp. to clone from
        if (ntrees[s.select] > 0) {
          break
        }
      }
      fix_adjust <- new.state[s] / mean.biomass.spp[mean.biomass.spp[, 1] == s.select, 2]
    }
    new.ntrees[s] <- as.numeric(ceiling(fix_adjust-.01))  #new number of ind. of each species
    if(new.ntrees[s]>200&!is.na(new.ntrees[s])){
      new.ntrees[s] = sample(size = 1, x = 50:150)
    } 
    print(s)
  }
  
  #making sure to stick with density dependence rules in linkages (< 198 trees per 800/m^2)
  #someday we could think about estimating this parameter from data
  if(sum(new.ntrees,na.rm = T) > 198) new.ntrees <- round((new.ntrees / sum(new.ntrees)) * runif(1,195,198))
  
  print(paste0("new.ntrees =", new.ntrees))
  
  new.n.index <- c(rep(1, new.ntrees[1]))
  for (i in 2:length(settings$pfts)) {
    new.n.index <- c(new.n.index, rep(i, new.ntrees[i]))
  }
  
  n.ind <- 200
  
  dbh.temp <- numeric(n.ind)
  iage.temp <- numeric(n.ind)
  nogro.temp <- numeric(n.ind)
  
  # sample from individuals to construct new states
  for (s in seq_len(nspec)) {
    if (new.ntrees[s] == 0) {
      next
    }
    if (new.ntrees[s] <= ntrees[s]) {
      # new are less than the old of the same spp.  print('new are less than the old of the same spp.')
      select <- sample(size = new.ntrees[s], x = which(n.index == s), replace = FALSE)
    } else {
      if (new.ntrees[s] > ntrees[s] & ntrees[s] >= 1) {
        # new are greater than the old of the same spp. and there are old trees to clone print('new are
        # greater than the old of the same spp. and there are old trees of same spp. to clone')
        select <- c(which(n.index == s), 
                    sample(size = (new.ntrees[s] - ntrees[s]), x = which(n.index == s), replace = TRUE))
      } else {
        # print(paste0('clone needed for spp. ',s))
        for (r in 1:(length(settings$pfts) - 1)) {
          s.select <- which(distance.matrix[s, ] == r)  #select a new spp. to clone from
          # print(paste0('r =',r))
          if (ntrees[s.select] > 0) {
            break
          }
        }
        # print(s.select)
        select <- sample(size = as.numeric(new.ntrees[s]), 
                         x = which(n.index == s.select), 
                         replace = TRUE)
      }
    }
    dbh.temp[which(new.n.index == s)] <- dbh[select]
    iage.temp[which(new.n.index == s)] <- iage[select]
    nogro.temp[which(new.n.index == s)] <- nogro[select]
  }
  
  # fix dbh of sampled individuals to match new.state
  nl <- 1  ## individual counter
  b_calc <- numeric(length(settings$pfts))  #biomass of sampled trees
  b_calc1 <- numeric(length(settings$pfts))  #biomass of sampled trees
  bcorr <- numeric(length(settings$pfts))  #biomass correction factor to new.state
  b_obs <- numeric(sum(new.ntrees))
  for (s in seq_len(nspec)) {
    if (new.ntrees[s] == 0) {
      next
    }
    nu <- nl + new.ntrees[s] - 1
    pft <- unique(spp.params$Spp_Name[new.n.index[nl:nu]])
    spp.biomass.params <- biomass_spp_params(new.params = new.params, 
                                             default.params = spp.params.default, 
                                             pft = pft)
    b_calc[s] <- sum(biomass_function(dbh.temp[nl:nu], 
                                      spp.biomass.params = spp.biomass.params)) * (1 / 833) * 0.48  # changing units to be kgC/m^2
    
    bcorr[s] <- new.state[s] / b_calc[s] #calculate biomass correction
    
    if (length(pft) > 1) {
      stop("error too many pfts assigned")
    }
    
    b_obs[nl:nu] <- biomass_function(dbh.temp[nl:nu], 
                                     spp.biomass.params = spp.biomass.params) * as.numeric(bcorr[s])
    bMax <- 200
    for (j in nl:nu) {
      dbh.temp[j] <- optimize(merit, c(1, bMax), b_obs = b_obs[j], 
                              spp.biomass.params = spp.biomass.params)$minimum
    }
    
    b_calc1[s] <- sum(biomass_function(dbh.temp[nl:nu],
                                       spp.biomass.params = spp.biomass.params)) * (1 / 833) * 0.48
    nl <- nu + 1
  }
  
  dbh <- dbh.temp
  iage <- iage.temp
  nogro <- nogro.temp  # numeric(200)#hack
  
  #nogro[nogro < 1] <- 1
  
  ntrees <- new.ntrees
  
  # print(dbh[1:ntrees[1]])
  
  # translate agb to dbh
  
  #dbh_spp[s] <- optimize(merit, c(0,200))$minimum bcorr = new.state[i,] /
  # agb.pft[,ncol(agb.pft),1] *(bcorr[s]/ntrees[s]) dbh.temp1[j] <- optimize(merit,
  # c(0,200))$minimum
  
  # for(n in 1:nspec){ slta <- spp.params$SLTA[n] sltb <- spp.params$SLTB[n] fwt <-
  # spp.params$FWT[n] frt <- spp.params$FRT[n] if (agb.pft[n,ncol(agb.pft),1]==0 &
  # new.state[i,n]>0){ abg.pft.temp <- sum(distance.matrix[,n]%*%t(agb.pft[n,ncol(agb.pft),1]))
  # ntrees.temp <- sum(distance.matrix[,n]%*%t(t(as.matrix(ntrees)))) dbh.temp <-
  # dbh[sum(ntrees[1:n])-1] for(j in 1:ntrees.temp){ b_obs <-
  # biomass_function(dbh[j],slta=slta,sltb=sltb,fwt=fwt,frt=frt)*bcorr[n] dbh.temp[j] <-
  # optimize(merit, c(0,200),b_obs=b_obs)$minimum } } nu <- nl + ntrees[n] - 1 nl <- nu + 1 }
  
  ##### SOIL
  if ("TotSoilCarb" %in% names(new.state.other)) {
    leaf.sum <- sum(tyl[1:12]) * 0.48
    if(new.state.other["TotSoilCarb"] > 1000) new.state.other["TotSoilCarb"] = rnorm(1,1000,10)
    soil.org.mat <- new.state.other["TotSoilCarb"] - leaf.sum
    soil.corr <- soil.org.mat / (sum(C.mat[C.mat[1:ncohrt, 5], 1]) * 0.48)
    #if(soil.corr > 1) soil.corr <- 1
    C.mat[C.mat[1:ncohrt, 5], 1] <- C.mat[C.mat[1:ncohrt, 5], 1] * as.numeric(soil.corr)
    C.mat[is.na(C.mat[1:ncohrt,1]),1] <- 0
    C.mat[C.mat[1:ncohrt,1] < 0,1] <- 0
  }
  
  if (RENAME) {
    file.rename(file.path(settings$rundir, runid, "linkages.restart.Rdata"), 
                file.path(settings$rundir, runid, paste0(start.time, "linkages.restart.Rdata")))  # save original output
  }
  restart.file <- file.path(settings$rundir, runid, "linkages.restart.Rdata")
  sprintf("%s", restart.file)
  
  
  save(dbh, tyl, ntrees, nogro, ksprt, iage, C.mat, ncohrt, bcorr, file = restart.file)
  
  # make a new settings with the right years min start date and end date - fail in informative way
  
  settings$run$start.date <- paste0(formatC(year(start.time + 1), width = 4, format = "d", flag = "0"), "/01/01")
  settings$run$end.date <- paste0(formatC(year(stop.time), width = 4, format = "d", flag = "0"), "/12/31")
  
  do.call(write.config.LINKAGES, 
          args = list(trait.values = new.params, settings = settings, run.id = runid, 
                      restart = TRUE, spinup = FALSE, inputs = inputs))
  
  # save original output
  if (RENAME) {
    file.rename(file.path(outdir, runid, "linkages.out.Rdata"), 
                file.path(outdir, runid, paste0(start.time, "linkages.out.Rdata")))
  }
} # write_restart.LINKAGES
