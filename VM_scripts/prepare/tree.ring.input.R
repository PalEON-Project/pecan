
## Tree_ring_input.R 
# Author: Ann Raiho, Marissa Kivi 
# Last modified: January 2020

# This script has two purposes: 
# (1) Determine the 98% species for the site based on the tree ring data. 
# (2) Process the estimated biomass estimates into species-level mean and covariance
# matrices for running SDA. 

################################################################################
################################################################################
################################################################################

## Part I: Data formatting

# set up working environment
library(tidyverse)
#library(PEcAn.settings)
#library(PEcAn.DB)
#library(PEcAn.data.land)

rm(list=ls())

### ADJUST VARIABLES HERE
site = 'NORTHROUND'

# input RDS data file location and name
input = readRDS(paste0('/Users/marissakivi/Desktop/PalEON/RW_model_results/',site,'/AGB_TAXA_STAN_',site,'_v2.0_082020.RDS'))

# 98% plot save location and name
plot.loc = paste0('/Users/marissakivi/Desktop/PalEON/RW_model_results/',site,'/sda_priority_spp.jpg')

# final data product save location
output = paste0('/Users/marissakivi/Desktop/PalEON/RW_model_results/',site,'/sda.obs.',site,'.Rdata')
###

# Before anything else, let's remove all NA values in order to reduce variable size.
input <- input %>% filter(!is.na(ab))

########### Which models are available for the site? ###########
# First, look at the number of models used. If more than one, which model do you want to use?
unique(input$model)
model.pick = 'Model RW'
# Then, reduce the input data to be just that model.
input = input %>% filter(model == model.pick) %>% select(-model)

########### Which plots do we want to include? ###########
# On rare occasions, we want to only include data for some of the plots
# with available data at a site. For example, North Round Pond has plots 
# with two distinct stand structures and disturbance histories, so we split 
# them up. If you want to do all the plots, you can skip this step.
# Choose which plots you want to include. 
unique(input$plot)
plts = c(3,4)
input <- input %>% filter(plot %in% plts)

################################################################################
################################################################################
################################################################################

## Part II: Determining 98% species

# set up variables for aggregation step
start_date = min(input$year)
end_date = max(input$year)
obs.mean <- obs.mean.tmp <- list()
obs.cov <- obs.cov.tmp <- list()
obs.times <- c(start_date:end_date)
time.type <- 'year'
biomass2carbon <- 0.48
var.names = 'AGB.pft'
variable = 'ab'

# convert all data points from Mg/ha of biomass to kg/m2 of carbon
input$ab <- udunits2::ud.convert(input$ab,'Mg/ha','kg/m^2') * biomass2carbon #* kgm2Mgha 

# check to make sure no NA-taxon entries persist
input <- input[!is.na(input$taxon),]

# find average biomass for each species for each year and iteration across plots 
input <- input %>% 
  group_by(taxon, year, iter) %>%
  summarize(ab = mean(ab, na.rm = TRUE)) %>% 
  ungroup()

# find average biomass for each species for each year across iterations
mean_mat <- input %>% 
  group_by(taxon, year) %>% 
  summarize(ab = mean(ab,na.rm = TRUE)) %>% 
  ungroup()

# determine overall cumulative biomass contribution of each species across all species 
prior_mat <- mean_mat %>% 
  dplyr::group_by(taxon) %>% 
  dplyr::summarize(contr = sum(ab)) %>%
  arrange(desc(contr))
prior_mat$perc = prior_mat$contr/sum(prior_mat$contr,na.rm=T)
prior_mat$cumsum = cumsum(prior_mat$perc)
prior_mat$taxon = factor(prior_mat$taxon, levels = prior_mat$taxon)
pl = ggplot(prior_mat) +
  geom_point(aes(x=taxon, y=cumsum)) + 
  geom_hline(yintercept = 0.98, col = 'red') + 
  labs(title = 'Overall Cumulative Proportion of Biomass by Species', 
       y = 'Cumulative Proportion of Biomass')
pl

######################## STOP ########################

# looking at the plot, record in the following vector the codes of the species whose points fall
# under the red line 
#species = c('QURU','ACRU','TSCA','QUVE','BEAL') # HF 
#species = c('FRAM','FAGR','BEAL','ACSA','TSCA') # NRP2
#species = c('QURU','PIST','QUAL','ACRU','ACSA') # GE
#species = c('QURU','PIST','ACRU','PCRU') #RH
species = c('TSCA','QURU','PIST','BELE') # NRP34

# save plot
ggsave(pl, filename = plot.loc)

################################################################################
################################################################################
################################################################################

## Part III: Reduce and reformat results for running SDA in PEcAn

# First, reduce data frame to only include top 98% taxa. 
input = input %>% filter(taxon %in% species)
unique(input$taxon)

# Now, we need to label the data for each species in a way that PEcAn will recognize so we 
# can successfully match the data to the model results at the species-level.
# If you're unsure what species the code refers to, you can checkout the ReadMe file for the RDS.

# We will label each species with a PFT name that will remain consistent for each species
# throughout your analysis. Traditionally, for LINKAGES, they follow the pattern: 
# Genus.Species_Common.Name. However, as you will learn later, we need to slightly vary this 
# pattern so as to keep your PFTs/priors distinct from the work of others. A few possibilities 
# include Genus.Species_Common.Name.YOURINITIALS or Genus.Species_Common.Name.VERSIONNUMBER. 
# Please note, I already used Genus.Species_Common.Name.2 for my analysis. You will see this later.

# Once you choose your pattern, you will need to fill in the following vector with the chosen PFT
# name for each of your species IN THE SAME ORDER AS THE SPECIES VECTOR ABOVE. 

species # for order 
x = c('Tsuga.Canadensis_Hemlock.2', 
      'Quercus.Rubra_Northern.Red.Oak.2',
      'Pinus.Strobus_White_Pine.2',
      'Betula.Lenta_Black.Birch.2'
      )
  
#      
#      'Acer.Rubrum_Red.Maple.2',
#      'Picea.Rubens_Red.Spruce.2'

#      'Quercus.Velutina_Black.Oak.2',

#     'Acer.Saccharum_Sugar.Maple.2',
#     'Fraxinus.Americana_White.Ash.2',
#     'Quercus.Alba_White.Oak.2',
#)



 

x = paste0(var.names, '.',x)
names(x) = species
x

################################################################################
################################################################################
################################################################################

## Part IV. Finish up processing of species-level biomass mean and covariance matrices

# adjust names of taxa to be with PFT name from BETY 
input$pft.cat <- x[as.character(input$taxon)]

# put biomass estimates in a matrix (iteration X species X year)
iter_mat <- reshape2::acast(input, iter ~ pft.cat ~ year, value.var = 'ab')
cov.test <- apply(iter_mat,3,function(x){cov(x)})

# find average biomass for each species for each year across iterations
mean_mat <- input %>%
  dplyr::select(-taxon) %>%
  group_by(pft.cat, year) %>% 
  summarize(ab = mean(ab,na.rm = TRUE)) %>% 
  ungroup()

# Organize results into SDA usable format. Note: you will get some warnings after running
# this section. Don't worry that's normal. 
for(t in seq_along(obs.times)){
  
  # first mean matrices
  obs.mean.tmp[[t]] <- rep(NA, length(unique(input$pft.cat)))
  names(obs.mean.tmp[[t]]) <- sort(unique(input$pft.cat))
  for(r in seq_along(unique(input$pft.cat))){
    k <- as.numeric(mean_mat[mean_mat$year==obs.times[t] & mean_mat$pft.cat==names(obs.mean.tmp[[t]][r]), variable])
    if(!is.na(k)){
      obs.mean.tmp[[t]][r] <- k
    }
  }
  
  
  # then, covariance matrices
  obs.cov.tmp[[t]] <- matrix(cov.test[,which(colnames(cov.test) %in% obs.times[t])],
                             ncol = sqrt(dim(cov.test)[1]),
                             nrow = sqrt(dim(cov.test)[1]))
  colnames(obs.cov.tmp[[t]]) <- names(iter_mat[1,,t]) 
  
}

obs.mean <- obs.mean.tmp
obs.cov <- obs.cov.tmp

# This data format is necessary for the output to work with the SDA workflow
names(obs.mean) <- paste0(obs.times,'-12-31')
names(obs.cov) <- paste0(obs.times,'-12-31')

# Save the final obs.list 
obs.list <- list(obs.mean = obs.mean, obs.cov = obs.cov)
save(obs.list,file=output)

################################################################################
################################################################################
################################################################################

plot_data = input %>%
  dplyr::select(-taxon) %>%
  group_by(pft.cat, year) %>% 
  summarize(mean = mean(ab,na.rm = TRUE), 
            lower = quantile(ab, 0.025, na.rm = TRUE),
            upper = quantile(ab, 0.975, na.rm = TRUE)) %>% 
  ungroup()
plot_data$species = plyr::mapvalues(plot_data$pft.cat, 
                                    from = c("AGB.pft.Acer.Rubrum_Red.Maple.2",
                                             "AGB.pft.Betula.Lenta_Black.Birch.2",       
                                             "AGB.pft.Pinus.Strobus_White_Pine.2",
                                             "AGB.pft.Quercus.Velutina_Black.Oak.2",
                                             "AGB.pft.Betula.Alleghaniensis_Yellow.Birch.2",
                                             "AGB.pft.Tsuga.Canadensis_Hemlock.2",
                                             "AGB.pft.Quercus.Rubra_Northern.Red.Oak.2"), 
                                    to = c("Red Maple", 
                                           "Black Birch",
                                           "White Pine",
                                           "Black Oak", 
                                           "Yellow Birch", 
                                           "Hemlock",
                                           "Red Oak"))

## Visualize mean AGB estimates for all species at site over time
ggplot(plot_data, aes(x = year, group = species)) +
  geom_line(aes(y = mean, color = species)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = species), alpha = 0.25, show.legend = FALSE) +
  labs(title = 'Harvard Forest : Biomass by PFT', 
       x = 'Year', 
       y = 'Aboveground Biomass (kgC/m2)', 
       col = 'Species') +
  theme(text = element_text(size = 20))


