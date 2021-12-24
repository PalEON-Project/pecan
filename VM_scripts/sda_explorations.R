# Some site-level visualizations

# available light over time 
ggplot(data.full) + 
  geom_point(aes(x = year, y = gap, color = ensemble)) + 
  theme(legend.position = 'none')

dist.melt = data$dist.melt 
dist.melt$species = plyr::mapvalues(dist.melt$species, from = c(1:length(spp)), to = spp)

# DBH skewness of distributions
# strong positive => higher mass to the left 
ggplot(dist.melt) + 
  geom_point(aes(x = year, y = dbh.skew, color = as.factor(species)), alpha = 0.1) + 
  facet_wrap(~species) + 
  geom_hline(yintercept = 0) + 
  labs(color = 'species', y = 'Skewness of DBH Distribution')

# DBH median
ggplot(dist.melt) + 
  geom_point(aes(x = year, y = dbh.median, color = as.factor(species)), alpha = 0.1) + 
  facet_wrap(~species) + 
  geom_hline(yintercept = 0) + 
  labs(color = 'species', y = 'Median DBH (cm)')


# precipitation and temperatures
climate.means = data.full %>% group_by(year) %>% 
  summarize(s.temp = mean(summer.temp), s.precip = mean(summer.precip), 
            w.temp = mean(winter.temp), w.precip = mean(winter.precip))
climate.means$s.temp = scale(climate.means$s.temp)
climate.means$s.precip = scale(climate.means$s.precip)
climate.means$w.temp = scale(climate.means$w.temp)
climate.means$w.precip = scale(climate.means$w.precip)


pl1 = ggplot() + 
  geom_line(data = climate.means, aes(x = year, y = s.temp/2), col = 'salmon') + 
  geom_histogram(data=sing.terms %>% filter(varname1 %in% c('summer.temp')), 
                 aes(x = year, weight = slope), binwidth = 1) + 
  labs(x = 'year', y = 'summer temperature', title = 'Summer Temperature')  

pl2 = ggplot() + 
  geom_line(data = climate.means, aes(x = year, y = s.precip/2), col = 'dodgerblue1') + 
  geom_histogram(data=sing.terms %>% filter(varname1 %in% c('summer.precip')), 
                 aes(x = year, weight = slope), binwidth = 1) + 
  labs(x = 'year', y = 'summer precipitation', title = 'Summer Precipitation')  

pl3 = ggplot() + 
  geom_line(data = climate.means, aes(x = year, y = w.precip/2), col = 'skyblue2') + 
  geom_histogram(data=sing.terms %>% filter(varname1 %in% c('winter.precip')), 
                 aes(x = year, weight = slope), binwidth = 1) + 
  labs(x = 'year', y = 'winter precipitation', title = 'Winter Precipitation')  

pl4 = ggplot() + 
  geom_line(data = climate.means, aes(x = year, y = w.temp/2), col = 'slateblue') + 
  geom_histogram(data=sing.terms %>% filter(varname1 %in% c('winter.temp')), 
                 aes(x = year, weight = slope), binwidth = 1) + 
  labs(x = 'year', y = 'winter temperature', title = 'Winter Temperature')  


grid.arrange(pl1, pl2, pl3, pl4, ncol = 2)


#### quick 

## avln in Ann's runs at HF

spp = c('sugar maple','yellow birch','american beech','white ash','hemlock')
nspec = length(spp)
# output directory
out.path = '/save/workflows/PEcAn_14000000149/out'
run.path = '/save/workflows/PEcAn_14000000149/run'
ens = list.dirs(out.path, full.names = FALSE)[-1]
ens = ens[201:400]

avln.save = matrix(NA, length(ens), 162)
gf.vec.mat = array(NA, dim = c(4,nspec,length(ens),162))
algf.array = array(NA, dim = c(200, nspec, 162, length(ens)))

#fwt.save = matrix(NA, 4, length(ens))
#sltb.save = matrix(NA, 4, length(ens))
#slta.save = matrix(NA, 4, length(ens))
#g.save = matrix(NA, 4, length(ens))

for (i in seq_along(ens)){
  
  e = ens[i]
  files = list.files(file.path(out.path,e))
  files = files[grep(x = files, pattern = 'linkages.out.Rdata')]
  
  #load(file.path(run.path,e,'linkages.input.Rdata'))
  
  #fwt.save[,i] = spp.params$FWT
  #sltb.save[,i] = spp.params$SLTB
  #slta.save[,i] = spp.params$SLTA
  #g.save[,i] = spp.params$G
  
  for (j in seq_along(files)){
    
    load(file.path(out.path,e,files[j]))
    
    if (j == 1){
      avln.save[i,1:101] = avln[,1]
      gf.vec.mat[,,i,1:101] = gf.vec.save[,,,1]
      algf.array[,,1:101,i] = algf.save.keep[,,,1]
    }else{
      avln.save[i,j+101] = avln[,1]
      gf.vec.mat[,,i,j+101] = gf.vec.save[,,1,1]
      algf.array[,,j,i] = algf.save.keep[,,1,1]
    }
  }
}

# visualize 

library(reshape2)
library(ggplot2)

avln.melt = melt(avln.save)
colnames(avln.melt) = c('ensemble','year','avln')

algf.melt = melt(algf.array)
colnames(algf.melt) = c('tree','species','year', 'ensemble', 'algf')
algf.melt = algf.melt %>% filter(!is.na(algf))

smgf.melt = melt(gf.vec.mat[,2,,])
colnames(smgf.melt) = c('species','ensemble','year','smgf')

sngf.melt = melt(gf.vec.mat[,3,,])
colnames(sngf.melt) = c('species','ensemble','year','sngf')

degdgf.melt = melt(gf.vec.mat[,4,,])
colnames(degdgf.melt) = c('species','ensemble','year','degdgf')

gf.melt = left_join(algf.melt, smgf.melt, by = c('ensemble','year','species')) %>%
  left_join(sngf.melt, by = c('ensemble','year','species')) %>%
  left_join(degdgf.melt, by = c('ensemble','year','species'))

gf.melt.2 = gf.melt %>% group_by(tree, species, year, ensemble) %>%
  summarize(lgf.ind = which.min(c(algf, smgf, sngf, degdgf)),
            lgf.value = min(algf, smgf, sngf, degdgf, na.rm = TRUE))

gf.melt.2$lgf.ind = plyr::mapvalues(gf.melt.2$lgf.ind, 
                                  from = c(1:4), 
                                  to = c('algf','smgf','sngf','degdgf'))
gf.melt.2$species = plyr::mapvalues(gf.melt.2$species, 
                                    from = c(1:nspec), 
                                    to = spp)

# avln 

ggplot(avln.melt) + 
  geom_line(aes(x = year, y = avln, group = ensemble, color = ensemble)) + 
  theme(legend.position = 'none') + 
  labs(title = 'My HF Run')

# lowest growth factor
gf.colors = c('algf'='yellow','smgf'='blue','sngf'='brown','degdgf'='green')
ggplot(gf.melt.2) + 
  geom_point(aes(x = year, y = lgf.value, color = lgf.ind)) + 
  facet_wrap(~species) + 
  scale_color_manual(values = gf.colors) + 
  labs(x = 'year', y = 'lowest growth factor value', 
       title = 'North Round 1/2')
  
# fwt 
fwt.melt$species = plyr::mapvalues(fwt.melt$species, from = c(1:5), to = as.vector(spp.params$Spp_Name))
ggplot(fwt.melt) + 
  geom_histogram(aes(x = fwt, fill = as.factor(species)), binwidth = 1) + 
  xlim(0,200) +
  facet_wrap(~species, scales = 'free') + 
  theme(legend.position = 'none')

# SLTA + SLTB 
mean.SLTA = apply(slta.save, 1, mean)
mean.SLTB = apply(sltb.save, 1, mean)

plot(NULL, xlim = c(0,40), ylim = c(0,25),
     xlab = 'DBH', ylab = 'CANOPY DIAMETER')
abline(b = mean.SLTB[1], a = mean.SLTA[1], col = 'red') 
abline(b = mean.SLTB[2], a = mean.SLTA[2], col = 'blue') 
abline(b = mean.SLTB[3], a = mean.SLTA[3], col = 'green') 
abline(b = mean.SLTB[4], a = mean.SLTA[4], col = 'black') 
abline(b = mean.SLTB[5], a = mean.SLTB[5], col = 'orange')





