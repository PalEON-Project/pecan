
### The following script was written to organize only the needed files from each SDA/free run so as to reduce the size of exported files 
### for large runs. In this script, we move all output, input, and restart files for LINKAGES for each year and SDA output. 


## M Kivi 
## October 2021


rm(list=ls())

# Move all needed files to one directory
id = 14000000292
# HARVARD: 14000000279
# ROOSTER: 14000000280
# GOOSE: 14000000288
# NRP: 14000000292
outdir = paste0('/save/workflows/PEcAn_',toString(id),'/')
setwd(outdir)
ids = list.dirs(file.path(outdir,'run'), full.names = FALSE)[-1]
ids = ids[1001:2000]

settings = read.settings('pecan.SDA.xml')
save(settings, file = 'SDAsettings.Rdata')
#sda.years = c(lubridate::year(settings$state.data.assimilation$start.date):lubridate::year(settings$state.data.assimilation$end.date))
sda.years = c(1960:2010)
sda.years = sda.years[-1] # remove first year because the bias is generally so large
nyr = length(sda.years)

# create some folders for saving
dir.create('export')
dir.create('export/ids')
for (cid in ids){
  dir.create(paste0('export/ids/',cid))
}

system(paste0('mv ', outdir,'/SDA/sda.output.Rdata ',outdir,'export/'))
system(paste0('mv ', outdir, '/SDA/outconfig.Rdata ',outdir,'export/'))
system(paste0('mv ', outdir, 'pecan.SDA.xml ',outdir,'export/'))
system(paste0('mv ', outdir, 'SDAsettings.Rdata ',outdir,'export/'))

for (cid in ids){
  
  print(cid)
  
  cdir = paste0('export/ids/',cid)
  
  # get input
  system(paste0('mv ', file.path(outdir,'run',cid,'linkages.input.Rdata'),' ',cdir))
  
  # get data for each year 
  for (j in 1:nyr){
    
    yr = sda.years[j]
    
    if (j == nyr){
      system(paste0('mv ', file.path(outdir,'run',cid,'linkages.restart.Rdata'),' ',cdir))
      system(paste0('mv ', file.path(outdir,'out',cid,'linkages.out.Rdata'),' ',cdir))
    }else{
      file.copy(file.path(outdir,'run',cid,paste0(toString(yr),'-12-31 23:59:59linkages.restart.Rdata')),
                file.path(cdir, paste0('restart_', toString(yr),'.Rdata')))
      try(file.remove(file.path(outdir,'run',cid,paste0(toString(yr),'-12-31 23:59:59linkages.restart.Rdata'))))
      
      file.copy(file.path(outdir,'out',cid,paste0(toString(yr),'-12-31 23:59:59linkages.out.Rdata')),
                file.path(cdir, paste0('out_', toString(yr),'.Rdata')))
      try(file.remove(file.path(outdir,'out',cid,paste0(toString(yr),'-12-31 23:59:59linkages.out.Rdata'))))
    }
  }
}


# let's rename the restart and output files for 2010 across all ids (also get 1960!)
ind = 1

for (i in ids){
  print(ind)
  ind = ind + 1
  
  cdir = paste0('export/ids/',i)
  
  # rename the restart file
  system(paste0('mv ', cdir, '/linkages.restart.Rdata ',cdir, '/restart_2010.Rdata'))
  
  # rename the out file
  system(paste0('mv ', cdir, '/linkages.out.Rdata ',cdir, '/out_2010.Rdata'))
  
  # now copy over 1960 files too!
  file.copy(file.path(outdir,'out',i,'1960-12-31 23:59:59linkages.out.Rdata'),
            file.path(cdir, paste0('out_1960.Rdata')))
  try(file.remove(file.path(outdir,'out',i,'1960-12-31 23:59:59linkages.out.Rdata')))
  
}



######################
### Free model run export
######################

# Move all needed files to one directory
id = 14000000297
outdir = paste0('/save/workflows/PEcAn_',toString(id),'/')
setwd(outdir)
ids = list.dirs(file.path(outdir,'run'), full.names = FALSE)[-1]

settings = read.settings('pecan.CHECKED.xml')
save(settings, file = 'SDAsettings.Rdata')

# create some folders for saving
dir.create('export')
dir.create('export/ids')
for (cid in ids){
  dir.create(paste0('export/ids/',cid))
}

#system(paste0('mv ', outdir,'/SDA/sda.output.Rdata ',outdir,'export/'))
#system(paste0('mv ', outdir, '/SDA/outconfig.Rdata ',outdir,'export/'))
system(paste0('mv ', outdir, 'pecan.CHECKED.xml ',outdir,'export/'))
system(paste0('mv ', outdir, 'SDAsettings.Rdata ',outdir,'export/'))

for (cid in ids){
  
  print(cid)
  
  cdir = paste0('export/ids/',cid)
  
  # get input
  system(paste0('mv ', file.path(outdir,'run',cid,'linkages.input.Rdata'),' ',cdir))
  
  # get output 
  system(paste0('mv ', file.path(outdir,'out',cid,'linkages.out.Rdata'),' ',cdir))
  
}
