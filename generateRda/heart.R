setwd('/Users/shijiazhu/Dropbox/mypaper/GIGSEA/codes/data')
heartDisease.MetaXcan <- read.table('heart_metaXcan.csv',sep=',',header=T)

setwd('~/Documents/MyPackages/RStudio/GIGSEA/')
save(heartDisease.MetaXcan,file='data/heartDisease.MetaXcan.rda')
