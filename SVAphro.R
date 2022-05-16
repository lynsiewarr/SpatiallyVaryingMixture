# Run script for APHRODITE

source('SVMixFixed_SameComp.R')

system.time(SVMix("aphro",ncomp=8,vconst=1,keep.its=1000,burn.its=50000,thin.its=200,plot.fits=TRUE,plot.trace=TRUE,num.cores=15L))
