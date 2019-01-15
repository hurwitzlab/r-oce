#!/usr/bin/env Rscript

library(oce)

args = commandArgs(trailingOnly=TRUE)
file.name = args[1]
#file.name = "~/work/imicrobe/jesse/x00.input.csv"

d <- read.csv(file.name)

tmp <- split(d, d$latitude)
names(tmp) <- paste0('stn', seq_along(tmp))

#Import the stations as a CTD object
stns <- lapply(tmp, function(x) {
    o <- order(x$depth) # the depths are not always in order
    ctd <- as.ctd(x$salinity[o], x$temperature[o], swPressure(x$depth[o]),
                  longitude=x$longitude, latitude=x$latitude)
    ctd <- oceSetData(ctd, 'depth', x$depth[o])
    ctd <- oceSetData(ctd, 'eASV_Relative_Abundance', x$Relative_Abundance[o])
})

sec <- as.section(stns) #Create a section from the CTD data
sg <- sectionGrid(sec, p=seq(20, 200, 5)) #For better interpolation (though be cautious with these plots!)

basename <- sapply(strsplit(file.name, split='.', fixed=TRUE), function(x) (x[1]))
outname <- paste0(basename, 'eASV-plot-%02d.png')
png(outname, width=7, height=7, res=300, units='in')

abun <- unlist(lapply(stns, function(x) x[['eASV_Relative_Abundance']]))

#plot(sec, showstations=TRUE, ztype="contour", showBottom=FALSE)

###Dot plot with interpolated temperature up top, only colour proportional to the abundance
par(mfrow=c(2, 1))
plot(sg, which='temperature', xtype="latitude", showstations=TRUE, showBottom=FALSE, ztype='image', zcol=oceColorsTemperature, cex=2, pch=16)
cm <- colormap(abun, col=oceColorsViridis)
#drawPalette(colormap=cm)
plot(sec, which='eASV_Relative_Abundance', xtype="latitude", showstations=TRUE, grid=TRUE, showBottom=FALSE, ztype='points', zcol=oceColorsViridis, cex=2, pch=16)

#Change up the sizes
#cex <- max(log10(abun))/log10(abun) #Normalized sizes
cex <- abun/max(abun) #Normalized sizes

###T-S diagram
par(mfrow=c(1, 1))
cm <- colormap(abun, col=oceColorsViridis)
drawPalette(colormap=cm)
plotTS(sec, pch=19, col=cm$zcol, mar=c(3.5, 3.5, 2, 4), cex=3*cex)

###Dot plot with interpolated temperature up top, color and the size proportional to the abundance
par(mfrow=c(2, 1))
plot(sg, which='temperature', xtype="latitude", showstations=TRUE, showBottom=FALSE, ztype='image', zcol=oceColorsTemperature, cex=2, pch=16)
plot(sec, which='eASV_Relative_Abundance', showBottom=FALSE, ztype='points', zcol=oceColorsViridis, cex=0) # setup the plot
points(sec[['distance']], sec[['pressure']], pch=19, cex=3*cex, col=cm$zcol)

###All interpolated (sketchy!)
par(mfrow=c(2, 1))
plot(sg, which='temperature', xtype="latitude", showstations=TRUE, showBottom=FALSE, ztype='image', zcol=oceColorsTemperature, cex=2, pch=16)
plot(sg, which='eASV_Relative_Abundance', xtype="latitude", showstations=TRUE, grid=TRUE, showBottom=FALSE, ztype='image', zcol=oceColorsViridis, cex=2, pch=16)

dev.off()
