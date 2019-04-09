rm(list = ls())
library(raster)
library(rdgal)
x = raster("x", xmn =120, xmx = 140, ymn = 34, ymx =41, res = 1)
str(x)
res(x)
values(x) = 1:ncell(x)
inMemory(x) # ?inMemory ; ?fromDisk ;?readAll
plot(x)
str(x)


# real data (DEM)
# http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
# 90m DEM
setwd("~")
fn <- "srtm_62_05/srtm_62_05.tif"
r1 <- raster(fn)
plot(r1)
str(r1)

# if the res changes, then raster* does not work!
res(r1) = 1
ncell(r1)
plot(r1)

# crop
r2 = crop(r1, extent(126, 127, 37, 38))
plot(r2)

# merge
r3 = crop(r1, extent(127, 128, 37, 38))
m = merge(r2, r3)
plot(m)


f1 <- "srtm_62_05/srtm_62_05.tif"
f2 <- "srtm_62_06/srtm_62_06.tif"
r1 <- raster(f1)
r2 <- raster(f2)
m1 <- mosaic(r1, r2, fun=mean)
plot(m1)

# supp
res(r1)
r1@crs
a = union(extent(r1), extent(r2))
r = raster(res = res(r1), xmn = a@xmin, xmx = a@xmax, 
           ymn = a@ymin, ymx = a@ymax,
           crs = r1@crs)


# mosaic function
r <- raster()
r1 <- crop(r, extent(-10, 11, -10, 11))
r2 <- crop(r, extent(0, 20, 0, 20))
r3 <- crop(r, extent(9, 30, 9, 30))

r1[] <- 1:ncell(r1)
r2[] <- 1:ncell(r2)
r3[] <- 1:ncell(r3)
plot(r1)

m1 <- mosaic(r1, r2, r3, fun=mean)
plot(m1)

# if you have a list of Raster objects, you can use do.call
x <- list(r1, r2, r3)
names(x)[1:2] <- c('x', 'y')
x$fun <- mean
x$na.rm <- TRUE

y <- do.call(mosaic, x)



# resolution
meuse.raster.disaggregate <- disaggregate(meuse.raster, fact=4)
res(meuse.raster.disaggregate)

###
res(m1)
m1.aggregate <- aggregate(m1, fact=10000/90)
plot(m1.aggregate)
# see 'disaggregate(meuse.raster, fact=4)'

