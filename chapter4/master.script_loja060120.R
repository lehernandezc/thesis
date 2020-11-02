# read packages
library(raster)
library(rgdal)
library(gdistance)
library(ResistanceGA)

# Set Orion cluster directory
setwd("/shared2/salmosim/luis2/ecuadoriensis/landgen_analysis/ResistanceGA_Loja")
#Orion Directory to write .asc files and results
write.dir <-"./Loja_281119/"

#Read rasters
ra1<- raster("./dem.loja_new.tif")
ra2 <- raster("./land.loja_new.tif")
ra3<- raster("./roads.loja_new.tif")


# Stacks rasters
cat.stack <- stack(ra1, ra2, ra3)

# Change resolution to 250 m and utm projection
res <- c(250,250)
crs <- "+proj=utm +zone=17 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# dem 250 m res
dem.lowres <- projectRaster(cat.stack[[1]], res = res, crs = crs)

# land 250 m res
land.lowres <- projectRaster(cat.stack[[2]], crs = crs, res = res, method = "ngb")

# roads 250 m res
roads.lowres <- projectRaster(cat.stack[[3]], crs = crs, res = res, method = "ngb")

# stack new rasters
lowres.stack <- stack(dem.lowres, land.lowres, roads.lowres)

# Reclassify categorical and feature rasters based on intial hypothesis: managed and cropland have lowest resistance
  # Original                                                          #reorder #Final scale during optimisation
#     1     Tree Cover, broadleaved, evergreen                        (8)       70               
#     2     Tree Cover, broadleaved, deciduous, closed                (7)       60
#     7     Tree Cover, regularly flooded, fresh  water (& brackish)  (10)      90
#     12    Shrub Cover, closed-open, deciduous                      (4)        30
#     13    Herbaceous Cover, closed-open                            (5)        40
#     14    Sparse Herbaceous or sparse Shrub Cover                  (6)        50
#     15    Regularly flooded Shrub and/or Herbaceous Cover          (9)        80
#     16    Cultivated and managed areas                             (1)        1
#     17    Mosaic: Cropland / Tree Cover / Other natural vegetation (3)        20
#     18    Mosaic: Cropland / Shrub or Grass Cover                  (2)        10
#     20    Water Bodies (natural & artificial)                      (11)       100

land.reclass <- reclassify(lowres.stack[[2]], c(0, 1, 8,
                                                1.1, 2, 7,
                                                2.1, 7, 10,
                                                7.1, 12, 4,
                                                12.1, 13, 5,
                                                13.1, 14, 6,
                                                14.1, 15, 9,
                                                15.1, 16, 1,
                                                16.1, 17, 3,
                                                17.1, 18, 2,
                                                18.1, 20, 11))


# reclassify roads based on intial hypothesis: Highways and trails are the lowest resistances.
# Original values    #reorder #final scale for optimisation
# 0=Trails,           (2)       25
# 1=Highways,         (1)       1
# 2=secondary roads,  (3)       50
# 3=third roads       (4)       75
# 4=NA                (5)       100
#  1, 25, 50, 75, 100
roads.reclass <- reclassify(lowres.stack[[3]], c(-Inf, 0, 2,
                                                 0.1, 1, 1,
                                                 1.1, 2, 3,
                                                 2.1, 3, 4,
                                                 3.1,4, 5))


# Crop rasters extent to speed up upcoming analyses

ext <- extent(582149,688165,9500533,9560456) # set the new extent

dem.crop <- crop(lowres.stack[[1]], ext)
land.crop <- resample(land.reclass, dem.crop, method="ngb") # method "ngb" for categorical data
roads.crop <- resample(roads.reclass, dem.crop, method="ngb")

# New stack surfaces to optimise
loja.stack <- stack(dem.crop, land.crop, roads.crop)

# Read in points
XY <- read.csv("./UTMpops.csv",header=TRUE)
# Create a site object that is a spatial points object with XY coordinates
sites.loja <- SpatialPoints(XY[ ,c(1,2)])

# Read in genetic distance
genmat <- read.csv("./Pairwise_Gst_pops.csv", header = TRUE, row.names = 1)
genmat <- as.matrix(genmat)
Dgen <- as.dist(genmat)

#####################################################
#  For multisurface optimisation only

######################################################

# Prepapre GA inputs
# Maximum value to explore is 700 and minimum is default.
# Run in 20 cores
# Change loja.stack[[-3]] for elevation and land cover; loja.stack[[-2]] for elevation and roads; loja.stack[[-1]] for
# land cover and elevation
GA.inputs <- GA.prep(ASCII.dir = loja.stack,
                     Results.dir = write.dir,
                     method = "LL",
                     max.cat = 700,
                     max.cont = 700,
                     parallel = 20)

# Prepare inputs for commute distance
# response is the genetic distance matrix
gdist.inputs <- gdist.prep(length(sites.loja),
                           response = as.vector(Dgen),
                           samples = sites.loja,
                           method = 'commuteDistance') # Optimize using commute distance

# Select an transformation for the continous elevation variable. This is just an initial hypothesis.
#plot.t <- Plot.trans(PARM = c(5, 100),
#                     Resistance = loja.stack$dem.loja_new,
#                    transformation = "Monomolecular")


# Create a vector with the parameters for the resistance surfaces in the same order as in the loja.stack object.
# Delete parameters of the surface not include on each run.

PARM <- c( 5, 3, 100, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 ,1, 25, 50, 75, 100)


# Combine resistance surfaces
Resist <- Combine_Surfaces(PARM = PARM,
                           gdist.inputs = gdist.inputs,
                           GA.inputs = GA.inputs,
                           out = NULL,
                           rescale = TRUE,
                           p.contribution = TRUE)

# Create the true resistance/response surface
gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                r = Resist$combined.surface)

# run gdsit.prep to add response
gdist.inputs <- gdist.prep(n.Pops = length(sites.loja),
                           samples = sites.loja,
                           response = as.vector(gdist.response),
                           method = 'commuteDistance')

# Run the multisurface function for multisurface optimisation
Multi.Surface_optim <- MS_optim(gdist.inputs = gdist.inputs,
                                GA.inputs = GA.inputs)

##################################################################################################################
# For single surface optimisation only

######

############## ELEVATION single surface optimisation:

GA.inputs <- GA.prep(ASCII.dir = loja.stack[[1]],
                     Results.dir = write.dir,
                     method = "LL",
                     max.cat = 700,
                     max.cont = 700,
                     parallel = 20)

gdist.inputs <- gdist.prep(length(sites.loja),
                           response = as.vector(Dgen),
                           samples = sites.loja,
                           method = 'commuteDistance') # Optimize using commute distance


#Run optimization
elevation_optim <- SS_optim(gdist.inputs = gdist.inputs,
                             GA.inputs = GA.inputs)


########## LAND COVER single surface optimisation:

land.opt.reclas <- reclassify(loja.stack[[2]], c(0,1,1,
                                                 1.1, 2, 10,
                                                 2.1, 3, 20,
                                                 3.1, 4, 30,
                                                 4.1, 5, 40,
                                                 5.1, 6, 50,
                                                 6.1, 7, 60,
                                                 7.1, 8, 70,
                                                 8.1, 9, 80,
                                                 9.1, 10, 90,
                                                 10.1, 11, 100))


GA.inputs <- GA.prep(ASCII.dir = land.opt.reclas,
                     Results.dir = write.dir,
                     method = "LL",
                     max.cat = 700,
                     max.cont = 700,
                     parallel = 20)

gdist.inputs <- gdist.prep(length(sites.loja),
                           response = as.vector(Dgen),
                           samples = sites.loja,
                           method = 'commuteDistance') # Optimize using commute distance


#Run optimization
landcover_optim <- SS_optim(gdist.inputs = gdist.inputs,
                                 GA.inputs = GA.inputs)

################# ROADS single surface optimisation

roads.opt.reclas <- reclassify(loja.stack[[3]], c(0,1,1,
                                                  1.1, 2, 25,
                                                  2.1, 3, 50,
                                                  3.1, 4, 75,
                                                  4.1, 5, 100))


GA.inputs <- GA.prep(ASCII.dir = roads.opt.reclas,
                     Results.dir = write.dir,
                     method = "LL",
                     max.cat = 700,
                     max.cont = 700,
                     parallel = 20)

gdist.inputs <- gdist.prep(length(sites.loja),
                           response = as.vector(Dgen),
                           samples = sites.loja,
                           method = 'commuteDistance') # Optimize using commute distance


#Run optimization
roads_optim <- SS_optim(gdist.inputs = gdist.inputs,
                                 GA.inputs = GA.inputs)


##########################

#### Bootstrap analysis



distance <- read_csv("./distance.csv", col_names = FALSE)
distance <- as.matrix(distance)

elevation <- read_csv("./elevation.csv", col_names = FALSE)
elevation <- as.matrix(elevation)

elevation_land <- read_csv("./elevation_land.csv", col_names = FALSE)
elevation_land <- as.matrix(elevation_land)

elevation_land_roads <- read_csv("./elevation_land_roads.csv", col_names = FALSE)
elevation_land_roads <- as.matrix(elevation_land_roads)

elevation_roads <- read_csv("./elevation_roads.csv", col_names = FALSE)
elevation_roads <- as.matrix(elevation_roads)

land <- read_csv("./land.csv", col_names = FALSE)
land <- as.matrix(land)

land_roads <- read_csv("./land_roads.csv", col_names = FALSE)
land_roads <- as.matrix(land_roads)

roads <- read_csv("./roads.csv", col_names = FALSE)
roads <- as.matrix(roads)

distance_k <- as.data.frame(read_csv("./distance_k.csv", col_names = TRUE))
elevation_k <- as.data.frame(read_csv("./elevation_k.csv", col_names = TRUE))
elevation_land_k <- as.data.frame(read_csv("./elevation_land_k.csv", col_names = TRUE))
elevation_land_roads_k <- as.data.frame(read_csv("./elevation_land_roads_k.csv", col_names = TRUE))
elevation_roads_k <- as.data.frame(read_csv("./elevation_roads_k.csv", col_names = TRUE))
land_k <- as.data.frame(read_csv("./land_k.csv", col_names = TRUE))
land_roads_k <- as.data.frame(read_csv("./land_roads_k.csv", col_names = TRUE))
roads_k <- as.data.frame(read_csv("./roads_k.csv", col_names = TRUE))



#Make a list of cost/resistance distance matrices
mat.list <- list(distance, elevation, elevation_land, elevation_land_roads, elevation_roads, land, land_roads, roads)
names(mat.list) <- c("distance", "elevation", "elevation_land", "elevation_land_roads", "elevation_roads", "land", "land_roads", "roads")

k <- rbind(distance_k, elevation_k, elevation_land_k, elevation_land_roads_k, elevation_roads_k, land_k, land_roads_k, roads_k)


genmat <- read.csv("./Pairwise_Gst_pops.csv", header = TRUE, row.names = 1)
genmat <- as.matrix(genmat)

Dgen <- as.dist(genmat)


# Create square distance matrix for response for use with
# the bootstrap function
response <- matrix(0, 25, 25)
response[lower.tri(response)] <- as.vector(Dgen)


# Run bootstrap
(AIC.boot <- Resist.boot(mod.names = names(mat.list),
                         dist.mat = mat.list,
                         n.parameters = k[,2],
                         sample.prop = 0.75,
                         iters = 10000,
                         obs = 25,
                         genetic.mat = response
))


write.csv(AIC.boot, file = "AIC.Boot.csv")




setwd("C:/Users/quiqu/Dropbox/Rhodnius_project_2019/LandGen_analysis/R_scripts/ResistanceGA_Loja/all_comb/AIC.boot")

# read saved optimisation files

elevation <- readRDS("./dem.loja_new.rds")
land <- readRDS("./land.loja_new.rds")
roads <- readRDS("./roads.loja_new.rds")
elr <- readRDS("./dem.loja_new.land.loja_new.roads.loja_new.rds")
er <- readRDS("./dem.loja_new.roads.loja_new.rds")
el <- readRDS("./dem.loja_new.land.loja_new.rds")
lr <- readRDS("./land.loja_new.roads.loja_new.rds")









