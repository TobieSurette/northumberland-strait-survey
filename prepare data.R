library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)

years <- 2007:2019
x <- read.nssset(years, survey = "regular")
y <- read.nsscat(years, survey = "regular", species = "gaspereau")
z <- read.nsslen(years, survey = "regular", species = "gaspereau")
b <- read.nssbio(years, survey = "regular", species = "gaspereau")

# Attach information to set card:
x <- x[which((x$experiment != 3) & !is.na(x$distance)), ]
import(x, fill = 1, var = "ratio") <- aggregate(z["ratio"], by = z[key(x)], unique)
import(x, fill = 0) <- freq(z, by = key(x))
x$gear.type <- gear(x$gear)
x$longitude <- -dmm2deg(lon(x))
x$latitude  <- dmm2deg(lat(x))
x$depth     <- depth(x$longitude, x$latitude)
x <- x[x$depth > 0, ]
x$station   <- station(x, method  = "latlong")
tmp        <- deg2km(x$longitude, x$latitude)
names(tmp) <- c("xkm", "ykm")
x      <- cbind(x, tmp)
fvars       <- names(x)[gsub("[0-9]", "", names(x)) == ""]

# Prepare interpolation grid:
grid       <- read.gulf.spatial("nss stations")
grid       <- grid[grid$longitude <= -61.90, ]
tmp        <- deg2km(grid$longitude, grid$latitude)
names(tmp) <- c("xkm", "ykm")
grid       <- cbind(grid, tmp)
grid$depth <- depth(grid$longitude, grid$latitude)

# Prepare data for analysis:
x           <- x[, setdiff(names(x), as.character(c(0:12, 34:50)))] # Remove small and over-large fish.
fvars       <- names(x)[gsub("[0-9]", "", names(x)) == ""]
data <- data.frame(f         = as.vector(as.matrix(x[fvars])),
                   length    = as.factor(as.numeric(repvec(fvars, nrow = nrow(x)))),
                   distance  = rep(x$distance, each = length(fvars)),
                   year      = as.factor(rep(year(x), each = length(fvars))),
                   gear      = as.factor(rep(x$gear.type, each = length(fvars))),
                   xkm       = rep(x$xkm, each = length(fvars)),
                   ykm       = rep(x$ykm, each = length(fvars)),
                   depth     = rep(x$depth, each = length(fvars)),
                   log.depth = rep(log(x$depth), each = length(fvars)),
                   station   = as.factor(rep(x$station, each = length(fvars))),
                   off       = as.numeric(rep(log(x$distance / 0.625), each = length(fvars))))

save(data, file = "gaspereau 2007-2019.rdata")

