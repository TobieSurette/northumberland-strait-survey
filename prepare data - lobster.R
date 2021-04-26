library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)
library(glmmTMB)

years <- 2007:2020
x <- read.nssset(years, survey = "regular")
y <- read.nsscat(years, survey = "regular", species = "lobster")
z <- read.nsslen(years, survey = "regular", species = "lobster")
b <- read.nssbio(years, survey = "regular", species = "lobster")
y <- y[-which(y$size.class == 2), ]

# Fill-in missing set information:
ix <- which(is.na(x$distance) & (x$experiment != 3))
x$distance[ix] <- (1 / 1.852) * distance(-dmm2deg(x$longitude.start[ix]), dmm2deg(x$latitude.start[ix]), -dmm2deg(x$longitude.end[ix]), dmm2deg(x$latitude.end[ix]), pairwise = FALSE)
x <- x[which((x$experiment != 3) & !is.na(x$distance)), ]

ix <- match(z[key(y)], y[key(y)])
# Calculate sampling ratio:
y$ratio <- y$weight.sampled / y$weight.caught
iz <- which(is.na(z$ratio))
ix <- match(z[iz, key(y)], y[key(y)])
z$ratio[iz] <- y$ratio[ix]

ix <-
z$ratio[ix]

z$ratio[is.na(z$ratio)] <- 1

# Attach information to set card:
import(x, fill = 1, var = "ratio") <- aggregate(z["ratio"], by = z[key(x)], unique)
import(x, fill = 0) <- freq(z, by = key(x))
x$gear[which((x$gear == 14) & (year(x) == 2019))] <- 19
x$gear[which((x$gear == 13) & (year(x) == 2019))] <- 16
x$gear[which((x$gear == 13) & (year(x) == 2020))] <- 16
x$gear.type <- gear(x$gear)
x$longitude <- -dmm2deg(lon(x))
x$latitude  <- dmm2deg(lat(x))
x$depth     <- depth(x$longitude, x$latitude)
x <- x[x$depth > 0, ]
x$station   <- deblank(station(x, method  = "latlong"))
tmp         <- deg2km(x$longitude, x$latitude)
names(tmp)  <- c("xkm", "ykm")
x           <- cbind(x, tmp)
fvars       <- names(x)[gsub("[0-9]", "", names(x)) == ""]

# Prepare interpolation grid:
grid       <- read.gulf.spatial("nss stations")
grid       <- grid[grid$longitude <= -61.90, ]
tmp        <- deg2km(grid$longitude, grid$latitude)
names(tmp) <- c("xkm", "ykm")
grid       <- cbind(grid, tmp)
grid$depth <- depth(grid$longitude, grid$latitude)
grid$station <- deblank(grid$station)
grid$depth    <- round(depth(grid$longitude, grid$latitude)/5)*5 # Round depth by 5m bins.
grid <- grid[grid$depth > 0, ]
grid$depth <- numFactor(grid$depth)

# Determine integer grid coordinates:
dx <- mean(diff(sort(unique(grid$longitude))))
dy <- mean(diff(sort(unique(grid$latitude))))
fun <- function(x, y){
   v <- NA * x
   for (i in 1:length(x)) v[i] <- which.min(abs(x[i] - y))
   return(v)
}
grid$x <- fun(grid$longitude, seq(min(grid$longitude), max(grid$longitude), dx))
grid$y <- fun(grid$latitude, seq(min(grid$latitude), max(grid$latitude), dy))
grid$position <- numFactor(grid$x, grid$y)

ix <- match(x$station, grid$station.number)
x$x <- grid$x[ix]
x$y <- grid$y[ix]
x$position <- grid$position[ix]
x <- x[!is.na(x$x) & !is.na(x$y), ]  # Remove problem grid coordinates.

# Prepare data for analysis:
x           <- x[, setdiff(names(x), as.character(c(0:12, 34:50)))] # Remove small and over-large fish.
fvars       <- names(x)[gsub("[0-9]", "", names(x)) == ""]
data <- data.frame(f         = as.vector(as.matrix(x[fvars])),
                   length    = as.factor(as.vector(repvec(fvars, nrow = nrow(x)))),
                   year      = as.factor(rep(year(x), length(fvars))),
                   distance  = rep(x$distance, length(fvars)),
                   gear      = as.factor(rep(x$gear.type, length(fvars))),
                   xkm       = rep(x$xkm, length(fvars)),
                   ykm       = rep(x$ykm, length(fvars)),
                   x         = rep(x$x, length(fvars)),
                   y         = rep(x$y, length(fvars)),
                   position  = rep(x$position, length(fvars)),
                   station   = as.factor(rep(x$station, length(fvars))),
                   off       = as.numeric(rep(log(x$distance / 0.625), length(fvars))))

data          <- data[!is.na(data$x) & !is.na(data$y), ]       # Problems with grid coordinates.
data$offset   <- log(data$distance / 0.625)                    # Tow distance offset term.
data$group    <- factor(rep(1, nrow(data)))                    # Dummy variable to treat global effects.
data$depth    <- grid$depth[match(data$station, grid$station)] # Get depth from survey grid.
data$depth.meters <- rep(round(x$depth,1), length(fvars))
data$len <- as.numeric(as.character(data$length))
data$len2 <- data$len^2
data$sampling <- factor(paste0(data$year, "-", data$station))

save(data, grid, file = paste0("gaspereau ", min(years), "-", max(years), ".rdata"))
