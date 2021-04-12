library(gulf.data)
library(gulf.graphics)
library(mgcv)
library(glmmTMB)

# To do:
# - Number of stations on LF plot.
# - Superimpose bubble-plot over interpolated map.
# - Abundance indices.
# -

# Load data:
load("gaspereau 2007-2020.rdata")

data <- data[data$year %in% 2007:2019, ]
data <- data[data$gear != "Northumberland trawl", ]
data$gear <- factor(as.character(data$gear))

years <- as.numeric(as.character(data$year))

# Plot empirical data:
clg()
m <- kronecker(matrix(1:14, ncol = 2), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
n <- NULL
for (i in 1:length(years)){
   ix <- data$year == years[i]
   z <- data[ix, ]
   print(years[i])
   r <- aggregate(z$f, by = z["length"], sum)
   gbarplot(r[, 2], as.numeric(as.character(r[, 1])))
   n[i] <- sum(r[, 2])
}

# Fit model:
model <- glmmTMB(f ~ gear +
                     ar1(year + 0 | group) +
                     ar1(length + 0 | group) +
                     ar1(depth + 0 | group) +
                     ar1(length + 0 | gear) +
                     ar1(length + 0 | year) +
                     exp(position + 0 | group) +
                     exp(position + 0 | year),
                 offset = off, family = poisson, data = data, verbose = TRUE)


# Next add gear effects, gear x length, maybe, tow distance offsets:
# Station effects may be necessary to do comparative effects:
# Model, local, annual, length-frequencies.
years <- as.numeric(as.character(unique(data$year)))
ix <- which(grid$position %in% data$position)
lens <- 13:33
r <- matrix(NA, nrow = length(years), ncol = length(lens))
dimnames(r) <- list(year = years, length = lens)

mu <- array(NA, dim = c(length(ix), length(lens), length(years)))
dimnames(mu) <- list(1:length(ix), length = lens, year = years)
for (i in 2:length(years)){
   print(years[i])
   for (j in 1:length(lens)){
      print(lens[j])
      grid$length <- factor(lens[j], levels = levels(data$length))
      grid$year   <- factor(years[i], levels = levels(data$year))
      grid$group <- factor(1)
      grid$gear <- unique(data$gear)[1]
      grid$off <- 0
      mu[,j,i] <- predict(model, newdata = grid[ix,])
      r[i,j] <- mean(exp(mu[,j,i]))
   }
}

# Display map:
year <- 2019
map <- matrix(NA, nrow = max(grid$x), ncol = max(grid$y))
for (k in 1:length(ix)) map[grid$x[ix[k]], grid$y[ix[k]]] <- mu[k, 20, year-2006]
print(range(exp(map), na.rm = TRUE))
xx <- -64.86194 + ((1:max(grid$x))-1) * dx
yy <- 45.7 + ((1:max(grid$y))-1) * dy
image(xx, yy, exp(map), zlim = c(0, 0.002))
ux <- x$longitude[year(x) == year]
uy <- x$latitude[year(x) == year]
uz <- apply(x[fvars], 1, sum)[year(x) == year]
points(ux, uy, cex = 0.3 * sqrt(uz))
points(ux[uz == 0], uy[uz == 0], pch = "x", lwd = 2, cex = 0.75)
map("coast")

# Length-frequencies:
f <- apply(exp(mu), 2:3, mean)
clg()
m <- kronecker(matrix(1:14, ncol = 2), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:ncol(f)){
   gbarplot(f[, i], lens, grid = TRUE, xaxt = "n", yaxt = "n", ylim = c(0, 2.5))
   text(par("usr")[1] + 0.8 * diff(par("usr"))[1:2],
        par("usr")[3] + 0.8 * diff(par("usr"))[3:4],
        years[i])

   if (i == 1) axis(2)
   if (i %in% 2:(ncol(f)/2)) axis(2)
   if (i %in% round((ncol(f)/2)*(1:2))) axis(1)

   if (i == round(ncol(f)/4)) mtext("Number per tow", 2, 2.5)
   if (i == ncol(f)) mtext("Length (cm)", 1, 3, at = par("usr")[1])
   box()
}


# Index:
clg()
gbarplot(apply(f, 2, sum))

# Bathymetry map for the Northumberland Strait:
map.new(xlim = c(-65, -61.5), ylim = c(45.6, 47.25))
map("bathymetry")
map("coast")
map.axis(1:2)
box()

# Plot spatial effect:
z <- ranef(model)[[1]]$group
z <- as.matrix(z[grep("position", names(z))])[1,]
str <- gsub("position", "", names(z))
str <- gsub("[()]", "", str)
zx <- as.numeric(unlist(lapply(strsplit(str, ","), function(x) x[1])))
zy <- as.numeric(unlist(lapply(strsplit(str, ","), function(x) x[2])))
map <- matrix(NA, nrow = max(grid$x), ncol = max(grid$y))
d <- depth(grid$longitude, grid$latitude)
for (k in 1:length(z)){
   map[zx[k], zy[k]] <- z[k]
}
xx <- -64.86194 + ((1:max(grid$x))-1) * dx
yy <- 45.7 + ((1:max(grid$y))-1) * dy
image(xx, yy, map)
points(x$longitude[year(x) == 2007], x$latitude[year(x) == 2007])
map("coast")

# Plot depth effect:
z <- ranef(model)[[1]]$group
z <- as.matrix(z[grep("depth", names(z))])[1,]
str <- gsub("position", "", names(z))
str <- gsub("[()]", "", str)
zx <- as.numeric(unlist(lapply(strsplit(str, ","), function(x) x[1])))
zy <- as.numeric(unlist(lapply(strsplit(str, ","), function(x) x[2])))
map <- matrix(NA, nrow = max(grid$x), ncol = max(grid$y))
d <- depth(grid$longitude, grid$latitude)
for (k in 1:length(z)){
   map[zx[k], zy[k]] <- z[k]
}
xx <- -64.86194 + ((1:max(grid$x))-1) * dx
yy <- 45.7 + ((1:max(grid$y))-1) * dy
image(xx, yy, map)
points(x$longitude[year(x) == 2007], x$latitude[year(x) == 2007])
map("coast")

dx <- mean(diff(sort(unique(grid$longitude))))
dy <- mean(diff(sort(unique(grid$latitude))))
fun <- function(x, y){
   v <- NA * x
   for (i in 1:length(x)) v[i] <- which.min(abs(x[i] - y))
   return(v)
}
grid$x <- fun(grid$longitude, seq(min(grid$longitude), max(grid$longitude), dx))
grid$y <- fun(grid$latitude, seq(min(grid$latitude), max(grid$latitude), dy))


xx <- -64.86194 + grid$x- * dx
yy <-  45.70000 + grid$y * dy



# Plot empirical data:
clg()
m <- kronecker(matrix(1:14, ncol = 2), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:length(years)){
   gbarplot(r[i, ], lens, grid = TRUE)
}
# Put number of tows, numberof fish, year:



# 2019 gear code 14 to 19
# 2020 code 13 to 16

# 19  nouveau Northumberland trawl Modified Rock-hopper
# - smaller foot gear smaller rollers.
# - mesh of the whole netting is smaller.
# - Possibly wider wing spread
# - Bigger mouth
# - Flatfish trawl
# - Cod-end with small mesh (as was the case with rock-hopper)

# Alewife and blue-back herring were not separated.
# Alosa pseudoharengus
# Alosa aestivalis

# Truncated surveys (smaller survey areas):
# 2019 - Comparative 33 tows with both trawls (within 3 days) (everywhere)
# 2020 - Comparative 50 paired tows (zone 25 only, instead of a wider survey)

# Figure out how to specify simple random effects:
model <- glmmTMB(f ~ length, family = poisson, data = r, verbose = TRUE)
model <- glmmTMB(f ~ (length | group), family = poisson, data = r, verbose = TRUE)
model <- glmmTMB(f ~ (length | year), family = poisson, data = r, verbose = TRUE)
model <- glmmTMB(f ~ (year|group) + (length | year), family = poisson, data = r, verbose = TRUE)

# Depth effect
# Depth*year effect
# length effect
# length * year effect

# Fit various models:
model <- list()
model[[1]] <- gamm(f ~ 1 + offset(off), random = list(year = ~1, length = ~ 1), family = poisson, data = data)
model[[2]] <- gamm(f ~ 1 + s(xkm, ykm) + offset(off), random = list(year = ~1, length = ~ 1), family = poisson, data = data)
model[[3]] <- gamm(f ~ 1 + s(xkm, ykm) + offset(off), random = list(year = ~1, length = ~ 1), family = poisson, data = data)
model[[4]] <- gamm(f ~ 1 + offset(off), random = list(year = ~1, length = ~ 1, station = ~1), family = poisson, data = data)

# Draw map of samples:
map.new(xlim = c(-65.25, -61.75), ylim = c(45.5, 47.25))
map("coastline")
points(x$longitude, x$latitude, pch = 21, bg = "grey")
points(grid$longitude, grid$latitude, pch = 21, bg = "red")
map.axis(1:4)
box()
tmp <- deg2km(x$longitude, x$latitude)
names(tmp) <- c("xkm", "ykm")
x <- cbind(x, tmp)

# Define variables which define the data set:
xlim = c(-64.9, -61.75)
ylim = c(45.65, 47.13)
k <- 200   # Interpolation resolution.
ratio <- distance(c(xlim[1], xlim[1]), ylim)[1,2] / distance(xlim, c(ylim[1], ylim[1]))[1,2]
k[2] <- round(ratio * k)
block.numbers <- c(1:7, 10)
points <- TRUE
biomass <- NULL
area <- NULL
legal <- TRUE

windows()
m <- matrix(0, ncol = 8, nrow = 10)
m[2:5,2:7] <- 1
m[6:9,2:7] <- 2
layout(m)
par(mar = c(0, 0, 0, 0))
for (i in 1:length(years)){

   # Define vector of available years.
   x$depth2 <- round(x$depth/2)*2
   x$zero <- as.numeric((x$weight.caught == 0))
   x$log.depth <- log(x$depth)
   x$weight.caught[x$zero == 0]
   x$log.density <- NA
   x$log.density[x$zero == 0] <- log(x$weight.caught[x$zero == 0])

   # Convert coordinates to kilometers:
   temp <- deg2km(longitude(x), latitude(x))
   center <- c(mean(temp$x), mean(temp$y))
   x$x <- temp$x - center[1]
   x$y <- temp$y - center[2]

   # Fit zero proportion model:
   z <- list()
   z[[1]] <- gam(zero ~ 1, data = x, family = binomial)
   z[[2]] <- gam(zero ~ log.depth, data = x, family = binomial)
   z[[3]] <- gam(zero ~ s(log.depth), data = x, family = binomial)
   z[[4]] <- gam(zero ~ s(log.depth) + s(x, y), data = x, family = binomial)
   temp <- unlist(lapply(z, AIC))
   index <- which(temp == min(temp))[1]
   z <- z[[3]]

   # Fit log-densities:
   m <- list()
   m[[1]] <- gam(log.density ~ 1, data = x[x$zero == 0, ])
   m[[2]] <- gam(log.density ~ log.depth, data = x[x$zero == 0, ])
   m[[3]] <- gam(log.density ~ s(log.depth), data = x[x$zero == 0, ])
   #if (sum(x$zero == 0) > 60){
   m[[4]] <- gam(log.density ~ s(x, y) , data = x[x$zero == 0, ])
   m[[5]] <- gam(log.density ~ log.depth + s(x, y) , data = x[x$zero == 0, ])
   m[[6]] <- gam(log.density ~ s(log.depth) + s(x, y) , data = x[x$zero == 0, ])
   #  }
   temp <- unlist(lapply(m, AIC))
   index <- which(temp == min(temp))[1]
   m <- m[[3]]

   # Define interpolation grid:
   xi <- seq(xlim[1], xlim[2], len = k[1])
   yi <- seq(ylim[1], ylim[2], len = k[2])
   grid <- expand.grid(xi, yi)
   names(grid) <- c("longitude", "latitude")
   grid$depth <- abs(depth(grid$longitude, grid$latitude))
   grid$log.depth <- log(grid$depth)
   grid$block.number <- block.number(grid$longitude, grid$latitude)
   index <- !is.na(grid$block.number) & (grid$block.number %in% block.number)
   temp <- deg2km(grid$longitude, grid$latitude)
   grid$x <- temp[, 1] - center[1]
   grid$y <- temp[, 2] - center[2]
   grid$fishing.zone <- fishing.zone(grid$longitude, grid$latitude, species = 2550)
   grid$pi <- NA
   temp <- predict(z, newdata = grid[index,])
   grid$pi[index] <- 1/(1+exp(temp))
   grid$log.density <- NA
   temp <- predict(m, newdata = grid[index,], se = TRUE)
   grid$log.density[index] <- temp$fit
   grid$density <- NA
   grid$density[index] <- grid$pi[index] * exp(grid$log.density[index] + (summary(m)$scale)/2 + (temp$se^2)/2)
   grid$density.sd <- NA
   grid$density.sd[index] <- (exp(summary(m)$scale + (temp$se^2))-1) *
      exp(2*grid$log.density[index] + (summary(m)$scale) + (temp$se^2))

   # print(paste(sum(grid$density > zlim[2],na.rm = TRUE), "points in", year[i]))
   grid$density[grid$depth < 6] <- NA
   pp <- t(matrix(grid$density, ncol = k, byrow = TRUE))

   uu <- grid$density
   uu[!(grid$fishing.zone %in% c("25"))] <- NA
   sum(grid$fishing.zone %in% c("26A"))
   uu <- t(matrix(uu, ncol = k, byrow = TRUE))
   biomass[1] <- mean(uu, na.rm = TRUE)

   uu <- grid$density
   uu[!(grid$fishing.zone %in% c("26A"))] <- NA
   uu <- t(matrix(uu, ncol = k, byrow = TRUE))
   biomass[2] <- mean(uu, na.rm = TRUE)

   print(c(mean(x$weight.caught[(x$fishing.zone == "25")]), mean(x$weight.caught[(x$fishing.zone == "26A")])))
   print(c(sd(x$weight.caught[(x$fishing.zone == "25")]), sd(x$weight.caught[(x$fishing.zone == "26A")])))
   print(c(length(x$weight.caught[(x$fishing.zone == "25")]), length(x$weight.caught[(x$fishing.zone == "26A")])))
   print(biomass)

   vv <- grid$density
   area[i] <- sum(!is.na(vv) & (vv >= 400)) / sum(!is.na(vv))

   nlevels <- 70

   levels <- pretty(seq(0, 1000, len = 11), n = 15)
   #pp[pp > 1000] <- ]
   #pp[pp < zlim[1]] <- 0
   #pp[is.na(pp) & grid$block.number %in% c(1:6)] <- 0

   colors <- colorRampPalette(c("blue", "yellow", "red"))(length(levels) - 1)
   axis.side = NULL
   if (i == 1) axis.side = c(2,3)
   if (i == 2) axis.side = c(3,4)
   if (i == 3) axis.side = c(2)
   if (i == 4) axis.side = c(4)
   if (i == 5) axis.side = c(1,2)
   if (i == 6) axis.side = c(1,4)

   gulf.map(xlim = xlim, ylim = ylim, axis.side = axis.side)
   #.filled.contour(xi, yi, pp, levels = levels, col = colors)

   text(-62.5, 46.9, year[i], cex = 2)
   if (points){
      points(longitude(x), latitude(x), cex = 0.12*sqrt(x$weight.caught));
      points(longitude(x)[x$weight.caught == 0], latitude(x)[x$weight.caught == 0], pch = "x", lwd = 2, col = "red")
   }
   map.fishing.zones(species = 2550)
   coastline()
   #map.place.names()
   box()
   print(block.number)

   legend("bottomleft",
          legend = c("0", "1", "10", "100", "1000"),
          pch = c(4, 1, 1, 1, 1),
          col = c("red", "black", "black", "black", "black"),
          pt.cex = c(1, 0.12*sqrt(1), 0.12*sqrt(10), 0.12*sqrt(100), 0.12*sqrt(1000)),
          bg = "white")
}
# Draw colour bar:
#plot.new()
#r <- range(levels)
#plot.window(xlim = c(0, 1), ylim = c(r[1]-diff(r)/10, r[2]+diff(r)/10), xaxs = "i", yaxs = "i")
#rect(0.45, levels[-length(levels)], 0.55, levels[-1L], col = colors)
#str <- as.character(levels[-length(levels)])
#index <- setdiff(1:length(str), grep(".", str, fixed = TRUE))
#str[index] <- paste(str[index], "", sep = "")
#str[length(str)] <- paste(str[length(str)], "+", sep = "")
#text(0.4, (levels[-length(levels)] + levels[-1L])/2, str, cex = 0.8)
#text(0.5, levels[1]-0.035*diff(range(levels)), "kg/tow")
#
#print(biomass[order(year)])
#print(area[order(year)])
