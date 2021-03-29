library(mgcv)
library(glmmTMB)

# Load data:
load("gaspereau 2007-2019.rdata")

aggregate(data$f, by = data["year"], mean)

model <- glmmTMB(f ~ year, family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (year | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (year + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (depth + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(depth + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (length + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(length + 0 | group), family = poisson, data = data, verbose = TRUE)

# model <- glmmTMB(f ~ (station + 0 |group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(year + 0 | group) + exp(position + 0 | group), family = poisson, data = data, verbose = TRUE)

# Year x depth:
model <- glmmTMB(f ~ ar1(year + 0 | group) + ar1(depth + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (depth + 0 | year), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (year + 0 | depth), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(depth + 0 | year), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(year + 0 | depth), family = poisson, data = data, verbose = TRUE)

# Year x length:
model <- glmmTMB(f ~ ar1(year + 0 | group) + ar1(length + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (length + 0 | year), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ (year + 0 | length), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(length + 0 | year), family = poisson, data = data, verbose = TRUE)

# Quadratic length:
model <- glmmTMB(f ~ (len|year) + len2 + ar1(length + 0 | year), family = poisson, data = data, verbose = TRUE)
newdata <- data.frame(year = factor(2014, levels = levels(data$year)),
                      len = 13:33,
                      len2 = (13:33)^2,
                      length = factor(13:33, levels = levels(data$length)))
gbarplot(predict(model, newdata = newdata), 13:33)




model <- glmmTMB(f ~ year + ar1(length + 0 | year), family = poisson, data = data, verbose = TRUE)
newdata <- data.frame(year = factor(2013, levels = levels(data$year)),
                      len = 13:33,
                      len2 = (13:33)^2,
                      length = factor(13:33, levels = levels(data$length)))
gbarplot(exp(predict(model, newdata = newdata)), 13:33)

model <- glmmTMB(f ~ (length + 0 | year), family = poisson, data = data, verbose = TRUE)

image(2007:2019, 13:33, t(as.matrix(ranef(model)[[1]]$length)))

# Length, depth and year:
model <- glmmTMB(f ~ year + ar1(depth + 0 | group) + ar1(length + 0 | year), family = poisson, data = data, verbose = TRUE)
newdata <- data.frame(year = factor(2009, levels = levels(data$year)),
                      depth = sort(unique(data$depth))[6],
                      group = factor(1),
                      length = factor(13:33, levels = levels(data$length)))
gbarplot(exp(predict(model, newdata = newdata)), 13:33)

model <- glmmTMB(f ~ ar1(depth + 0 | group) + ar1(year + 0 | length), family = poisson, data = data, verbose = TRUE)

# Position models:
model <- glmmTMB(f ~ exp(position + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ exp(position + 0 | group), family = poisson, data = data, verbose = TRUE)
model <- glmmTMB(f ~ ar1(depth + 0 | group) + ar1(length + 0 | year) + exp(position + 0 | group), family = poisson, data = data, verbose = TRUE)

grid$group <- factor(1)
grid$year <- factor(2009, levels(data$year))
grid$length <- factor(25, levels(data$length))
mu <- predict(model, newdata = grid)

model <- glmmTMB(f ~ ar1(depth + 0 | group) + ar1(length + 0 | year) + exp(position + 0 | group), family = poisson, data = data, verbose = TRUE)

# This seems to converge:
formula <- f ~ ar1(depth + 0 | group) +      # AR(1) depth effects
               ar1(length + 0 | year) +      # AR(1) length effects
               exp(position + 0 | group) +   # Global spatial effect.
               exp(position + 0 | year)      # Anuual spatial effects.
model2 <- glmmTMB(formula, family = poisson, data = data, verbose = TRUE)

# Next add gear effects, gear x length, maybe, tow distance offsets:
# Station effects may be necessary to do comparative effects:
# Model, local, annual, length-frequencies.

predict(model, newdata = grid)

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


library(mvtnorm)
library(rockchalk)
simGroup <- function(g, n=6, rho=0.7) {
    x <- mvrnorm(mu = rep(0,n),
             Sigma = rho ^ as.matrix(dist(1:n)) )   ## Simulate the process
    y <- x + rnorm(n)                               ## Add measurement noise
    times <- factor(1:n)
    group <- factor(rep(g,n))
    data.frame(y, times, group)
}
dat1 <- do.call("rbind", lapply(1:1000, simGroup) )

(fit.ar1 <- glmmTMB(y ~ ar1(times + 0 | group), data=dat1))

dat1$times <- numFactor(dat1$times)
levels(dat1$times)

fit.mat <- glmmTMB(y ~ mat(times + 0 | group), data=dat1, dispformula=~0)
fit.exp <- glmmTMB(y ~ exp(times + 0 | group), data=dat1)


d <- data.frame(z = as.vector(volcano),
                x = as.vector(row(volcano)),
                y = as.vector(col(volcano)))

set.seed(1)
d$z <- d$z + rnorm(length(volcano), sd=15)
d <- d[sample(nrow(d), 100), ]

volcano.data <- array(NA, dim(volcano))
volcano.data[cbind(d$x, d$y)] <- d$z
image(volcano.data, main="Spatial data", useRaster=TRUE)

d$pos <- numFactor(d$x, d$y)
d$group <- factor(rep(1, nrow(d)))

f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data=d)

x <- sample(1:2, 10, replace=TRUE)
y <- sample(1:2, 10, replace=TRUE)
(pos <- numFactor(x,y))

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
