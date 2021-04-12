library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)
library(glmmTMB)

years <- c(2007:2009, 2012:2019)
x <- read.nssset(years, survey = "regular")
z <- read.nsslen(years, survey = "regular", species = "gaspereau")

# Fill-in missing set information:
x$longitude <- -dmm2deg(lon(x))
x$latitude  <- dmm2deg(lat(x))
ix <- which(is.na(x$distance) & (x$experiment != 3))
x$distance[ix] <- (1 / 1.852) * distance(-dmm2deg(x$longitude.start[ix]), dmm2deg(x$latitude.start[ix]), -dmm2deg(x$longitude.end[ix]), dmm2deg(x$latitude.end[ix]), pairwise = FALSE)
x <- x[which((x$experiment != 3) & !is.na(x$distance)), ]
z$ratio[is.na(z$ratio)] <- 1

# Attach information to set card:
import(x, fill = 1, var = "ratio") <- aggregate(z["ratio"], by = z[key(x)], unique)
import(x, fill = 0) <- freq(z, by = key(x))
x$gear[which((x$gear == 14) & (year(x) == 2019))] <- 19
x$gear[which((x$gear == 13) & (year(x) == 2019))] <- 16
x$gear[which((x$gear == 13) & (year(x) == 2020))] <- 16
x <- x[x$gear == 16, ]

# Prepare data for analysis:
x           <- x[, setdiff(names(x), as.character(c(0:12, 34:50)))] # Remove small and over-large fish.
fvars       <- names(x)[gsub("[0-9]", "", names(x)) == ""]

# Standardize catches:
x[fvars] <- repvec((1 / x$ratio) * (0.625 / x$distance), ncol = length(fvars)) * x[fvars]


# Calculate abundance index:
f <- aggregate(apply(x[fvars], 1, sum), by = list(year = year(x)), mean)
f <- f[, -1]
names(f) <- years
n <- aggregate(apply(x[fvars], 1, sum), by = list(year = year(x)), length)
n <- n[, -1]
names(n) <- years
s <- aggregate(apply(x[fvars], 1, sum), by = list(year = year(x)), sd)
s <- s[, -1]
names(s) <- years

clg()
gbarplot(f, grid = TRUE, ylim = c(0, 20))
error.bar(years, lower = f - 1.96 * s / sqrt(n), upper = f + 1.96 * s / sqrt(n))
mtext("Year", 1, 3, cex = 1.5)
mtext("Number per tow", 2, 2.5, cex = 1.5)
box()


# Calculate number per tow length-frequencies:
f <- aggregate(x[fvars], by = list(year = year(x)), mean)
f <- f[, -1]
rownames(f) <- years
n <- aggregate(x[fvars], by = list(year = year(x)), length)
n <- n[, -1]
rownames(n) <- years
s <- aggregate(x[fvars], by = list(year = year(x)), sd)
s <- s[, -1]
rownames(s) <- years

# Length-frequencies:
clg()
m <- kronecker(matrix(1:12, ncol = 2), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:nrow(f)){
   gbarplot(f[i, ], lens, grid = TRUE, xaxt = "n", yaxt = "n", ylim = c(0, 4))
   text(par("usr")[1] + 0.8 * diff(par("usr"))[1:2],
        par("usr")[3] + 0.8 * diff(par("usr"))[3:4],
        years[i], cex = 1.5)

   text(par("usr")[1] + 0.15 * diff(par("usr"))[1:2],
        par("usr")[3] + 0.8 * diff(par("usr"))[3:4],
        paste0(n[i,1], " tows"), cex = 1.25)

   if (i == 1) axis(2)
   if (i %in% 2:6) axis(2, at = 0:3)
   if (i %in% round((nrow(f)/2)*(1:2))) axis(1)

   if (i == round(nrow(f)/4)) mtext("Number per tow", 2, 2.5, at = 0, cex = 1.25)
   if (i %in% c(6, nrow(f))) mtext("Length (cm)", 1, 3, cex = 1.25)
   box()
}



# Plot maps:
clg()
m <- kronecker(matrix(1:12, ncol = 3), matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:length(years)){
   y <- x[year(x) == years[i], ]

   ny <- apply(y[fvars], 1, sum)
   iy <- which(ny > 0)

   map.new(xlim = c(-65, -61.75), ylim = c(45.6, 47.25))
   grid()
   map("coast")
   points(y$longitude[iy], y$latitude[iy], pch = 21, cex = 0.15 * sqrt(ny[iy]), bg = "tomato1")
   points(y$longitude[-iy], y$latitude[-iy], pch = 4, lwd = 1, cex = 0.65)

   text(par("usr")[1] + 0.8 * diff(par("usr"))[1:2],
        par("usr")[3] + 0.8 * diff(par("usr"))[3:4],
        years[i], cex = 1.5)

   if (i %in% 1:4) map.axis(2)
   if (i %in% c(4, 8, 11)) map.axis(1)

   box()
}
plot.new()
v <- c(1, 5, 10, 25, 100, 200)
legend("center",
       legend = c(0, v),
       pt.cex = c(1, 0.15 * sqrt(v)),
       pch = c(4, rep(21, length(v))),
       pt.bg = c("black", rep("tomato1", length(v))),
       title = "# / tow",
      )


