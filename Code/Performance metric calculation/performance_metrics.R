# load packages
library(ggplot2)
library(sp)
library(alphahull)
library(dplyr)

#data import
df0     <- read.csv("mean_noise_data.csv")
d_hull  <- data.frame(mean = log10(df0$mean), cv = log10(df0$cv), row.names = NULL)

#construct hull
alpha   <- 0.212               # choose alpha parameter
hull    <- ahull(d_hull, alpha = alpha) # find alpha hull of point set
indx    <- hull$arcs[, "end1"]          # index of hull vertices in point set
epoints <- d_hull[indx, ]               # match points with vertices from alpha hull object
min.x.ind <- first(which(epoints$mean == min(epoints$mean)))
epoints <- rbind(epoints[min.x.ind:length(epoints$mean),], epoints[1:min.x.ind,]) #begin dataframe with lowest x-value
max.x.ind <- first(which(epoints$mean == max(epoints$mean)))

plot(log10(df0$mean), log10(df0$cv))    # sanity check hull and point set overlay 
lines(epoints$mean, epoints$cv)
text(epoints$mean+0.01, epoints$cv+0.01, c(1:length(epoints$mean)))

#dynamic area
poly <- Polygon(epoints) # convert hull vertices to polygon object
poly@area                # compute area of polygon
F_A <- 10^poly@area #compute F_A
print(F_A)

#F_eta search
upp <- epoints[1:max.x.ind, ] #separate hull vertices in to "upper" and "lower"
low <- epoints[max.x.ind:length(epoints$mean), ] #separate hull vertices into "upper" and "lower" classes

low_lines <- matrix(nrow = length(low$mean)-1, ncol = 4) #matrix for lower hull edges
colnames(low_lines) <- c("m", "b", "beg", "end")
for (i in 1:(length(low$mean)-1)) {
  low_lines[i,1] <- (low$cv[i] - low$cv[i+1]) / (low$mean[i] - low$mean[i+1])
  low_lines[i,2] <- low$cv[i] - low_lines[i,1]*low$mean[i] 
  low_lines[i,3] <- low$mean[i+1]
  low_lines[i,4] <- low$mean[i]
}

upp_lines <- matrix(nrow = length(upp$mean)-1, ncol = 4) #matrix for upper hull edges
colnames(upp_lines) <- c("m", "b", "beg", "end")
for (i in 1:(length(upp$mean)-1)) {
  upp_lines[i,1] <- (upp$cv[i] - upp$cv[i+1]) / (upp$mean[i] - upp$mean[i+1])
  upp_lines[i,2] <- upp$cv[i]- upp_lines[i,1]*upp$mean[i] 
  upp_lines[i,3] <- upp$mean[i]
  upp_lines[i,4] <- upp$mean[i+1]
}

upp_lines <- data.frame(upp_lines) #convert matrix to data.frame
upp_lines <- upp_lines[!is.infinite(rowSums(upp_lines)),] #do not include infinite slopes
low_lines <- data.frame(low_lines) #convert matrix to data.frame
low_lines <- low_lines[!is.infinite(rowSums(low_lines)),] #do not include infinite slopes

lengths <- matrix(ncol = 4, nrow = length(epoints$mean)+1) #create matrix for vertical chord lengths
colnames(lengths) <- c("x", "y", "opposite", "diff")

for (i in 1:length(upp$mean)) { #compute length of chords originating from upper vertices 
  m <- subset(low_lines, beg <= upp$mean[i] & end >= upp$mean[i])$m[1] #choose the first, if multiple match
  b <- subset(low_lines, beg <= upp$mean[i] & end >= upp$mean[i])$b[1] #choose the first, if multiple match
  lengths[i,1] <- upp$mean[i]
  lengths[i,2] <- upp$cv[i]
  lengths[i,3] <- m*upp$mean[i] + b
  lengths[i,4] <- upp$cv[i] - (m*upp$mean[i] + b)
}

for (i in 1:length(low$mean)) { #compute length of chords originating from lower vertices
  m <- subset(upp_lines, beg <= low$mean[i] & end >= low$mean[i])$m[1] #choose the first, if multiple match
  b <- subset(upp_lines, beg <= low$mean[i] & end >= low$mean[i])$b[1] #choose the first, if multiple match
  lengths[(i+length(upp$mean)),1] <- low$mean[i]
  lengths[(i+length(upp$mean)),2] <- low$cv[i]
  lengths[(i+length(upp$mean)),3] <- m*low$mean[i] + b
  lengths[(i+length(upp$mean)),4] <- low$cv[i] - (m*low$mean[i] + b)
}

lengths <- data.frame(lengths)

max.length <- lengths[which(lengths$diff == max(lengths$diff)),] #find longest chord
F_eta <- 10^max.length$diff # F_eta
print(F_eta)

plot(upp, col = "blue", type = "o", ylim = c(-0.6, 0.5)) #visualize hull with F_eta 
points(low, col = "red", type = "o")
segments(lengths$x, lengths$y, lengths$x, lengths$opposite)
segments(max.length$x, max.length$y, max.length$x, max.length$opposite, col = "green", lwd = 3)



