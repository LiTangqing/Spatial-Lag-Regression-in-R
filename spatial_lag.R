 ## construct model
N <- 1000
set.seed(1)

## condo Location : 0 to 1, uniform distributed
condoLoc <- runif(N, min=0, max=1)

## assume footage is between 500 - 3000 (inclusive)
footage <- sample.int(3000-500+1, size = N, replace = T) # draw sample from 1 to 2501
footage <- footage + 499 ## now sample is from 500 to 3000

## Number of bedrooms - integer from 1 to 4, inclusive, uniform distributed
no_bedroom <- ceiling(runif(N, min=0, max=4)) ## convert to integer form

## histograms for location, square footage, and bed rooms
hist(condoLoc, main = "Histogram of Condo Location", xlab = "Condo Location")
hist(footage, main = "Histogram of Footage", xlab = "footage")
hist(no_bedroom, breaks=seq(min(no_bedroom) - 0.5, max(no_bedroom) + 0.5, by = 1),
     main = "Histogram of No. of Bedrooms", xlab = "No of Bedrooms")

## simulate price 
e = rnorm(N)
price = 1200 * footage + 10^6 * no_bedroom + e

df <- data.frame(condoLoc,footage,no_bedroom, price)
names(df) <- c("condoLoc","footage","no_bedroom","price")
## Plot the location on the x-axis, and price on the y-axis 
plot(df$condoLoc,df$price, main = "Price vs Location",
  xlab = "Location", ylab = "Price")

## Run a linear regression model to back out the parameters from the simulated data.
reg_model <- lm(df$price ~ df$footage + df$no_bedroom)
summary(reg_model)

## part b
## Now add a penalty parameter to your simulation from part (a) 
## where by the parameter on the distance squared from the CBD is -8 * 10^7
price_b <- 1200 * footage + 10 ^ 6 * no_bedroom - 8 * 10 ^ 7 * (condoLoc - 0.5) ^ 2 + e

plot(df$condoLoc,price_b, main = "Price vs Location",
  xlab = "Condo Location", ylab="Price")

df["price_b"] <- price_b
## Run an OLS linear regression model to back the parameter without spatial correlation
ols <- lm(price_b~df$footage + df$no_bedroom  )
summary(ols)

## implement Moran's I 
## Assumptions:
## 1. for one particular condo,
## neighbors are defined as all the condos except itself
distance_mat <- as.matrix(1 / abs(df$condoLoc - df[1,]$condoLoc))
for (i in 2:N){
  distance_mat <- cbind(distance_mat, as.matrix(abs(df$condoLoc - df[i,]$condoLoc)))
}

# possible weight functions
# 1. {function (d) 1/d}  ## weight is 1/distance
# 2. {function (d) (1/d)^2} ## weight is 1/distance^2
get_weighted_mat <- function(distance_mat, weight_function){
  # transform the distance to weighted matrix according to weight function
  weight_mat <- weight_function(distance_mat) 
  diag(weight_mat) <- 0  # set diag to be 0
  return (weight_mat)
}

## get weight matrix using different weight function
weight_mat_inv_dist <- get_weighted_mat(distance_mat,{function (d) 1/d})
weight_mat_inv_dist_squ <- get_weighted_mat(distance_mat,{function (d) (1/d)^2})

moranI <- function(weight_mat){
  # Row-standardized spatial weights matrix
  weight_mat.sums.rows = apply(weight_mat,1,sum)
  weight_mat = weight_mat/weight_mat.sums.rows
  y_mean = mean(df$price_b)
  weighted_sim = 0
  y_var = 0
  for (i in 1:N){
    for (j in 1:N){
      weighted_sim = weighted_sim + weight_mat[i,j]*(df[i,]$price_b - y_mean)*(df[j,]$price_b - y_mean)
    }
    y_var = y_var +(df[i,]$price_b - y_mean)^2
  }
  weight_sum = N # since we normalised the weight matrix, each row will have a sum of 1;
  nom = N*weighted_sim
  denom = weight_sum * y_var
  moran1 = nom /denom
  return(moran1)
}

moran1_inv_dist = moranI(weight_mat_inv_dist)
moran1_inv_dist_squ = moranI(weight_mat_inv_dist_squ)

#############################################
## part d 
# 1. create spatial weights matrix
# create absolute distance matrix first 
distance_mat = as.matrix(1/abs (df$condoLoc - df[1,]$condoLoc))
for (i in 2:N){
  distance_mat = cbind(distance_mat,as.matrix(abs(df$condoLoc - df[i,]$condoLoc)))
}
distance_mat[distance_mat > 0.1] <- NA 

# possible weight functions
# 1. {function (d) 1/d}  ## weight is 1/distance
# 2. {function (d) (1/d)^2} ## weight is 1/distance^2
get_weighted_mat <- function(distance_mat,weight_function){
  # transform the distance to weighted matrix according to weight function
  weight_mat = weight_function(distance_mat) 
  diag(weight_mat) <- 0  # set diag to be 0
  return (weight_mat)
}

#weight_mat <- get_weighted_mat(distance_mat,{function (d) 1/d})
weight_mat <- get_weighted_mat(distance_mat,{function (d) (1/d)^2})
weight_mat[is.na(weight_mat)] <- 0
# Row-standardized spatial weights matrix
weight_mat.sums.rows <- apply(weight_mat,1,sum)
weight_mat <- weight_mat/weight_mat.sums.rows

# Converting the spatial weights matrix to a weights list object
# and some inspection on lw and neighbors
library(spdep)
weight_mat.lw <- mat2listw(weight_mat, style="W")

weight_mat.nb <- weight_mat.lw$neighbours
# Spatial neighbour sparse representation
weight_mat.df <- listw2sn(weight_mat.lw)

# Spatial Lag Regression 
lag_reg <- lagsarlm(df$price_b ~ df$footage + df$no_bedroom,
  listw = weight_mat.lw, type = "lag", tol.solve = 1.0e-30)

summary(lag_reg)




