library(tidyverse)
library(NMOF)

data <- list(mtx = matrix(c(0.00,  0.00,  0.90,  0.00,  0.10, 0.25,  0.50,  0.00,   
                            0.25,  0.00, 0.60,  0.20,  0.00,  0.20,  0.00), ncol = 3),
             pop = c(2600, 16200, 13400),
             target = c(7190, 4650, 5920))

# True mtx
mtx.true <- matrix(c(0.00,  0.00,  0.75,  0.00,  0.25, 0.35,  0.45,  0.00,
                     0.20,  0.00, 0.65,  0.15,  0.00,  0.20,  0.00), ncol = 3)

# Function to optimise
households <- function(x, data) {

  # Estimate living arrangements
  z <- x %*% data$pop 
  
  # Estimate number of households using parent cohorts (1,2 and 4)
  z <- c(z[1,1]/2, z[2,1]/2, z[4,1]) - data$target
  sum(abs(z))
  
}

# Local search function to perturb propensities
neighbour <- function(x, data) {

  # Choose random column from mtx
  i <- sample(1:ncol(x), 1)
  # Select two non-zero propensities from mtx column
  j <- which(x[, i] != 0) %>% sample(2, replace = FALSE)

  # Randomnly select one to perturb positively 
  x[j[1], i] <- 0.1 * (1 - x[j[1], i]) + x[j[1], i]
  # Perturb second propensity to ensure mtx column adds to 1
  x[j[2], i] <- x[j[2], i] + (1 - sum(x[,i]))

  x

}

# Local search algorithm inputs 
localsearch <- list(x0 = data$mtx,             
                    neighbour = neighbour,
                    nS = 1000,                
                    printBar = TRUE)

# Execute 
now <- Sys.time()
results <- vector("list", 50)
for (i in 1:100) {
  results[[i]] <- LSopt(OF = households, algo = localsearch, data)
}

Sys.time() - now

# Inspect best hit propensity matrices
y <- map_df(results, `[`, c('OFvalue')) %>% with(which.min(OFvalue))
results[[y]]$xbest %>% round(2)
mtx.true
(results[[y]]$xbest - mtx.true) %>% round(2)

# Test

# Estimate living arrangements
z <- mtx.true %*% data$pop  #results[[y]]$xbest %>% round(2) %*% data$pop 

# Estimate number of households using parent cohorts (1,2 and 4)
z <- c(z[1,1]/2, z[2,1]/2, z[4,1]) - data$target
z

