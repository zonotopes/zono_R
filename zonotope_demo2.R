# zonotope_demo2.R
# test more datasets

source("zonotope.R")

star <- read.table("zono_data/R5gen30.txt", header=F, skip=3)
#star <- read.table("zono_data/R5gen50_Le.txt", header=F, skip=3)
#star <- read.table("zono_data/R5gen100_Le.txt", header=F, skip=3)
star <- as.matrix(star)
v = zonotope(star); # Elapsed time: 8 min 02 sec
