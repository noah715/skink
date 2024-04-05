library(boot)

fecundity.fun <- \(N, H, c1 = 1, c2 = 1, c3 = 1){
  f <- log(c2*N) * (1 - N / (c1 * H)) + c3

  p <- inv.logit(f)
  
  return(p)
}


##### Finding average mortality ######

mu = c(.65, .6, .35, .1)

muZ.red = 20
muO.red = 20
muT.red = 10
muM.red = 20
avg.eggs = 2.10

muZ = (1 - muZ.red / 100) * mu[1]
muO = (1 - muO.red / 100) * mu[2]
muT = (1 - muT.red / 100) * mu[3]
muM = (1 - muM.red / 100) * mu[4]
average_eggs = avg.eggs


leslie <- matrix(c(   0,       0,       0, average_eggs,
                      1 - muZ,       0,       0,            0,
                      0, 1 - muO,       0,            0,
                      0,       0, 1 - muT,      1 - muM),
                 nrow = 4, ncol = 4, byrow = T)

leslie.eig <- eigen(leslie)

mean.mu <- Re(weighted.mean(mu, leslie.eig$vectors[,1]))

##### Setting up system to solve #####

eta = 4
n = 2500
carrying.cap = 12000
patch_areas = c(130, 80, 20, 140)
c1 = carrying.cap / sum(patch_areas)

phi.bar = average_eggs / 4

n = 2500
h = sum(patch_areas)

c3 <- logit(mean.mu)
c2 <- exp((logit(phi.bar) - c3) / (1 - n/carrying.cap)) / 2500

c.vals <- c(c1, c2, c3)

test_vals <- seq(0, 30000, length = 200)
test_p <- sapply(test_vals, function(x) pmax(fecundity.fun(x, patch_areas, c1, c2, c3), 0))[1,]

df <- data.frame(Abundance = rep(test_vals, 2), Fecundity = c(test_p, pmax(pmin(((mean.mu - 1)/carrying.cap)*(test_vals) + 1, 2.17/4), 0)),
                 Allee = factor(rep(c("Allee effects", "No allee effects"), each = length(test_vals)), levels = c("Allee effects", "No allee effects")))





