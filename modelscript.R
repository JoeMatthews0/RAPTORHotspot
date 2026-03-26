init = read.table("halle_nodes.txt", header=TRUE)
head(init)

write.csv(init, file="halle_nodes.csv", sep=",", row.names=FALSE, col.names=TRUE)

colls = init$Accidents
year = init$Year

library(MASS)
library(dplyr)

mod = glm.nb(Accidents ~ LogVolume + Year + Urban + Intersection, data = init)

site = 102

mod_data = init |> 
  select(ID, Accidents, Year) |> 
  mutate(mu = fitted(mod)) |> 
  filter(ID == site)

data = list(y = mod_data$Accidents,
            year = mod_data$Year - max(mod_data$Year),
            mu = mod_data$mu,
            predmu = mod_data$mu[nrow(mod_data)] * exp(mod$coefficients["Year"]),
            n_past = nrow(mod_data) - 1)

modelstring = "
  model{
    for(i in 1 : n_past){
      y[i] ~ dnegbin(1 / c[i], lambda[i] / (c[i] - 1))
      lambda[i] <- exp(log(mu[i]) + sigma + (alpha * year[i]))
      c[i] <- exp(-year[i] * tau)
    }
    
    y[n_past + 1] ~ dpois(lambda[n_past + 1])
    lambda[n_past + 1] <- exp(log(mu[n_past + 1]) + sigma)
      
    pred ~ dnegbin(1 / exp(tau), lambda_pred / exp(tau))
    lambda_pred = exp(log(predmu) + sigma + alpha)
      
    sigma ~ dnorm(0, 0.1)
    alpha <- alpha_n * alpha_z
    alpha_n ~ dnorm(0, 1)
    alpha_z ~ dbern(p)
    p ~ dunif(0, 1)
    tau ~ dgamma(2, 20)
  }
"
library(rjags)
model = jags.model(textConnection(modelstring), data = data, n.adapt = 1000)
update(model, n.iter = 1000)
output = coda.samples(model = model,
                      variable.names = c("lambda", "pred"),
                      n.iter = 50000,
                      thin = 5)
