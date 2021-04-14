source("simulation.R")
candidate = seq(0, 1, 0.1)
nsim = 100

for (i in 1:12){
  power = c()
  for (k in candidate){
    string = paste("setting", i, "(k, nsim)", sep='')
    power = c(power, eval(parse(text=string)))
  }
  png(filename=paste("../figure/setting", i, ".png", sep=''))
  plot(candidate, power)
  dev.off()
}




