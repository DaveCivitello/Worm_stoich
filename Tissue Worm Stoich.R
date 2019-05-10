library(deSolve)

Worm_stoich<-function (t, y, parameters){  
  N = y[1]; A = y[2]; QA = y[3]; H = y[4]; QH = y[5]; P = y[6]; QP = y[7]
  with(    
    as.list(parameters),
    {  

      # Assimilation efficiencies
      e_A = (QAx - QA)/(QAx - QAn) # nutrient assimilation efficiency of Autotrophs
      e_H = (QHx - QH)/(QHx - QHn) # nutrient assimilation efficiency of Herbivores
      e_P = (QPx - QP)/(QPx - QPn) # nutrient assimilation efficiency of Parasites
      
      # Probability of infection
      P_inf = ((sigma*f*H/(hA+A))/(dE + f*H/(hA+A)))

      dNdt = s - l*N - v*N/(hN + N)*e_A*A +
        dA*A*QA + dH*H*QH + (dH+alpha+dP)*P*QP  +  alpha*P*QH + alpha*P^2/H*(k+1)/k*QP +
        f*A*QA*H/(hA+A)*(1 - e_H) + f_p*P*H*QH*(1 - e_P)/(hH + H) #+
        #(1 -P_inf)*f_p*P*H*QP/(hH + H)*(1 - QPn/QP)

      dAdt = mu_A*(1-A/K)*(1 - QAn/QA)*A - dA*A - f*A*H/(hA + A)
      dQAdt = v*(N/(hN + N))*e_A - mu_A*(1-A/K)*(1 - QAn/QA)*QA
      
      dHdt = f*A*H/(hA + A)*(1 - QHn/QH) - dH*H - alpha*P - f_p*P*H/(hH + H)
      dQHdt = f*A*QA*e_H/(hA + A) - f*A/(hA + A)*(1 - QHn/QH)*QH 
      
      dPdt = P_inf*f_p*P*H/(hH + H)*(1 - QPn/QP) - (dH+alpha+dP)*P - alpha*P^2/H*(k+1)/k
      dQPdt = f_p*H*QH*e_P/(hH + H) - P_inf*f_p*H/(hH + H)*(1 - QPn/QP)*QP

      res = c(dNdt,dAdt,dQAdt,dHdt,dQHdt,dPdt, dQPdt) 
      list(res)
    } 
  )
}

params = c(# nutrient pars
  s = 0, l=0,
  # autotroph pars
  v = 1, mu_A=0.1, K=100, hN = 1, QAx = 0.1, QAn = 0.01, dA = 0.001, 
  # Host pars
  f =0.3, QHx = 0.15, QHn = 0.1, dH = 0.02, hA = 50,
  # Parasite pars
  f_p =2, hH = 2, QPx = 0.2, QPn = 0.17, dP = 0.01, dE = 0.2, sigma=0.8, 
  alpha = 0.01, k = 10)

inits = c(N = 5, A = 1, QA = 0.02, H = 0.1, QH = 0.15, P = 0.1, QP = 0.172)
run = lsoda(y = inits, times=0:2000, parms = params, func=Worm_stoich, atol=1e-9, rtol=1e-9)
plot(run)

Total_N = run[,"N"] + run[,"A"]*run[,"QA"] + run[,"H"]*run[,"QH"] + run[,"P"]*run[,"QP"]
var(Total_N)
plot(run[,"time"], Total_N, typ="l")


N.out = numeric()
A.out = numeric()
QA.out = numeric()
H.out = numeric()
QH.out = numeric()
P.out = numeric()
QP.out =  numeric()

ns = seq(from=0.0001, to=0.2, length.out = 100)
for (n in 1:200){
  inits["N"] = n/40
  output = lsoda(y = inits, times=0:2000, parms = params, func=Worm_stoich)
  N.out[n] = mean(output[1001:2001, 2])
  A.out[n] = mean(output[1001:2001, 3])
  QA.out[n] = mean(output[1001:2001, 4])
  H.out[n] = mean(output[1001:2001, 5])
  QH.out[n] = mean(output[1001:2001, 6])
  P.out[n] = mean(output[1001:2001, 7])
  QP.out[n] = mean(output[1001:2001, 8])
}
par(mfrow=c(3, 3))
plot(1:200, N.out, typ="l")
plot(1:200, A.out, typ="l")
plot(1:200, QA.out, typ="l")
plot(1:200, H.out, typ="l")
plot(1:200, QH.out, typ="l")
plot(1:200, P.out, typ="l")
plot(1:200, QP.out, typ="l")

run.dt = data.frame( "N" = rep((1:200)/4, times=7), "Output" = c(N.out, A.out, QA.out, H.out, QH.out, P.out, QP.out), "Label" = rep(c("N", "A", "QA", "H", "QH", "P", "QP"), each=200))

N.out = numeric()
A.out = numeric()
QA.out = numeric()
H.out = numeric()
QH.out = numeric()
P.out = numeric()
QP.out =  numeric()

ns = seq(from=0.0001, to=0.2, length.out = 100)
for (n in 1:200){
  inits["N"] = n/40
  params["f_p"] = 0
  output = lsoda(y = inits, times=0:2000, parms = params, func=Worm_stoich)
  N.out[n] = mean(output[1001:2001, 2])
  A.out[n] = mean(output[1001:2001, 3])
  QA.out[n] = mean(output[1001:2001, 4])
  H.out[n] = mean(output[1001:2001, 5])
  QH.out[n] = mean(output[1001:2001, 6])
  P.out[n] = mean(output[1001:2001, 7])
  QP.out[n] = mean(output[1001:2001, 8])
}
par(mfrow=c(3, 3))
plot(1:200, N.out, typ="l")
plot(1:200, A.out, typ="l")
plot(1:200, QA.out, typ="l")
plot(1:200, H.out, typ="l")
plot(1:200, QH.out, typ="l")
plot(1:200, P.out, typ="l")
plot(1:200, QP.out, typ="l")

run.dt2 = data.frame( "N" = rep((1:200)/4, times=7), "Output" = c(N.out, A.out, QA.out, H.out, QH.out, P.out, QP.out), "Label" = rep(c("N", "A", "QA", "H", "QH", "P", "QP"), each=200))

library(ggplot2)

par(mfrow=c(1,2))

fig1 = ggplot(subset(run.dt, Label %in% c("A", "H", "P")),aes(x = N, y = Output, fill = Label)) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "None")+
        scale_y_continuous(limits=c(0,65)) +
        scale_fill_manual(values=c("dark green", "magenta", "brown")) +
        geom_area(position = 'stack')
fig2 = ggplot(subset(run.dt2, Label %in% c("A", "H", "P")), aes(x = N, y = Output, fill = Label)) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "None")+
        scale_y_continuous(limits=c(0,65)) +
        scale_fill_manual(values=c("dark green", "magenta", "brown")) +
        geom_area(position = 'stack')

library(cowplot)
plot_grid(fig1, fig2)

QA_P = as.numeric(unlist(subset(run.dt, Label == "QA", select=Output) - params["QAn"])/(params["QAx"] - params["QAn"]))
plot(1:200, QA_P, typ="l", col="dark green", lwd=3, ylab="Relative nutrient quota")
QH_P = as.numeric(unlist(subset(run.dt, Label == "QH", select=Output) - params["QHn"])/(params["QHx"] - params["QHn"]))
lines(1:200, QH_P, col="brown", lwd=3)
QP_P = as.numeric(unlist(subset(run.dt, Label == "QP", select=Output) - params["QPn"])/(params["QPx"] - params["QPn"]))
lines(1:200, QP_P, col="magenta", lwd=3)

QA = as.numeric(unlist(subset(run.dt2, Label == "QA", select=Output) - params["QAn"])/(params["QAx"] - params["QAn"]))
plot(1:200, QA, typ="l", col="dark green", lwd=3, ylab="Relative nutrient quota")
QH = as.numeric(unlist(subset(run.dt2, Label == "QH", select=Output) - params["QHn"])/(params["QHx"] - params["QHn"]))
lines(1:200, QH, col="brown", lwd=3)
#QP = as.numeric(unlist(subset(run.dt2, Label == "QP", select=Output) - params["QPn"])/(params["QPx"] - params["QPn"]))
#lines(1:200, QP, col="magenta", lwd=3)

QA_diff = QA_P - QA
plot(1:200, QA_diff, typ="l", col="dark green", ylim = c(-1, 1), lwd=3, ylab="Effect of parasitism on nutrient quota")
QH_diff = QH_P - QH
lines(1:200, QH_diff, typ="l", col="brown", lwd=3)
abline(h = 0, col="grey", lty=2, lwd=3)

N_A = subset(run.dt2, Label == "QA", select=Output)*subset(run.dt2, Label == "A", select=Output)
N_H = subset(run.dt2, Label == "QH", select=Output)*subset(run.dt2, Label == "H", select=Output)

N_A_P = as.numeric(unlist(subset(run.dt, Label == "QA", select=Output)*subset(run.dt, Label == "A", select=Output)))
N_H_P = as.numeric(unlist(subset(run.dt, Label == "QH", select=Output)*subset(run.dt, Label == "H", select=Output)))
N_P_P = as.numeric(unlist(subset(run.dt, Label == "QP", select=Output)*subset(run.dt, Label == "P", select=Output)))

run_N = data.frame("N" = rep((1:200)/4, times=3), "Output" = c(N_A_P, N_H_P, N_A_P), "Label" = rep(c("A", "H", "P"), each=200))

fig1 = ggplot(run_N,aes(x = N, y = Output, fill = Label)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "None")+
  scale_y_continuous(limits=c(0,10)) +
  scale_fill_manual(values=c("dark green", "magenta", "brown")) +
  geom_area(position = 'stack')



N_A = as.numeric(unlist(subset(run.dt2, Label == "QA", select=Output)*subset(run.dt2, Label == "A", select=Output)))
N_H = as.numeric(unlist(subset(run.dt2, Label == "QH", select=Output)*subset(run.dt2, Label == "H", select=Output)))
N_P = as.numeric(unlist(subset(run.dt2, Label == "QP", select=Output)*subset(run.dt2, Label == "P", select=Output)))

run_N2 = data.frame("N" = rep((1:200)/4, times=3), "Output" = c(N_A, N_H, N_P), "Label" = rep(c("A", "H", "P"), each=200))

fig2 = ggplot(run_N2,aes(x = N, y = Output, fill = Label)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position = "None")+
  scale_y_continuous(limits=c(0,10)) +
  scale_fill_manual(values=c("dark green", "magenta", "brown")) +
  geom_area(position = 'stack')

plot_grid(fig1, fig2)
