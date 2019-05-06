library(deSolve)

Worm_stoich<-function (t, y, parameters){  
  N = y[1]; A = y[2]; QA = y[3]; H = y[4]; QH = y[5]; P = y[6]
  with(    
    as.list(parameters),
    {  
      # Simplifying functions for calculating ingestion, intersection of saturating function and parasite aggregation
      mean_i_p = a_p*(P/H)/(a_h + a_p*P/H) # proportion of ingestate stolen by parasites
      d2_i_p = -2*a_p^2*a_h/(a_h + a_p*P/H)^3
      var_HP = P/H + (P/H)^2/k
      i_p = mean_i_p  + 0.5*d2_i_p*var_HP # Jensen's inequality
      # Assimilation efficiencies
      e_A = (QAx - QA)/(QAx - QAn) # nutrient assimilation efficiency of 
      e_H = (QHx - QH)/(QHx - QHn)
      
      dNdt = s - l*N - v*N/(hN + N)*e_A*A +
             dA*A*QA + dH*H*QH + (dH+alpha+dP)*P*qP +  alpha*P*QH + alpha*P^2/H*(k+1)/k*qP +
             f*A*H*QA*(1 -i_p)*(1 - e_H) + 
             qP*(1 - sigma*f*H/(dE + f*H))*f*A*H*i_p*min(1,QA/qP)
      
      dAdt = mu_A*(1-A/K)*(1 - QAn/QA)*A - dA*A - f*A*H
      dQAdt = v*(N/(hN + N))*e_A - mu_A*(1-A/K)*(1 - QAn/QA)*QA
      
      dHdt = f*A*H*(1 - i_p)*(1 - QHn/QH) - dH*H - alpha*P
      dQHdt = f*A*QA*(1 -i_p)*e_H - f*A*(1 - i_p)*(1 - QHn/QH)*QH 
                            
      dPdt = (sigma*f*H/(dE + f*H))*f*A*H*i_p*min(1,QA/qP) - (dH+alpha+dP)*P - alpha*P^2/H*(k+1)/k
      
      res = c(dNdt,dAdt,dQAdt,dHdt,dQHdt,dPdt) 
      list(res)
    } 
  )
}

params = c(# nutrient pars
          s = 0, l = 0,
          # autotroph pars
          v = 1, mu_A=1, K=100, hN = 1, QAx = 0.05, QAn = 0.01, dA = 0.01, 
          # Host pars
          f =0.1, QHx = 0.1, QHn = 0.05, dH = 0.05,
          # Parasite pars
          a_p = 5, a_h = 1, qP = 0.12, dP = 0.01, dE = 0.1, sigma=0.9, alpha = 0.0001, k = 10)

inits = c(N = 10, A = 1, QA = 0.02, H = 1, QH = 0.1, P = 1)
run = lsoda(y = inits, times=0:1000, parms = params, func=Worm_stoich)
plot(run)

Total_N = run[,"N"] + run[,"A"]*run[,"QA"] + run[,"H"]*run[,"QH"] + run[,"P"]*params["qP"]
var(Total_N)
plot(run[,"time"], Total_N, typ="l")
# 
N.out = numeric()
A.out = numeric()
QA.out = numeric()
H.out = numeric()
QH.out = numeric()
P.out = numeric()

for (n in 1:100){
  inits["N"] = n/10
  output = lsoda(y = inits, times=0:5000, parms = params, func=Worm_stoich)
  N.out[n] = output[5001, 2]
  A.out[n] = output[5001, 3]
  QA.out[n] = output[5001, 4]
  H.out[n] = output[5001, 5]
  QH.out[n] = output[5001, 6]
  P.out[n] = output[5001, 7]
}
par(mfrow=c(2, 3))
plot(1:100, N.out, typ="l")
plot(1:100, A.out, typ="l")
plot(1:100, QA.out, typ="l")
plot(1:100, H.out, typ="l")
plot(1:100, QH.out, typ="l")
plot(1:100, P.out, typ="l")



