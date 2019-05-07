library(deSolve)

Worm_stoich<-function (t, y, parameters){  
  N = y[1]; A = y[2]; QA = y[3]; H = y[4]; QH = y[5]; P = y[6]
  with(    
    as.list(parameters),
    {  
      #Simplifying functions for calculating ingestion, intersection of saturating function and parasite aggregation
      mean_i_p = a_p*(P/H)/(1 + a_p*P/H) # proportion of ingestate stolen by parasites
      # 2nd, 3rd, and 4th derivatives
      d2 = -2*a_p^2/(1 + a_p*P/H)^3
      d3 = 6*a_p^3/(1 + a_p*P/H)^4
      d4 = -24*a_p^4/(1 + a_p*P/H)^5
      # variance, skew, kurtosis
      var_HP = P/H + (P/H)^2/k
      skew_HP = (1 + 2*(P/H)/k)/sqrt(var_HP)#2*sqrt(var_HP)/(P/H) - 1/sqrt(var_HP)
      kurt_HP = 3 + 1/var_HP  + 6*(var_HP - P/H)/((P/H)^2)
      i_p = mean_i_p  #+ (1/2)*d2*var_HP + (1/6)*d3*skew_HP + (1/24)*d4*kurt_HP # Taylor series for Jensen's inequality
      #print(mean_i_p  + (1/2)*d2*var_HP + (1/6)*d3*skew_HP) #+ (1/24)*d4*kurt_HP)
       #p_burdens = dnbinom(x= 0:1000, mu=P/H, size=k) # brute force dist
       #i_p = sum(p_burdens*a_p*(0:1000)/(1 + a_p*(0:1000))) # Jensen's inequality

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
          v = 1, mu_A=1, K=400, hN = 1, QAx = 0.05, QAn = 0.01, dA = 0.001, 
          # Host pars
          f =0.001, QHx = 0.1, QHn = 0.05, dH = 0.01,
          # Parasite pars
          a_p =3, qP = 0.125, dP = 0.01, dE = 0.1, sigma=0.9, 
          alpha = 0.0001, k = 2)

inits = c(N = 40, A = as.numeric(params["K"]), QA = 0.02, H = 0.1, QH = 0.06, P = 0.01)
run = lsoda(y = inits, times=0:5000, parms = params, func=Worm_stoich)
plot(run)

Total_N = run[,"N"] + run[,"A"]*run[,"QA"] + run[,"H"]*run[,"QH"] + run[,"P"]*params["qP"]
var(Total_N)
plot(run[,"time"], Total_N, typ="l")

N.out = numeric()
A.out = numeric()
QA.out = numeric()
H.out = numeric()
QH.out = numeric()
P.out = numeric()

ns = seq(from=0.0001, to=0.2, length.out = 100)
for (n in 1:100){
  inits["N"] = n*2
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