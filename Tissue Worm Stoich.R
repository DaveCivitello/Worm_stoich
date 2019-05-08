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
  v = 1, mu_A=0.1, K=100, hN = 1, QAx = 0.05, QAn = 0.01, dA = 0.001, 
  # Host pars
  f =0.3, QHx = 0.2, QHn = 0.1, dH = 0.01, hA = 50,
  # Parasite pars
  f_p =1, hH = 10, QPx = 0.25, QPn = 0.20, dP = 0.01, dE = 0.1, sigma=0.9, 
  alpha = 0.001, k = 10)

inits = c(N = 10, A = 1, QA = 0.02, H = 0.1, QH = 0.15, P = 0.1, QP = 0.22)
run = lsoda(y = inits, times=0:2000, parms = params, func=Worm_stoich, atol=1e-9, rtol=1e-9)
plot(run)

Total_N = run[,"N"] + run[,"A"]*run[,"QA"] + run[,"H"]*run[,"QH"] + run[,"P"]*run[,"QP"]
var(Total_N)
plot(run[,"time"], Total_N, typ="l")

# 
# N.out = numeric()
# A.out = numeric()
# QA.out = numeric()
# H.out = numeric()
# QH.out = numeric()
# P.out = numeric()
# 
# ns = seq(from=0.0001, to=0.2, length.out = 100)
# for (n in 0:100){
#   inits["N"] = n
#   output = lsoda(y = inits, times=0:100000, parms = params, func=Worm_stoich)
#   N.out[n] = output[100001, 2]
#   A.out[n] = output[100001, 3]
#   QA.out[n] = output[100001, 4]
#   H.out[n] = output[100001, 5]
#   QH.out[n] = output[100001, 6]
#   P.out[n] = output[100001, 7]
# }
# par(mfrow=c(2, 3))
# plot(1:100, N.out, typ="l")
# plot(1:100, A.out, typ="l")
# plot(1:100, QA.out, typ="l")
# plot(1:100, H.out, typ="l")
# plot(1:100, QH.out, typ="l")
# plot(1:100, P.out, typ="l")