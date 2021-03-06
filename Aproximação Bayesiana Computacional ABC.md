<p>
Estimador de máxima verosimilhança aproximado por Aproximação Bayesiana
Computacional ABC
</p>

Monte Carlo MC
--------------

    erro <- 0.01
    amostra <- rexp(100, 1)#Gerando amostra
    priori <- runif(1000,min=0,max=2)#Priori para theta
    resultados <- c()#Resultados de cada iteração

    for(i in 1:1000){
      amostra_sim <- rexp(100,rate = priori[i])
      if(abs(1/mean(amostra) - 1/mean(amostra_sim)) < erro){
        resultados[i] <- "Valor(es) aceito(s)"
      }else{
        resultados[i] <- "Valor(es) não aceito(s)"
      }
    }
    table(resultados)

    ## resultados
    ##     Valor(es) aceito(s) Valor(es) não aceito(s) 
    ##                       6                     994

Movimento Browniano
-------------------

    Mi = -1 # Valor para Mi
    X = rep(0,1000) # Armazena a cadeia
    X[1] = 0 # X inicial
    for(i in 2:1000){ # t = 1000
      X[i] = abs(X[i-1]+Mi+rnorm(1,0,1)) # Gerando a Cadeia de Markov
    }

    amostra <- X[501:1000]
    erro <- 0.01
    priori <- runif(1000,min=0,max=2)#Priori para theta
    resultados <- c()#Resultados de cada iteração

    for(i in 1:1000){
      amostra_sim <- rexp(100,rate = priori[i])
      if(abs(1/mean(amostra) - 1/mean(amostra_sim)) < erro){
        resultados[i] <- "Valor(es) aceito(s)"
      }else{
        resultados[i] <- "Valor(es) não aceito(s)"
      }
    }
    table(resultados)

    ## resultados
    ##     Valor(es) aceito(s) Valor(es) não aceito(s) 
    ##                       7                     993

Metropolis-Hastings
-------------------

    target = function(x){
      if(x<0){
        return(0)}
      else {
        return( exp(-x))
      }
    }
    X = rep(0,1000)
    X[1] = 3     #Valor inicial
    for(i in 2:1000){
      currentx = X[i-1]
      proposedx = currentx + rnorm(1,mean=0,sd=1)
      A = target(proposedx)/target(currentx) 
      if(runif(1)<A){
        X[i] = proposedx       # Probabilidade de ser aceito min(1,A)
      } else {
        X[i] = currentx        # Rejeitar
      }
    }
    amostra <- X
    erro <- 0.01
    priori <- runif(1000,min=0,max=2)#Priori para theta
    resultados <- c()#Resultados de cada iteração

    for(i in 1:1000){
      amostra_sim <- rexp(100,rate = priori[i])
      if(abs(1/mean(amostra) - 1/mean(amostra_sim)) < erro){
        resultados[i] <- "Valor(es) aceito(s)"
      }else{
        resultados[i] <- "Valor(es) não aceito(s)"
      }
    }
    table(resultados)

    ## resultados
    ##     Valor(es) aceito(s) Valor(es) não aceito(s) 
    ##                       7                     993

Gibbs Sampler
-------------

    X<-rep(0,501) # auxiliar
    X[1]<-3 # Valor inicial
    for(i in 2:501){
      y<-runif(1,0,exp(-X[i-1]))
      X[i]<-runif(1,0,-log(y))
    }

    amostra <- X
    erro <- 0.01
    priori <- runif(1000,min=0,max=2)#Priori para theta
    resultados <- c()#Resultados de cada iteração

    for(i in 1:1000){
      amostra_sim <- rexp(100,rate = priori[i])
      if(abs(1/mean(amostra) - 1/mean(amostra_sim)) < erro){
        resultados[i] <- "Valor(es) aceito(s)"
      }else{
        resultados[i] <- "Valor(es) não aceito(s)"
      }
    }
    table(resultados)

    ## resultados
    ##     Valor(es) aceito(s) Valor(es) não aceito(s) 
    ##                       8                     992
