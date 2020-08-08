O estimador de máxima verosimilhança aproximado por Monte Carlo MCEMV
<p>
Gerando dados a partir de:
</p>

Monte Carlo MC
--------------

<p>
Utilizando a função vista em sala de aula para otimizar a
Log-verossimilhança
</p>

    xamostra <- rexp(100, 1)  # amostra inicial de tamanho 100
    thetaestimado <- seq(0, 2, length.out = 100)  # possíveis valores para theta 
    verossimilhanca <- c()# vetor para armazenar a verossimilhança de cada theta estimado
    for (i in 1:length(thetaestimado)) {
      x_eta <- rexp(500, thetaestimado[i])# Gerando 500 valores seguindo distribuição exponencial(theta estimado)
      teta_eta <- 1/mean(x_eta) # solução analítica para achar theta
      verossimilhanca[i] <- log(sum(exp(-xamostra))) -100*log((1/500)*sum((teta_eta*exp(-teta_eta*x_eta))/exp(-x_eta)))
    }
    thetaestimado[which.max(verossimilhanca)] #valor máximo da verossimilhança

    ## [1] 0.9494949

<p>
Utilizando a função optim para otimizar a Log-verossimilhança
</p>

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-3-1.png)

Movimento Browniano
-------------------

    Mi = -1 # Valor para Mi
    X = rep(0,1000) # Armazena a cadeia
    X[1] = 0 # X inicial
    for(i in 2:1000){ # t = 1000
      X[i] = abs(X[i-1]+Mi+rnorm(1,0,1)) # Gerando a Cadeia de Markov
    }
    hist(X, main = "Histograma de Xt",ylab = "Frequência",ylim = c(0,400)) # Histograma da Cadeia
    abline(v=1,col="red")

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-4-1.png)

    X = data.frame(X = X[501:1000]) #Escolhendo as 500 observações
    ggplot(data = X, mapping = aes(sample = X)) +
      stat_qq_band(distribution = "exp", dparams = list(rate = 1)) +
      stat_qq_line(distribution = "exp", dparams = list(rate = 1)) +
      stat_qq_point(distribution = "exp", dparams = list(rate = 1)) +
      labs(x = "Quantis definidos", y = "Quantis amostrados")

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-4-2.png)

<p>
Utilizando a função vista em sala de aula para otimizar a
Log-verossimilhança
</p>

    xamostra <- X 
    thetaestimado <- seq(0, 2, length.out = 100)  # possíveis valores para theta 
    verossimilhanca <- c()# vetor para armazenar a verossimilhança de cada theta estimado
    for (i in 1:length(thetaestimado)) {
      x_eta <- rexp(500, thetaestimado[i])# Gerando 500 valores seguindo distribuição exponencial(theta estimado)
      teta_eta <- 1/mean(x_eta) # solução analítica para achar theta
      verossimilhanca[i] <- log(sum(exp(-xamostra))) -100*log((1/500)*sum((teta_eta*exp(-teta_eta*x_eta))/exp(-x_eta)))
    }
    thetaestimado[which.max(verossimilhanca)] #valor máximo da verossimilhança

    ## [1] 1.070707

<p>
Utilizando a função optim para otimizar a Log-verossimilhança
</p>

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-6-1.png)

Metropolis-Hastings (Seguindo referência)
-----------------------------------------

    target = function(x){
      if(x<0){
        return(0)}
      else {
        return( exp(-x))
      }
    }
    X = rep(0,1000)
    X[1] = 3     #this is just a starting value, which I've set arbitrarily to 3
    for(i in 2:1000){
      currentx = X[i-1]
      proposedx = currentx + rnorm(1,mean=0,sd=1)
      A = target(proposedx)/target(currentx) 
      if(runif(1)<A){
        X[i] = proposedx       # accept move with probabily min(1,A)
      } else {
        X[i] = currentx        # otherwise "reject" move, and stay where we are
      }
    }
    X<-data.frame(X=X)
    ggplot(data = X, mapping = aes(sample = X)) +
      stat_qq_band(distribution = "exp", dparams = list(rate = 1)) +
      stat_qq_line(distribution = "exp", dparams = list(rate = 1)) +
      stat_qq_point(distribution = "exp", dparams = list(rate = 1)) +
      labs(x = "Quantis definidos", y = "Quantis amostrados")

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-7-1.png)

<p>
Utilizando a função vista em sala de aula para otimizar a
Log-verossimilhança
</p>

    xamostra <- X 
    thetaestimado <- seq(0, 2, length.out = 100)  # possíveis valores para theta 
    verossimilhanca <- c()# vetor para armazenar a verossimilhança de cada theta estimado
    for (i in 1:length(thetaestimado)) {
      x_eta <- rexp(500, thetaestimado[i])# Gerando 500 valores seguindo distribuição exponencial(theta estimado)
      teta_eta <- 1/mean(x_eta) # solução analítica para achar theta
      verossimilhanca[i] <- log(sum(exp(-xamostra))) -100*log((1/500)*sum((teta_eta*exp(-teta_eta*x_eta))/exp(-x_eta)))
    }
    thetaestimado[which.max(verossimilhanca)] #valor máximo da verossimilhança

    ## [1] 1.070707

<p>
Utilizando a função optim para otimizar a Log-verossimilhança
</p>

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-9-1.png)

Gibbs Sampler
-------------

    X<-rep(0,501) # auxiliar
    X[1]<-3 # Valor inicial
    for(i in 2:501){
      y<-runif(1,0,exp(-X[i-1]))
      X[i]<-runif(1,0,-log(y))
    }
    hist(X, main = "Histograma de Xt",ylab = "Frequência",ylim = c(0,400)) # Histograma da Cadeia
    abline(v=1,col="red")

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-10-1.png)

    X = data.frame(X=X)
    ggplot(data = X, mapping = aes(sample = X)) +
      stat_qq_band(distribution = "exp", dparams = list(rate = 1)) +
      stat_qq_line(distribution = "exp", dparams = list(rate = 1)) +
      stat_qq_point(distribution = "exp", dparams = list(rate = 1)) +
      labs(x = "Quantis definidos", y = "Quantis amostrados")

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-10-2.png)
<p>
Utilizando a função vista em sala de aula para otimizar a
Log-verossimilhança
</p>

    xamostra <- X 
    thetaestimado <- seq(0, 2, length.out = 100)  # possíveis valores para theta 
    verossimilhanca <- c()# vetor para armazenar a verossimilhança de cada theta estimado
    for (i in 1:length(thetaestimado)) {
      x_eta <- rexp(500, thetaestimado[i])# Gerando 500 valores seguindo distribuição exponencial(theta estimado)
      teta_eta <- 1/mean(x_eta) # solução analítica para achar theta
      verossimilhanca[i] <- log(sum(exp(-xamostra))) -100*log((1/500)*sum((teta_eta*exp(-teta_eta*x_eta))/exp(-x_eta)))
    }
    thetaestimado[which.max(verossimilhanca)] #valor máximo da verossimilhança

    ## [1] 0.8888889

<p>
Utilizando a função optim para otimizar a Log-verossimilhança
</p>

![](EMVmarkdown_files/figure-markdown_strict/unnamed-chunk-12-1.png)
