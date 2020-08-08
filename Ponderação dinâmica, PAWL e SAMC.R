#Ponderação dinâmica, PAWL e SAMC

## f(x) = (1/7,2/7,4/7)
```
if(!require(mvtnorm)){ install.packages("mvtnorm"); require(mvtnorm) }
if(!require(PAWL)){ install.packages("PAWL"); require(PAWL) }
if(!require(SAMCpack)){ install.packages("SAMCpack"); require(SAMCpack) }
if(!require(ggplot2)){ install.packages("ggplot2"); require(ggplot2) }

DistribuicaoX <- function(x){ #Função para a distribuição dos estados de X
  ifelse(x == 1, return(1/7), ifelse(x == 2, return(2/7),return(4/7)))
}

MatrizTransicaoA <- matrix(c(0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow = 3, ncol = 3)
MatrizTransicaoB <- matrix(c(0,0.5,4/7,1/3,0,3/7,2/3,0.5,0), nrow = 3, ncol = 3)
```
## Ponderação Dinâmica

#### Matriz de transição A
```
X_a <- c(1, rep(NA, 2999))
PesoW_a <- c(1, rep(NA, 2999))
for(i in 2:3000){
  x_prop <- sample(c(1,2,3), size = 1, prob = MatrizTransicaoA[X_a[i-1],])
  TaxaDinamicaA <- PesoW_a[i-1]*(DistribuicaoX(x_prop)/DistribuicaoX(X_a[i-1]))
  a <- TaxaDinâmicaA/(1+TaxaDinamicaA)
  u <- runif(n = 1, min = 0, max = 1)
  if(u <= a){
    X_a[i] <- x_prop
    PesoW_a[i] <- TaxaDinamicaA/a
  }else{
    X_a[i] <- X_a[i-1]
    PesoW_a[i] <- PesoW_a[i-1]/(1-a)
  }
}
```
#### Matriz de transição B
```
X_b <- c(1, rep(NA, 2999))
PesoW_b <- c(1, rep(NA, 2999))
for(i in 2:3000){
  x_prop <- sample(c(1,2,3), size = 1, prob = MatrizTransicaoB[X_b[i-1],])
  TaxaDinamicaB <- PesoW_b[i-1]*(DistribuicaoX(x_prop)/DistribuicaoX(X_b[i-1]))
  b <- TaxaDinamicaB/(1+TaxaDinamicaB)
  u <- runif(n = 1, min = 0, max = 1)
  if(u <= a){
    X_b[i] <- x_prop
    PesoW_b[i] <- TaxaDinamicaB/b
  }else{
    X_b[i] <- X_b[i-1]
    PesoW_b[i] <- PesoW_b[i-1]/(1-b)
  }
}
```
# Comparando os erros padrões
```
cat("Matriz de transição A: E(X) = ",mean(X_a), "Erro padrão da estimativa = ",var(X_a),"\n")
cat("Matriz de transição B: E(X) = ",mean(X_b), "Erro padrão da estimativa = ",var(X_b),"\n")
par(mfrow=c(1,2))
hist(X_a, main="Matriz de transição A", xlab = "X", ylab = "Frequência")
hist(X_b, main="Matriz de transição B", xlab = "X", ylab = "Frequência")
```
# Amostra da variável X
```
DistribuicaoAlvo <- c(1/7,2/7,4/7)
```
#Log da probabilidade alvo
```
parametros <- list(logpi = log(DistribuicaoAlvo))
log_densidade <- function(x, parametros){
  return(parametros$logpi[x])
}

proposal_param_A <- list(matriz_transicao = MatrizTransicaoA, card = 3)
proposal_param_B <- list(matriz_transicao = MatrizTransicaoB, card = 3)
```
#Função para gerar as propostas
```
rproposal <- function(states, parametro_proposto){
  for(index in 1:length(states)){
    states[index] <- sample(x = 1:parametro_proposto$card, size = 1,
                            prob = parametro_proposto$matriz_transicao[states[index], ])
  }
  return(list(states = states))
}
```
# Algoritmo - PAWL
```
dproposal <- function(states, ys, parametro_proposto){#Função para calcular a densidade do núcleo proposto
  for (index in 1:(length(states))){
    states[index] <- log(parametro_proposto$matriz_transicao[states[index], ys[index]])
  }
  return(states)
}

proposal_instance_A <- proposal(rproposal = rproposal,
                                dproposal = dproposal,
                                proposalparam = proposal_param_A)

proposal_instance_B <- proposal(rproposal = rproposal,
                                dproposal = dproposal,
                                proposalparam = proposal_param_B)
```
#Função para criar os pontos de partidos do algoritmo MCMC
```
rinit <- function(size) return(rep(1, size))
```
#Definindo o alvo
```
discretetarget <- target(name = "discrete toy example", dimension = 1, type = "discrete",
                         rinit = rinit, logdensity = log_densidade, parameters = parametros)
```
#Especificando os tuning paramethers do M-H
```
mhparameters <- tuningparameters(nchains = 1, niterations = 10000, storeall = TRUE)

getPos <- function(points, logdensity) 2 - (points <= 2)
positionbinning <- binning(position = getPos,
                           name = "position",
                           bins = c(1, 2),
                           desiredfreq = c(0.8, 0.2),
                           useLearningRate = FALSE)

pawlresults_A <- pawl(discretetarget, binning = positionbinning,
                      AP = mhparameters, proposal_instance_A)
pawlresults_B <- pawl(discretetarget, binning = positionbinning,
                      AP = mhparameters, proposal_instance_B)

pawlchains_A <- ConvertResults(pawlresults_A)
pawlchains_B <- ConvertResults(pawlresults_B)

plot(pawlchains_A, main="Matriz de transição A")
plot(pawlchains_B, main="Matriz de transição B")
```
# Exemplo utilizando o pacote SAMC
```
# Step 1 : Define negative log-density function as an R function
func_r = function(x){
  x1 = x[1]; x2 = x[2];
  val1 = (-(x1*sin(20*x2)+x2*sin(20*x1))^2)*cosh(sin(10*x1)*x1);
  val2 = (-(x1*cos(10*x2)-x2*sin(10*x1))^2)*cosh(cos(20*x2)*x2);
  return(val1+val2);
}

## Step 2 : Prepare a setting option
myoption = list()
myoption$partition = c(-Inf,seq(from=-8,to=0,length.out=41))
myoption$tau = 1.0
myoption$domain = c(-1.1,1.1)
myoption$vecpi = as.vector(rep(1/41,41))
myoption$niter = 20000
myoption$stepsize = c(0.25, 10)
## Step 3 : Run The Code
res = SAMC(2,func_r,options=myoption)
## Step 4 : Visualize
select = seq(from=101,to=myoption$niter,by=100) # 100 burn-in, 1/100 thinning
par(mfrow=c(1,2))
plot(res$samples[select,1], res$samples[select,2],xlab='x',ylab='y',main='samples')
barplot(as.vector(res$frequency/sum(res$frequency)),
        main="visiting frequency by energy partition",
        names.arg=myoption$partition[-1], xlab="energy")
```
