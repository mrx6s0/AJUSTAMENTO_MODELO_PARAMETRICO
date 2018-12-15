############################################################
#                                                          #
#  Universidade Federal do Pampa                           #
#                                                          #
#  Ajustamento de Observações Goedésicas                   #
#                                                          #
#  Prova Prática              
#                                                          #
#                                                          #
#                         
############################################################

# distâncias em km; 

dn1 <- 450.31/1000
dn2 <- 422.12/1000
dn3 <- 413.14/1000
dn4 <- 256.22/1000
dn5 <- 261.48/1000
dn6 <- 222.98/1000

# cota M 

M <- 348.567

# Observações

n <- 6

# incógnitas 

u <- 3

#########################################################
## Vetor das distâncias (transformadas em quilomêtros)

D <- c(dn1, dn2, dn3, dn4, dn5, dn6) 
D

# Matriz peso usando identidade ou seja, peso igual
# para todas as observações

P1 <- diag( rep(1,6) )
P1

# Matriz peso das distâncias usando Schmidtz

P2 <-  diag(1/(sqrt(4*100*D)))
P2

# Vetor das observações (desníveis);
##############################################
### Vetor Com os valores Observados (Não Ajustados - medidos)

lb <- c(1.111,1.111,-2.222, 345.234 , 346.345, 347.456)
lb
##############################################
##############################################
# Matriz dos coeficientes das incógnitas 

# Equações

# Ordem em função de lb; 
# 
# l1a = B - A
# l2a = C - B
# l3a = A - C
# l4a = M - A 
# l5a = M - B 
# l6a = M - C  

# Ordem - BA CB AC MA MB MC 

# Comprovando os resultados:

# BA = -(-683.801) + (-674.912)+0 = 8.889
# CB = (
# AC = 8.33 * 1 + 0 + (-27.22) * 1 = -18.89 
BA <- c(-1, 1, 0)
CB <- c(0, -1, 1)   
AC <- c(1, 0, -1)    
MA <- c(-1, 0, 0)    
MB <- c(0, -1, 0)   
MC <- c(0, 0, -1)  

A <- rbind(BA, CB, AC, MA, MB, MC)
A

#########################################################
# TRATA-SE... PARA RESPONDER A PRIMEIRA QUESTÃO TEÓRICA, ONDE PEDE 
# SE A MATRIZ DE PESOS MELHOROU O RESULTADO DO AJUSTAMENTO
# EM COMPARAÇÃO COM A MATRIZ DE PESOS UNITÁRIAS 
# 
##########################################################
# VETOR  DOS VALORES OBERSVADOS PESOS UNITÁRIOS 
X1 <- solve( t(A) %*% P1 %*% A  ) %*% t(A) %*% P1 %*% lb
X1
##########################################################
# VETOR DOS VALORES OBSERVADOS USANDO PESOS USANDO DISTÃNCIAS
X2 <- solve( t(A) %*% P2 %*% A  ) %*% t(A) %*% P2 %*% lb
X2
##########################################################
##############################################
cov(X1,X2) > 0 

## Vetor dos parametros iniciais ("chutados")......

H0.A <- 4
H0.B <- 1
H0.C <- 2

x0 <- c(H0.A, H0.B, H0.C)
x0

##############################################
## Vetor L0 das observações calculadas em função dos valores iniciais ("chutados")

l0.1 <- H0.B - H0.A
l0.2 <- H0.C - H0.B
l0.3 <- H0.A - H0.C
l0.4 <- M - H0.A
l0.5 <- M - H0.B
l0.6 <- M - H0.C

l0 <- c(l0.1, l0.2, l0.3, l0.4, l0.5, l0.6)
l0


##############################################
## Vetor L dos valores observados (l = l0 - lb)

l <- l0-lb
l

##############################################
## Matriz dos pesos P (nxn) = var.pri * solve(pp1)
# Variância de unidade de peso a priori

var.pri <- 1
pp1 <- diag((vector(mode = "numeric", length = n))+1)
P = var.pri * solve(pp1)
pp1
P
A
P
l
##############################################
## Vetor das correções x
N <- t(A) %*% P %*% A
U <- t(A) %*% P %*% l
x <- solve(N) %*% U
N

U
x
##############################################
## Vetor dos parâmetros observados ajustados
xa <- x0 + x
xa

sd(xa)
##############################################
## Vetor dos resíduos das observações (v)
v <- (A%*%x)+l
v
##############################################
## Vetor dos valores observados ajustados (la)
la <- lb + v
lb
la
##############################################
## MEDIDAS DE QUALIDADE DO AJUSTAMENTO GEODÉSICO

## 1º) Variância da unidade de peso a posteriori (V^T*P*v)/GL
GL <- n-u
var.pos <- (t(v)%*%P%*%v)/GL
GL
var.pos
## 2º) Matriz Variância-Covariância dos parâmetros observados ajustados (MVCxa)

MVCxa <- var.pri * solve(N)
MVCxa

## 3º) Matriz Variância-Covariância dos valores observados ajustados (MVCla)
MVCla <- var.pri * (A%*%solve(N)%*%t(A))
MVCla


## 4º) Matriz Variância-Covariância dos resíduos (MVCv)

MVCv <- var.pri * (solve(P))- MVCla
MVCv
cov_between_MVCla_MVCv <- cov(MVCla,MVCv)  
cov_between_MVCla_MVCv
cov_between_MVCla_MVCv >= 0

## 5º) Comparação da Variância a priori com posteriori
qui.calculado <- (var.pos/var.pri)*(GL)
qui.calculado
alfa.1 <- 0.025
tabelado.a.1 <- 0.216
alfa.2 <- 0.995
tabelado.a.2 <- 12.838

tabelado.a.1 < qui.calculado # Se Falso, H1, se verdadeiro testar a condição abaixo
tabelado.a.2 > qui.calculado # Se primeira condição TRUE e Segunda condição TRUE, 
# H0, caso contrário H1
## De uma forma geral, TRUE e TRUE para ser aceita a hipótese H0.

## LOCALIZAÇÃO DOS ERROS NAS OBSERVAÇÕES PELO TESTE DATA SNOOPING DE BAARDA

# NÚMEROS DE REDUNDÂNCIA
R <- (1/var.pri)*MVCv%*%P
R

trR <-  (R[1,1] + R[2,2] + R[3,3] + R[4,4] + R[5,5] + 
           R[6,6] ) 

# Estatística do teste
# Criando a matriz onde serao armazenados os resultados

W <- matrix(data = 0, nrow = n, ncol = 3)
W[,2] <- 1.96
colnames(W) <- c("Estatistica", "K(95%)","Resultado")
colnames(W)

for(i in 1:n){
  
  W[i,1] <- abs (v[i] / (sqrt(MVCla[i,i])*sqrt(R[i,i])))
  if (W[i,1] > W[i,2]) W[i,3] ="H1"   else W[i,3] ="H0" }
W
