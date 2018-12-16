############################################################
#                                                          #
#  Universidade Federal do Pampa                           #
#                                                          #
#  Ajustamento de Observações Goedésicas                   #
#                                                          #
#  Prova Prática  Recuperatória II - II sem/2018           #
#                                                          #
############################################################

## Organização dos dados 

# cota M 

M <- 348.567

# Observações

n <- 6

# incógnitas 

u <- 3

# distâncias em km; 

dn1 <- 450.31/1000
dn2 <- 422.12/1000
dn3 <- 413.14/1000
dn4 <- 256.22/1000
dn5 <- 261.48/1000
dn6 <- 222.98/1000

#########################################################
## Vetor das distâncias (transformadas em quilomêtros)

d <- c(dn1, dn2, dn3, dn4, dn5, dn6) 
d
#########################################################
#########################################################

########### Vetor das observações (desníveis);
##############################################
### Vetor Com os valores Observados (Não Ajustados - medidos)

l <- c(1.111,1.111,-2.222, 3.333,  2.222, 1.111)
l
#########################################################
#########################################################
# Matriz dos coeficientes das incógnitas 

# Equações

# l1a = B - A
# l2a = C - B
# l3a = A - C
# l4a = M - A 
# l5a = M - B 
# l6a = M - C  

#########################################
#########################################
# Ordem - BA CB AC MA MB MC 

BA <- c(-1, 1, 0)
CB <- c(0, -1, 1)   
AC <- c(1, 0, -1)    
MA <- c(-1, 0, 0)    
MB <- c(0, -1, 0)   
MC <- c(0, 0, -1)  

A <- rbind(BA, CB, AC, MA, MB, MC)
colnames(A) <- c( "A", "B", "C")
A
## Vetor dos parametros iniciais ("chutados")......

H0.A <- 12
H0.B <- 4
H0.C <- 1

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

#Vetor L dos valores observados (l = l0 - lb)
l <- l0-lb
l
#Matriz dos pesos P (nxn) = var.pri * solve(pp1)

#Variancia de unidade de peso a priori

var.pri <- 1

#######################################
#######################################

#### Matriz de pesos com valores iguais

pesos_iguais <- diag((vector(mode = "numeric", length = n))+1)
pesos_iguais
#### Matriz com valores diferentes

pesos1 <- diag((vector(mode = "numeric", length = n))+1/sqrt(4*100*d))
pesos1

pesos_iguais==pesos1

Pesos=var.pri * solve(pesos_iguais)
Pesos

Pesos_1 = var.pri * solve(pesos1)
Pesos_1

Pesos==Pesos_1

#Vetor das correcoes (x):

N <- t(A)%*%Pesos%*%A
N

U <- t(A)%*%Pesos%*%l
U

x <- -(solve(N))%*%U
x 

#### Vetor dos parametros observados ajustados
xa <- x0 + x
xa

#### Vetor dos residuos das observacoes (v)
v <- (A%*%x)+l
v

#### Vetor dos valores observados ajustados (la)
la <- lb + v
la

######################################
######################################

#### MEDIDAS DE QUALIDADE DO AJUSTAMENTO
#1?) Variancia da unidade de peso a posteriori (V^T*P*v)/GL
GL <- n-u
GL

var.pos <- (t(v)%*%Pesos%*%v)/GL
var.pos

# 2) Matriz Variancia-Covariancia dos parametros observados ajustados (MVCxa)

MVCX_ajustado <- var.pri * solve(N)
MVCX_ajustado

sd(MVCX_ajustado)

# 3) Matriz Variancia-Covariancia dos valores observados ajustados (MVCla)

MVC_l_ajustado <- var.pri * (A%*%solve(N)%*%t(A))
MVC_l_ajustado
sd(MVC_l_ajustado)

## Não houve covariância 
cov(MVC_l_ajustado) 

# 4) Matriz Variancia-Covariancia dos residuos (MVCv)

MVC_var_cov <- var.pri * (solve(Pesos))- MVC_l_ajustado
MVC_var_cov
cov(MVC_var_cov)

# 5) Comparacao da Variancia a priori com posteriori
qui.calculado <- (var.pos/var.pri)*(GL)
alfa.1 <- 0.025
tabelado.a.1 <- 0.484
alfa.2 <- 0.975
tabelado.a.2 <- 14.860

tabelado.a.1 < qui.calculado 
tabelado.a.2 > qui.calculado 

#Como tabelado.a.1 deu false, aceita-se a hipotese H1.


#LOCALIZAÇÃO DOS ERROS NAS OBSERVACOES PELO TESTE DATA SNOOPING DE BAARDA:

#NUMEROS DE REDUNDÂNCIA:
R <- (1/var.pri)*MVC_l_ajustado%*%Pesos
R

#Traco da Matriz R:
trR <-  (R[1,1] + R[2,2] + R[3,3] + R[4,4] + R[5,5] + 
           R[6,6]) 
trR

### Tr deu igual Gl, portanto os calculos estao corretos. Entretanto, o ajustamento nao foi satisfatorio.


#Estatistica do teste:
#Criando a matriz onde serao armazenados os resultados

W <- matrix(data = 0, nrow = n, ncol = 3)
W[,2] <- 1.96
colnames(W) <- c("Estatistica", "K(95%)","Resultado")
for(i in 1:n){W[i,1] <- abs (v[i] / (sqrt(MVC_l_ajustado[i,i])*sqrt(R[i,i])))
if (W[i,1] > W[i,2]) W[i,3] ="H1"   else W[i,3] ="H0"}
W
