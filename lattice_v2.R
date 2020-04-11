# HACEMOS Lo mismo que en lattice_v1.R, pero esta vez, hacemos
# metropolis sweep. Un metropolis sweep es pasar por cada nodo 
# adyacente del nodo elegido (incluyendo est ultimo), para determinar
# el estado final del sistema, antes de elegir aleatoriamente otro nodo.


#File name: lattice_v2.R
#Creation date: 07-abr-18


# Notas:
# 07-abr-20: creation

# resources:
# https://repl.it/languages/python3
# https://rajeshrinet.github.io/blog/2014/ising-model/


# paremeters of the system
m = 10
J = 1 # acoples 
Temp = 100  # initial temperature: must finish in T=1
kb = 1 # Boltzmann constant
#



# Generacion de un lattice de m X m
library(igraph)
#g <- make_lattice(length = 4, dim = 2, circular = TRUE)
g <- make_lattice(length = m, dim = 2)
plot(g)

# assign V(g) -1,1 for each node at random
#V(g)$spin <- sample(c(-1,1), length(V(g)), replace=T, prob=c(0.5,0.5))
# choose random state
V(g)$spin <- rep(1, length(V(g)))
s <- V(g)$spin


# parameters of the simulation
sims <- 1000 #number of metropolis iterations
system_energy <- numeric(length=sims)
magn <- numeric(length = sims)
probs_transition <- numeric(length = sims)
system_states <- matrix(NA, ncol=length(V(g)), nrow=sims)

temps <- sort(seq(from=1, to=Temp, length.out = sims), decreasing = T)
beta = temps*kb


# # # START SIMULATION---------------------------------------
require(svMisc)
for (i in 1:sims) {
  progress(i)
  
  # para hacer mas rapido el proceso, puedo evitar almacenar los estados del sistema.
  #system_states[i, ] <- s
  
  # calculo de la energia del sistema
  U <- get_energySys(s=s, J = J, N = length(V(g)), g = g )
  system_energy[i] <- U
  
  # calculo de la magnetizacion del sistema
  #magn[i] <- mean(s)
  magn[i] <- get_magnetization(s)
  
  # choose node at random 
  node <- as.numeric(sample(V(g), 1))
  
  # buscamos los nodos adyacentes de node
  ady <- unlist(adjacent_vertices(g, node ))

  
  # METROPOLIS SWEEP over the adjacent nodes included selected node
  id_involved <- c(node, ady)
  n_inv <- length(id_involved)
  for (j in 1:n_inv) {
    
    # calculo de diferencia de energia
    dE <- get_energyDif(s=s, id = id_involved[-j], node = id_involved[j], J = J )
    
    # calculo del nuevo estado del sistema
    r <- transition_metropolis(s=s, node = id_involved[j], delta_U = dE , beta = beta[i])
    s <- r$state
  }
  
  # no guardamos la probabilidad de transisicion
  probs_transition[i] <- r$transition_prob
  
}
plot(1:sims, system_energy, type = "l", col = "red", xlab = "epoch", ylab = "U")
plot(1:sims, magn, type = "l", pch = 19, col = "red", xlab = "epoch", ylab = "M")
plot(1:sims, probs_transition, type = "l", pch = 19, col = "black", xlab = "epoch", ylab = "P->")

plot(probs_transition, magn, type = "l", pch = 19, col = "black", xlab = "epoch", ylab = "P->")

#https://rajeshrinet.github.io/blog/2014/ising-model/











# FUNCTIONS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # ENERGY CMPUTATION OF THE SYSTEM
# entrada :
# s = estado del sistema de cada nodo 
# J = valor de los acoples
# N = numero de nodos por ejemplo: length(V(g))
# g = objeto igraph red
# salida
# E energia del sistema
get_energySys <- function(s, J, N, g) { 
  sumenergy = 0
  for (i in 1:N ) {
    # buscamos los nodos adyacentes de node i
    ady <- unlist(adjacent_vertices(g, i ))
    energy <- -1*sum(s[ady]*s[i])/2 # se divide por dos para considerar el otro nodo que tiene el mismo edge
    sumenergy <- sumenergy + energy
  }
  return(sumenergy)
}
# Example
# get_energySys(s=s0, J = J, N = length(V(g)), g= g)
# # # # # ENERGY CMPUTATION OF THE SYSTEM




# # # # # ENERGY DIFFERENCE COMPUTATION FUNCTION
# # nota:
# El calculo de la energia se calcula solo sobre los vecinos del nodo seleccionado k.
# La diferencia de energia E(S)-E(S') siendo S=estado antes de los vecinos del nodo
# k, y S' estado despues de los vecinos de k, incluyendo nodo k. 
# Esta diferencia DE se calcula solo tomando encuenta las conecciones de nodo k con
# sus vecino, porque las relaciones entre nodos vecinales se cancelan.
# O sea: E = -J * sum_{sobre los vecinos de k} s_k*s_i  donde i son todos los vecinos de k.

# entrada :
# s = estado del sistema de cada nodo 
# id = id de los nodos que son vecinos del nodo k seleccionado
# node = nodo seleccionado
# J = valor de los acoples
# salida
# deltaE diferencia de energia
get_energyDif <- function(s, id, node, J) { 
  # num_sums <- length(ady)
  # get the states of the nodes
  s_node <- s[node]
  # s_node_flipped <- s_node*(-1) 
  s_neig <- s[c(id)]
  # energia inicial
  #E0 <- -J*sum(s_node*s_neig)
  #E1 <- -J*sum(s_node_flipped*s_neig)
  #deltaE <- (E1 - E0) 
  # las tres linean anteriores es lo mismo que:
  deltaE <-  2*s_node*s_neig 
  
  return(deltaE)
}
# Example
# get_energyDif(s=s0, id = ady, node = node, J = J )
# # # # # ENERGY DIFFERENCE COMPUTATION FUNCTION




# # # # # STATE CHANGE SIMULATION
# entrada :
# s = estado del sistema de cada nodo 
# node = nodo seleccionado
# delta_U: energy difference computed from get_energyDif
# beta = valores de beta = temperatura en el momento t * kb
# salida
# tptob = transition probability
# sout = nuevo estado del sistema
transition_metropolis <- function(s, node, delta_U, beta) { 
  # num_sums <- length(ady)
  # get the states of the node
  s_node <- s[node]
  s_node_flipped <- s_node*(-1)
  # system state in case of aprove:
  s_out_changed <- s
  s_out_changed[node] <- s_node_flipped 
  
  # metropolis
  if (delta_U <= 0) { 
    tprob <- 1
    sout <- s_out_changed 
  } else {
    #tprob <- exp(-delta_U/beta)
    tprob <- min(1, exp(-delta_U/beta) )
    u <- runif(1)
    if (u < tprob) {
      #sout <- s
      sout <- s_out_changed
    } else {
      #sout <- s_out_changed
      sout <- s
    }
  }
  
  return(list(state = sout, transition_prob = tprob))
}
# Example
# r <- transition_metropolis(s=s0, node = node, delta_U = dE , beta = beta[2])
# r$state
# r$transition_prob
# # # # # ENERGY DIFFERENCE COMPUTATION FUNCTION


# # # # # MAGNETIZATION COMPUTATION FUNCTION
# input: s = state vector of the system
get_magnetization <- function(s) {
  mag <- sum(s)
  return(mag)
}
# # # # # MAGNETIZATION COMPUTATION FUNCTION