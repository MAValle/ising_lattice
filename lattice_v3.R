# HACEMOS Lo mismo que en lattice_v2.R, pero esta vez, hacemos
# metropolis sweep muchas veces a una misma temperatura antes de 
# pasar a la siguiente iteracion. Ademas agregamos muchas iteraciones
# de metropolis sweep de precalentamiento.
# Se hace de manera similar a lo indicado en # https://rajeshrinet.github.io/blog/2014/ising-model/

# Notar que ahora, no se selecciona un nodo al azar, si no que
# se van recorriendo cada uno de los nodos uno por uno, y a cada
# uno se le hace un numero de MH sweeps a una temperatura determinada.

#File name: lattice_v3.R
#Creation date: 09-abr-18


# Notas:
# 09-abr-20: creation

# resources:
# https://repl.it/languages/python3 # emulador de python online
# https://rajeshrinet.github.io/blog/2014/ising-model/


# paremeters of the system
m = 10
J = 1 # acoples 
Temp_start = 3.28
Temp_end = 1.53
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
N <- length(V(g)) # Numeber of nodes



# parameters of the simulation
sims <- 88 # number of main iterations
eqSims <- 1020 # 1020 number of MH sweeps for equilibrium
mcSims <- 1020 # 1020 number of MH sweeps for calculation

system_energy <- numeric(length=sims)
magn <- numeric(length = sims)
#probs_transition <- numeric(length = sims)
#system_states <- matrix(NA, ncol=length(V(g)), nrow=sims)

temps <- sort(seq(from=Temp_start, to=Temp_end, length.out = sims), decreasing = T)
beta = temps*kb

# # # START SIMULATION---------------------------------------
require(svMisc)
for (i in 1:sims) {
  progress(i)
  # create randon spin
  s <- sample(c(-1,1), length(V(g)), replace=T, prob=c(0.5,0.5))
  E1 <- M1 <- 0 #almacenamiento de energia y magnetizacion
  
  # MH sweeps for equilibrium
  for (h in 1:eqSims) {
    s <- mh_sweeps(g=g, N=N, beta = beta[i], s=s)
  }
  
  # MH sweeps for computation
  for (k in 1:mcSims) {
    s <- mh_sweeps(g=g, N=N, beta = beta[i], s=s)
    Ene <- get_energySys(s=s, J = J, N = length(V(g)), g = g )
    Mag <- get_magnetization(s)
    
    E1 <- E1 + Ene
    M1 <- M1 + Mag
  }
  
  # calculo del promedio de energia
  system_energy[i] <- E1/(mcSims*N*N)
  # calculo del promedio de magnetizacion
  magn[i] <- M1/(mcSims*N*N)
}
plot(1:sims, system_energy, type = "l", col = "red", xlab = "epoch", ylab = "U")
plot(1:sims, magn, type = "l", pch = 19, col = "red", xlab = "epoch", ylab = "M")
plot(1:sims, probs_transition, type = "l", pch = 19, col = "black", xlab = "epoch", ylab = "P->")
plot(probs_transition, magn, type = "l", pch = 19, col = "black", xlab = "epoch", ylab = "P->")
# # # END SIMULATION------------------------------------------


# # # # # METROPOLIS HASTING SWEEPS
# FOR ALL NODES ON A NETWORK
# inputs: 
# * g = objeto igraph
# * N = V(g)  todos los nodos de g.
# * s = estado del sistema de cada nodo 
# * beta = valores de beta = temperatura en el momento t * kb
# salida:
# s: estado del sistema
mh_sweeps <- function(g, N, beta, s) {
  for (nn in 1:N) {
    
    # choose node at random 
    node <- as.numeric(sample(V(g), 1))
    node_state_original <- s[node]
    node_state <- s[node]
    
    # buscamos los nodos adyacentes de node
    ady <- unlist(adjacent_vertices(g, node ))
    
    # calculo de diferencia de energia
    dE <- get_energyDif(s=s, id = ady, node = node, J = J )
    
    # calculo del nuevo estado del sistema
    #r <- transition_metropolis(s=s, node = node, delta_U = dE , beta = beta[i])
    #new_s <- r$state
    # Esto es nueva mas simple, y reemplaza la funcion transition_metropolis
    if (dE < 0) {
      node_state <- -1*node_state
    } else {
      if (runif(1) < exp(-dE/beta) ) {
        node_state <- -1*node_state
      }
    }
    s[node] <- node_state
    
    
    # cat("Node:",node," dE:", dE,  "\r")
    # cat("Node:",node," dE:", dE,  "\b")
    #cat("Node: ",node, " state or: ", node_state_original, " state end: ", node_state,   " dE: ", dE,  "\n")
  }
  return(s)
}
# ejemplo
# mh_sweeps(g=g, N=N, beta = beta[i], s=s)
# # # # # METROPOLIS HASTING SWEEPS






# # # # # STATE CHANGE SIMULATION
# ONLY FOR ONE NODE!!!
# entrada :
# * s = estado del sistema de cada nodo 
# node = nodo seleccionado
# delta_U: energy difference computed from get_energyDif
# * beta = valores de beta = temperatura en el momento t * kb
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
  
  #cat("Node: ",node," dE: ", dE,  "\n")
  cat("dU: ",delta_U, " Node: ", node, " beta: ", beta, "\n")
  
  # metropolis
  if (delta_U < 0) { 
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
# r <- transition_metropolis(s=s, node = node, delta_U = dE , beta = beta[i])
# r$state
# r$transition_prob
# # # # # ENERGY DIFFERENCE COMPUTATION FUNCTION























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
  deltaE <-  sum(2*s_node*s_neig )
  
  return(deltaE)
}
# Example
# get_energyDif(s=s0, id = ady, node = node, J = J )
# # # # # ENERGY DIFFERENCE COMPUTATION FUNCTION







# # # # # MAGNETIZATION COMPUTATION FUNCTION
# input: s = state vector of the system
get_magnetization <- function(s) {
  mag <- sum(s)
  return(mag)
}
# # # # # MAGNETIZATION COMPUTATION FUNCTION