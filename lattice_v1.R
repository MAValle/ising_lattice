# Este script es un primer intento para genera runa grilla que representa
# un lattice Ising de N X N .

#File name: lattice_v1.R
#Creation date: 30-dic-18

# Notas:
# 30-dic-18: cration
# 07-apr-20: Actualmente, en cada iteracion el procedimiento consiste
# en seleccionar aleatoriamente un nodo de la red. Luego se identifican
# los nodos adyacentes al nodo seleccionados y se calcula la energia 
# del considerando solo el nodo sleccionado y sus vecinos.
# con los spin actuales. Luego se hace lo mismo paero cambiando la orientacion
# del spin del nodo seleccionado. Luego se calcula la diferencia de energia
# y se determina la proababilidad de transiscion. 
# En cada uno de estos pasos, se rescata Y guarda el estado del sistema entero.
# MI DUDA ES SI DEBO HACER UN METROPOLIS SWEEP A TODOS LOS NODOS ADYACENTES
# DEL NODO SELECCIONADO INCLUYENDO ESTE ULTIMO ANTES DE RESCATAR Y GUARDAR EL ESTADO 
# DEL SISTEMA.



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
V(g)$spin <- sample(c(-1,1), length(V(g)), replace=T, prob=c(0.5,0.5))

# choose random state
# s <- rep(1, length(V(g) ) )
s <- V(g)$spin


# parameters of the simulation
sims <- 3000 #number of metropolis iterations
system_energy <- numeric(length=sims)
magn <- numeric(length = sims)
probs_transition <- numeric(length = sims)
system_states <- matrix(NA, ncol=length(V(g)), nrow=sims)
 
temps <- sort(seq(from=1, to=Temp, length.out = 1000), decreasing = T)
beta = temps*kb


# # # ESTO SE REPITE---------------------------------------
require(svMisc)
for (i in 1:sims) {
  progress(i)
  
  # para hacer mas rapido el proceso, puedo evitar almacenar los estados del sistema.
  #system_states[i, ] <- s
  
  # calculo de la energia del sistema
  U <- get_energySys(s=s, J = J )
  system_energy[i] <- U
  
  # calculo de la magnetizacion del sistema
  #magn[i] <- mean(s)
  magn[i] <- sum(s)
  
  # choose node at random 
  node <- as.numeric(sample(V(g), 1))
  
  # buscamos los nodos adyacentes de node
  ady <- unlist(adjacent_vertices(g, node ))
  
  # calculo de diferencia de energia
  dE <- get_energyDif(s=s, id = ady, node = node, J = J )
  
  # calculo del nuevo estado del sistema
  r <- transition_metropolis(s=s, node = node, delta_U = dE , beta = beta[i])
  s <- r$state
  probs_transition[i] <- r$transition_prob
  
}

plot(1:sims, system_energy, type = "l", col = "red", xlab = "epoch", ylab = "U")
plot(1:sims, magn, type = "l", pch = 19, col = "red", xlab = "epoch", ylab = "M")
plot(1:sims, probs_transition, type = "l", pch = 19, col = "black", xlab = "epoch", ylab = "P->")

plot(temps, probs_transition, type = "l", pch = 19, col = "black", xlab = "Temperature*Kb", ylab = "P->")

# # # # fin # # # # # # # # # # # # # # # # # # # # 


Resources
https://rajeshrinet.github.io/blog/2014/ising-model/
https://criticathink.wordpress.com/2018/07/15/ising-model-simulation-in-r-using-the-metropolis-monte-carlo-algorithm/
https://www.youtube.com/watch?v=kvf7aUPZCWk
https://www.hermetic.ch/compsci/tranprob.htm
http://farside.ph.utexas.edu/teaching/329/lectures/node110.html
http://phd.fisica.unimi.it/assets/Comp_Phys-esercizio3-1.pdf
https://www.asc.ohio-state.edu/braaten.1/statphys/Ising_MatLab.pdf
https://notendur.hi.is/jeg1/Ising.pdf






# # # # # ENERGY CMPUTATION OF THE SYSTEM
# entrada :
# s = estado del sistema de cada nodo 
# J = valor de los acoples
# salida
# E energia del sistema
get_energySys <- function(s, J) { 
  all_combinations <- combn(V(g),2)
  iter <- ncol(all_combinations)
  spin_prod <- numeric(length=iter)
  for (i in 1:iter) {
    id1 <- all_combinations[1,i]
    id2 <- all_combinations[2,i]
    spin_prod[i] <-  s[id1]*s[id2] 
  }
  sume <- -J*sum(spin_prod)

  return(sume)
}
# Example
# get_energySys(s=s0, J = J )
# # # # # ENERGY CMPUTATION OF THE SYSTEM




# # # # # ENERGY DIFFERENCE COMPUTATION FUNCTION
# # nota:
# El calculo de la energia se calcula solo sobre los vecinos del nodo seleccionado k.
# La diferencia de energia E(S)-E(S') siendo S=estado antes de los vecinos del nodo
# k, y S' estado despues de los vecinos de k, incluyendo nodo k. 
# Esta diferencia SE se calcula solo toma encuenta las conecciones de nodo k con
# sus vecino, porque las relaciones entre nodos vecinales se cancelan.
# O sea: DE = -J * sum sobre los vecinos de k s_k*s_i  donde i son todos los vecinos de k.

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
  s_node_flipped <- s_node*(-1)
  s_neig <- s[c(id)]
  # energia inicial
  E0 <- -J*sum(s_node*s_neig)
  E1 <- -J*sum(s_node_flipped*s_neig)
  deltaE <- (E1 - E0) 
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
      tprob <- exp(-delta_U/beta)
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
