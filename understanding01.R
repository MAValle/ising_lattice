# Del apunte "The metropolis Hastign algorithm by example"
# de John Kerl, tomo la pag.10, S4.2 para simular
# un muestro independiente.

# 07 abril, 2020
# SAMPLING WITH THE INDEPENDENCE MODEL.
# Parameters
N = 3 # numero de nodos
M = 100 # number of trials
tprob = 0.5 # umbral de probabilidad para saber si cabia de direccion el spin
q <- numeric(length=N) # vector que indica para cada estado (2^N) el numero
# de veces que ocurre.
s <- rep(-1, N) # vector de estado inicial
system_states <- matrix(NA, ncol=N, nrow=M)


for (m in 1:M) { # number of trials
  for (j in 1:N) { #numero de nodos
    u <- runif(1)
    if (u < tprob) {
      s[j] <- -1*s[j]
    } else {
      s[j] <- s[j]
    }
    
  }
  system_states[m, ] <- s
}