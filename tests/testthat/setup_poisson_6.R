set.seed(123)
N = 5
G = 60

P <- read.csv(
    "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt",
    sep = '\t'
)
P <- as.matrix(P[,2:ncol(P)])

# change to 6x6 only center mutation. collapse the cosmic matrix
P6 <- matrix(nrow = 6, ncol = ncol(P))
start <- 1
end <- 16
for (i in (1:6)){
    P6[i,] = as.matrix(colSums(P[start:end,]))
    start <- end + 1
    end <- start + 15
}

# DELTED BELOW (low cossim)
sim <- pairwise_sim(P6, P6)
row_maxes <- sapply(1:nrow(sim), function(i) {
    sum(sim[i,-i] > 0.7)
})

P6 <- P6[,-which(row_maxes > 0)]


# choosing 5 random signatures
sigs = 2:(2+N - 1) # c(2,3,4,5,10)
P6 <- P6[,sigs]

K = nrow(P6)

E <- matrix(rexp(N*G, 0.001), nrow = N, ncol = G)

M <- matrix(nrow = K, ncol = G)
for (k in 1:K) {
    M[k,] <- rpois(G, P6[k,]%*%E)
}
M <- round(M)

true_P <- P6
true_E <- E

# remove these variables from the environment after source
rm(list = c('P6','E','N','G','sigs'))
