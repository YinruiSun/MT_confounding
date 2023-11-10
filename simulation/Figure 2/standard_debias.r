
# package loading
library(doParallel) # doparallel

# file loading
source("functions.r")

# parameter
n = 600
p = 800
s_seq = c(5, 10, 20, 30, 40, 50)
q = 5
nu = 3
mu = 3
sigma = 1
type = "identity" # "ER", "band"

alpha = 0.1

nrep = 64 * 4
ncores = 64 * 1

# parallel computing
registerDoParallel(cores = ncores)

results = list()
for (k in 1:length(s_seq))
{
    s = s_seq[k]
    print(paste("s = ", s, sep = ""))

    start = Sys.time()
    print(start)

    result = foreach(i = 1:nrep, .combine = "rbind") %dopar% {

        prob = 10 / p
        band_size = 2
        data = data_generation(n = n, p = p, s = s, q = q, sigma = sigma, type = type, nu = nu, mu = mu, prob = prob, band_size = band_size)
        result1 = test_stat_standard(Y = data$Y, X = data$X, centralize = FALSE)
        data$X = NULL
        data$Y = NULL
        result2 = multiple_testing(statistics = result1$statistics, signal = data$signal, alpha = alpha)

        list(c(data, result1, result2))
    }
    results[[k]] = result

    end = Sys.time()
    print(end)
    print(end - start)
}
rm(result)

# summary
fdp = matrix(0, nrow = length(s_seq), ncol = nrep)
power = matrix(0, nrow = length(s_seq), ncol = nrep)
fdrs = matrix(0, nrow = length(p_seq), ncol = 1)
powers = matrix(0, nrow = length(p_seq), ncol = 1)

for (k in 1:length(s_seq))
{
    fdp[k, ] = sapply(results[[k]][, 1], a <- function(x) {return(x$fdp)})
    power[k, ] = sapply(results[[k]][, 1], a <- function(x) {return(x$power)})

    fdrs[k, 1] = mean(fdp[k, ])
    powers[k, 1] = mean(power[k, ])
}

savefile = paste("standard_debias_type_", type, ".RData", sep = "")
save(n, p, s_seq, q, sigma, type, nu, mu, alpha, nrep, results, fdp, power, fdrs, powers, file = savefile)