
# package loading
library(doParallel) # doparallel

# file loading
source('functions.r')

# parameter
n = 600
p_seq = c(4, 6, 8, 10, 12) * 100
s = 30
q = 5
nu = 3
mu = 3
sigma = 1
type = "identity" # "ER", "band"

initial = "decorrelation"
rho = 0.3
alpha = 0.1

nrep = 64 * 4
ncores = 64 * 1

# parallel computing
registerDoParallel(cores = ncores)

results = list()
for (k in 1:length(p_seq))
{
    p = p_seq[k]
    print(paste("p = ", p, sep = ""))

    start = Sys.time()
    print(start)

    result = foreach(i = 1:nrep, .combine = "rbind") %dopar% {

        prob = 10 / p
        band_size = 2
        data = data_generation(n = n, p = p, s = s, q = q, sigma = sigma, type = type, nu = nu, mu = mu, prob = prob, band_size = band_size)
        result1 = test_stat_decorrelation(Y = data$Y, X = data$X, initial = initial, rho = rho, centralize = FALSE)
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
fdp = matrix(0, nrow = length(p_seq), ncol = nrep)
power = matrix(0, nrow = length(p_seq), ncol = nrep)
fdrs = matrix(0, nrow = length(p_seq), ncol = 1)
powers = matrix(0, nrow = length(p_seq), ncol = 1)

for (k in 1:length(p_seq))
{
    fdp[k, ] = sapply(results[[k]][, 1], a <- function(x) {return(x$fdp)})
    power[k, ] = sapply(results[[k]][, 1], a <- function(x) {return(x$power)})

    fdrs[k, 1] = mean(fdp[k, ])
    powers[k, 1] = mean(power[k, ])
}

savefile = paste("decorrelate_debias_dc_type_", type, ".RData", sep = "")
save(n, p_seq, s, q, sigma, type, nu, mu, initial, rho, alpha, nrep, results, fdp, power, fdrs, powers, file = savefile)