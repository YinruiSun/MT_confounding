
# codes for A DECORRELATING AND DEBIASING APPROACH TO SIMULTANEOUS INFERENCE FOR HIGH-DIMENSIONAL CONFOUNDED MODELS

# ---------- package loading ----------

library(MASS) # sampling from multivariate normal distribution
library(RSpectra) # svd/eigen decomposition
library(Matrix) # sparse matrix
library(glmnet) # lasso fit

# ---------- data generation ----------

data_generation = function(n, p, s, type, nu = 3, mu = 3, q = 5, sigma = 1,
                           prob = 0.01, band_size = 2, con_num = 50)
{
    # n: sample size
    # p: dimension of covariance / precision matrix
    # s: sparsity of regression coefficient
    # type: "identity", "ER" or "band"
    # nu: parameter for coefficient signal strength
    # mu: parameter for vector phi
    # q: number of confounders
    # sigma: sd of noise
    # prob: scalar in (0,1), probability of random graph, supply if type == "ER"
    # band_size: integer, band width of band graph, supply if type == "band"
    # con_num: scalar, supply if type == "ER" or "band"

    cov_result = cov_function(p, type, prob, band_size, con_num)
    Sigma = cov_result$Sigma
    Omega_diag = cov_result$Omega_diag
    rm(cov_result)

    signal = sample(1:p, s, replace = FALSE)
    beta = rep(0, p)
    beta[signal] = 1.2^(-nu) * sigma * sqrt(8 * Omega_diag[signal] * log(p) / n)
    beta[signal] = sample(c(-1, 1), s, replace = TRUE) * beta[signal]
    beta = Matrix::Matrix(beta, ncol = 1, sparse = TRUE)
    signal = sort(signal)

    Gamma = matrix(runif(q * p, min = -2, max = 2), nrow = q, ncol = p) # iid uniform entry
    phi = rnorm(q, mean = mu, sd = 1)
    H = matrix(rnorm(n * q, mean = 0, sd = 1), nrow = n, ncol = q)
    E = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    xi = rnorm(n, mean = 0, sd = sigma)
    rm(Sigma)

    X = H %*% Gamma + E
    Y = X %*% beta + H %*% phi + xi
    Y = as.vector(Y)
    rm(Gamma, phi, H, E, xi)

    result = list(Y = Y, X = X, beta = beta, signal = signal, Omega_diag = Omega_diag)
    return(result)
}

cov_function = function(p, type, prob, band_size, con_num)
{
    # p: dimension of covariance / precision matrix
    # type: "identity", "ER" or "band"
    # prob: scalar in (0,1), probability of random graph, supply if type == "ER"
    # band_size: integer, band width of band graph, supply if type == "band"
    # con_num: scalar, supply if type == "ER" or "band"

    if (type == "identity")
    {
        Sigma = diag(1, p)
    }else
    {
        if (type == "ER") Omega = matrix(rbinom(p^2, 1, prob), ncol = p, nrow = p)

        if (type == "band")
        {
            Omega = matrix(0, ncol = p, nrow = p)
            for (i in 1:band_size)
            {
                diag(Omega[1:(p - i), (i + 1):p]) = 1 / i
            }
        }

        Omega = Omega * matrix(runif(p^2, 0.5, 1), ncol = p, nrow = p)
        Omega = Omega * (matrix(rbinom(p^2, 1, 0.5), ncol = p, nrow = p) * 2 - 1)
        diag(Omega) = 0
        Omega = as.matrix(forceSymmetric(Omega))

        sigma_min = RSpectra::eigs_sym(Omega, k = 1, which = "SA", opts = list(retvec = FALSE))$values
        sigma_max = RSpectra::eigs_sym(Omega, k = 1, which = "LA", opts = list(retvec = FALSE))$values
        u = (sigma_max - sigma_min) / (con_num - 1)
        diag(Omega) = abs(sigma_min) + u
        Omega = Omega / (abs(sigma_min) + u)

        Sigma = solve(Omega)
        rm(Omega)
    }

    return(list(Sigma = Sigma, Omega_diag = rep(1, p)))
}

# ---------- multiple testing ----------

multiple_testing = function(statistics, signal = NULL, alpha = 0.1, num_discrete = 100)
{
    # statistics: vector of studentized statistics
    # alpha: scalar, prespecified error rate
    # signal: vector of locations of nonzero signals. Supplied in simulations.
    # num_discrete: scalar.

    p = length(statistics)
    t_p = sqrt(2 * log(p) - 2 * log(log(p)))
    # discrete points in interval [0,t_p]
    t = t_p * (1:num_discrete) / num_discrete

    fdp = rep(0, num_discrete) # estimated fdp
    for (i in 1:num_discrete)
    {
       num_discovery = sum(abs(statistics) >= t[i])
       fdp[i] = 2 * p * (1 - pnorm(t[i])) / max(num_discovery, 1)
    }

    if (min(fdp) <= alpha)
    {
        threshold = t[min(which(fdp <= alpha))]
    }else
    {
        threshold = sqrt(2 * log(p))
    }

    index_discovery = which(abs(statistics) >= threshold)
    num_discovery = length(index_discovery)

    fdp = NULL
    power = NULL
    index_Tdiscovery = NULL
    num_Tdiscovery = NULL

    if (is.null(signal) == 0)
    {
        index_Tdiscovery = intersect(index_discovery, signal)
        num_Tdiscovery = length(index_Tdiscovery)

        power = num_Tdiscovery / length(signal)
        fdp = (num_discovery - num_Tdiscovery) / max(num_discovery, 1)
    }

    result = list(fdp = fdp, power = power, alpha = alpha, threshold = threshold,
                  num_discovery = num_discovery, index_discovery = index_discovery,
                  num_Tdiscovery = num_Tdiscovery, index_Tdiscovery = index_Tdiscovery)
    return(result)
}

# ---------- nodewise lasso regression (for decorrelate & debias and standard debias) ----------

nodewise_lasso = function(X, nlambda = 100)
{
    # X: matrix, size = n * p
    # nlambda: scalar, number of candidate parameters

    p = ncol(X)
    n = nrow(X)

    # nodewise lasso procedure
    Gamma = NULL
    denominator = NULL
    sequence = (nlambda:1) / nlambda * 3 * sqrt(log(p) / n)
    for (j in 1:p)
    {
        lambda_j = sqrt(sum(X[, j]^2) / n) * sequence
        fit_j = glmnet::glmnet(x = X[, -j], y = X[, j], family = "gaussian", lambda = lambda_j, intercept = FALSE)

        Gamma = rbind(Gamma, fit_j$beta) # matrix, size = (p - 1) * nlambda
        temp = as.vector(crossprod(X[, j] - predict(fit_j, newx = X[, -j]), X[, j]))
        denominator = rbind(denominator, temp)  # vector, length = nlambda
    }
    rm(lambda_j, fit_j, temp)

    # tuning parameter selection
    # lambda: from large to small
    index = 1
    IC = Inf

    for (i in 1:nlambda)
    {
        Omega = matrix(0, ncol = p, nrow = p)
        gamma = matrix(Gamma[, i], nrow = (p - 1), ncol = p)

        Omega[upper.tri(Omega, diag = FALSE)] = - gamma[upper.tri(gamma, diag = FALSE)]
        Omega[lower.tri(Omega, diag = FALSE)] = - gamma[lower.tri(gamma, diag = TRUE)]
        diag(Omega) = 1
        rm(gamma)

        Omega = t(t(Omega) / (denominator[, i] / n))

        # symmetrization
        Omega = Omega * (abs(Omega) <= abs(t(Omega))) + t(Omega) * (abs(Omega) > abs(t(Omega)))
        Omega = Matrix::Matrix(Omega, sparse = TRUE)

        # IC
        eig_min = RSpectra::eigs_sym(as(Omega, "matrix"), k = 1, which = "SA", opts = list(retvec = FALSE))$values
        if (eig_min <= 0)
        {
            next
        }else
        {
            IC_new = - determinant(Omega)$modulus[1] + log(n) / n * (sum(Omega != 0) - p) / 2 +
                     sum(Omega * crossprod(X) / n)
            if (IC_new <= IC)
            {
                IC = IC_new
                index = i
            }
        }
    }
    gamma = Matrix::Matrix(Gamma[, index], nrow = (p - 1), ncol = p, sparse = TRUE)
    denominator = denominator[, index]
    rm(Gamma, Omega)

    return(list(gamma = gamma, denominator = denominator))
}

# ---------- method: decorrelate & debias ----------

test_stat_decorrelation = function(Y, X, initial = "trim", centralize = TRUE,
                                   q = NULL, nlambda = 100, rho = 0.3, kmax = 20)
{
    # Y: vector, length = p; or matrix, size = p * 1
    # X: matrix, size = n * p
    # initial: "trim" or "decorrelation"
    # centralize: TRUE or FALSE, Y and X would be centralized if centralize = TRUE
    #       in simulations data are generated with zero population mean and centralize = FALSE
    # q: scalar, number of factors. If q = NULL, q would be determined by eigen_ratio method
    #       in practical data analysis, a larger q is preferred
    # nlambda: scalar, number of candidates parameters in nodewise lasso procedures
    # rho: tuning parameter in intial trim estimation procedure
    # kmax: max number of factors in eigen_ratio method

    # sample size and dimension
    n = dim(X)[1]
    p = dim(X)[2]

    Y = scale(Y, centralize, FALSE)
    X = scale(X, centralize, FALSE)

    # svd
    num_sv = max(c(floor(rho * min(n, p)), kmax + 1, q))

    UDV_list = RSpectra::svds(A = X, k =  num_sv, nu = num_sv, nv = 0)
    U_truncated = UDV_list$u # matrix, singular vectors
    D_truncated = UDV_list$d # vector, left singular values

    # construction of transform matrix
    if (initial == "trim")
    {
        # construction of trim transform matrix
        D_tilde =  D_truncated[floor(rho * min(n, p))] / D_truncated[1:(floor(rho * min(n, p)))]
        D_tilde = sqrt(1 - D_tilde)
        F_trim = U_truncated[, 1:(floor(rho * min(n, p)))] %*% diag(D_tilde)
        F_trim = diag(n) - F_trim %*% t(F_trim)
    }

    # construction of decorrelation transform matrix
    if (is.null(q))
    {
        # eigenvalue ratio method to estimate q
        D_1_kmax = D_truncated[1:kmax]
        D_2_kmaxp1 = D_truncated[2:(kmax + 1)]
        q = which.max(D_1_kmax / D_2_kmaxp1)
    }
    F_dc = diag(n) - U_truncated[, 1:q] %*% t(U_truncated[, 1:q])

    rm(UDV_list, U_truncated, D_truncated)

    # initial estimator for coefficient & estimator for variance of noise
    F1 = F_dc
    if (initial == "trim") F1 = F_trim
    F1_Y = F1 %*% Y
    F1_X = F1 %*% X

    sequence = (nlambda:1) / nlambda * 3
    lambda = sqrt(sum(F1_Y^2) / n) * sqrt(log(p) / n) * sequence
    fit = glmnet::cv.glmnet(x = F1_X, y = F1_Y, lambda = lambda, intercept = FALSE, family = "gaussian", nfolds = 10)
    beta_init = coef(fit, s = fit$lambda.min)[-1] # vector, length = p
    beta_init = Matrix::Matrix(beta_init, ncol = 1, sparse = TRUE) # matrix, size = p * 1

    signal_hat = which(beta_init != 0)
    sigma_hat = sum((F1_Y - F1_X %*% beta_init)^2) # scalar
    if (length(signal_hat) == 0)
    {
        sigma_hat = sigma_hat / sum(F1 * F1)
        sigma_hat = sqrt(sigma_hat)
    }else {
        UDV_list = svd(F1_X[, signal_hat], nv = 0)
        rank = sum(UDV_list$d != 0)
        projection = tcrossprod(UDV_list$u[, 1:rank])
        projection = (diag(1, n) - projection) %*% F1
        sigma_hat = sigma_hat / sum(projection * projection)
        sigma_hat = sqrt(sigma_hat)
    }
    rm(F1, F1_X, F1_Y, fit)
    if (length(signal_hat) != 0) rm(signal_hat, UDV_list, rank, projection)
    if (initial == "trim") rm(F_trim)

    # debiased estimator for coefficient
    statistics = rep(0, p)
    beta_debias = rep(0, p)

    F2_Y = F_dc %*% Y
    F2_X = F_dc %*% X
    rm(F_dc)
    rm(X, Y)

    nodewise_fit = nodewise_lasso(F2_X, nlambda)
    gamma = nodewise_fit$gamma
    denominator = nodewise_fit$denominator
    rm(nodewise_fit)

    for (j in 1:p)
    {
        z_j = as.vector(F2_X[, j] - F2_X[, -j] %*% gamma[, j]) # vector, length = n
        z_j_norm = sqrt(sum(z_j^2))

        # debiased estimators & studentized statistics
        beta_debias[j] = beta_init[j] + sum(z_j * (F2_Y - F2_X %*% beta_init)) / denominator[j]
        statistics[j] = beta_debias[j] * z_j_norm / sigma_hat
    }
    rm(F2_X, F2_Y, z_j, gamma, denominator)

    result = list(statistics = statistics, beta_debias = beta_debias, beta_init = beta_init, 
                  q = q, sigma_hat = sigma_hat)
    return(result)
}

# ---------- method: doubly debias ----------

test_stat_doubly = function(Y, X, centralize = TRUE, nlambda = 100, rho = 0.3, rho_j = 0.3)
{
    # Y: vector, length = p; or matrix, size = p * 1
    # X: matrix, size = n * p
    # centralize: TRUE or FALSE, Y and X would be centralized if centralize = TRUE
    #       (in simulations data are generated with zero population mean and centralize = FALSE)
    # nlambda: scalar, number of candidates parameters in lasso procedures
    # rho: tuning parameter in intial trim estimation procedure
    # rho_j: tuning parameter in debiased procedure

    # sample size and dimension
    n = dim(X)[1]
    p = dim(X)[2]

    Y = scale(Y, centralize, FALSE)
    X = scale(X, centralize, FALSE)

    # svd
    num_sv = floor(rho * min(n, p))
    UDV_list = RSpectra::svds(A = X, k =  num_sv, nu = num_sv, nv = 0)
    U_truncated = UDV_list$u # matrix, singular vectors
    D_truncated = UDV_list$d # vector, left singular values

    # construction of trim transform matrix
    D_tilde =  D_truncated[floor(rho * min(n, p))] / D_truncated[1:(floor(rho * min(n, p)))]
    D_tilde = sqrt(1 - D_tilde)
    F_trim = U_truncated[, 1:(floor(rho * min(n, p)))] %*% diag(D_tilde)
    F_trim = diag(n) - F_trim %*% t(F_trim)

    rm(UDV_list, U_truncated, D_truncated)

    # initial estimator for coefficient & estimator for variance of noise
    F1_Y = F_trim %*% Y
    F1_X = F_trim %*% X

    sequence = (nlambda:1) / nlambda * 3
    lambda = sqrt(sum(F1_Y^2) / n) * sqrt(log(p) / n) * sequence
    fit = glmnet::cv.glmnet(x = F1_X, y = F1_Y, lambda = lambda, intercept = FALSE, family = "gaussian", nfolds = 10)
    beta_init = coef(fit, s = fit$lambda.min)[-1] # vector, length = p
    beta_init = Matrix::Matrix(beta_init, ncol = 1, sparse = TRUE) # matrix, size = p * 1

    signal_hat = which(beta_init != 0)
    sigma_hat = sum((F1_Y - F1_X %*% beta_init)^2) # scalar
    if (length(signal_hat) == 0)
    {
        sigma_hat = sigma_hat / sum(F_trim * F_trim)
        sigma_hat = sqrt(sigma_hat)
    }else {
        UDV_list = svd(F1_X[, signal_hat], nv = 0)
        rank = sum(UDV_list$d != 0)
        projection = tcrossprod(UDV_list$u[, 1:rank])
        projection = (diag(1, n) - projection) %*% F_trim
        sigma_hat = sigma_hat / sum(projection * projection)
        sigma_hat = sqrt(sigma_hat)
    }
    rm(F_trim, F1_X, F1_Y, fit)
    if (length(signal_hat) != 0) rm(signal_hat, UDV_list, rank, projection)

    # debiased estimator for coefficient
    statistics = rep(0, p)
    beta_debias = rep(0, p)

    for (j in 1:p)
    {
        # svd
        num_sv = floor(rho_j * min(n, p - 1))
        UDV_list_j = RSpectra::svds(A = X[, -j], k = num_sv, nu = num_sv, nv = 0)
        U_truncated_j = UDV_list_j$u # matrix
        D_truncated_j = UDV_list_j$d # vector

        # trim transformation for X_negj
        D_tilde_j =  D_truncated_j[num_sv] / D_truncated_j[1:num_sv]
        D_tilde_j = sqrt(1 - D_tilde_j)
        F_trim_j = U_truncated_j %*% diag(D_tilde_j)
        F_trim_j = diag(n) - F_trim_j %*% t(F_trim_j)

        # preconditioned repsonse and design
        Fj_Y = F_trim_j %*% Y
        Fj_X = F_trim_j %*% X

        # lasso regression
        lambda_j = sqrt(sum(Fj_X[, j]^2) / n) * sqrt(log(p) / n) * sequence
        fit_j = glmnet::cv.glmnet(x = Fj_X[, -j], y = Fj_X[, j], lambda = lambda_j, intercept = FALSE, family = "gaussian", nfolds = 5)

        gamma_j = coef(fit_j, s = fit_j$lambda.min)[-1]
        gamma_j = Matrix::Matrix(gamma_j, ncol = 1, sparse = TRUE)
        z_j = as.vector(Fj_X[, j] - Fj_X[, -j] %*% gamma_j)
        denominator_j = sum(z_j * Fj_X[, j])
        Fj_zj_norm = sqrt(sum((F_trim_j %*% z_j)^2))
        V_inflate = sqrt(1.25) * Fj_zj_norm / denominator_j

        # tuning parameter selection
        lambda_seq = sort(fit_j$lambda) # increasing order
        for (lam in lambda_seq)
        {
            gamma_j = coef(fit_j, s = lam)[-1]
            gamma_j = Matrix(gamma_j, ncol = 1, sparse = TRUE)
            z_j = as.vector(Fj_X[, j] - Fj_X[, -j] %*% gamma_j)
            denominator_j = sum(z_j * Fj_X[, j])
            Fj_zj_norm = sqrt(sum((F_trim_j %*% z_j)^2))
            tau_j = Fj_zj_norm / denominator_j
            if (tau_j <= V_inflate) break
        }

        # debiased estimators & studentized statistics
        beta_debias[j] = beta_init[j, 1] + sum(z_j * (Fj_Y - Fj_X %*% beta_init)) / denominator_j
        statistics[j] = beta_debias[j] / (sigma_hat * tau_j)
    }
    rm(UDV_list_j, U_truncated_j, D_truncated_j, D_tilde_j, F_trim_j,
       Fj_X, Fj_Y, fit_j, gamma_j, z_j, denominator_j, Fj_zj_norm)
    rm(X, Y)

    result = list(statistics = statistics, beta_debias = beta_debias, beta_init = beta_init, sigma_hat = sigma_hat)
    return(result)
}

# ---------- method: standard debias ----------

test_stat_standard = function(Y, X, centralize = TRUE, nlambda = 100)
{
    # Y: vector, length = p; or matrix, size = p * 1
    # X: matrix, size = n * p
    # centralize: TRUE or FALSE, Y and X would be centralized if centralize = TRUE
    #       (in simulations data are generated with zero population mean and centralize = FALSE)
    # nlambda: scalar, number of candidates parameters in lasso procedures

    # sample size and dimension
    n = dim(X)[1]
    p = dim(X)[2]

    Y = scale(Y, centralize, FALSE)
    X = scale(X, centralize, FALSE)

    # initial estimator for coefficient & estimator for variance of noise
    sequence = (nlambda:1) / nlambda * 3
    lambda = sqrt(sum(Y^2) / n) * sqrt(log(p) / n) * sequence
    fit = glmnet::cv.glmnet(x = X, y = Y, lambda = lambda, intercept = FALSE, family = "gaussian", nfolds = 10)
    beta_init = coef(fit, s = fit$lambda.min)[-1] # vector, length = p
    beta_init = Matrix::Matrix(beta_init, ncol = 1, sparse = TRUE) # matrix, size = p * 1

    signal_hat = which(beta_init != 0)
    sigma_hat = sum((Y - X %*% beta_init)^2) # scalar
    if (length(signal_hat) == 0)
    {
        sigma_hat = sigma_hat / n
        sigma_hat = sqrt(sigma_hat)
    }else {
        UDV_list = svd(X[, signal_hat], nv = 0)
        rank = sum(UDV_list$d != 0)
        sigma_hat = sigma_hat / (n - rank)
        sigma_hat = sqrt(sigma_hat)
    }
    rm(fit)
    if (length(signal_hat) != 0) rm(signal_hat, UDV_list, rank)

    # debiased estimator for coefficient
    statistics = rep(0, p)
    beta_debias = rep(0, p)

    nodewise_fit = nodewise_lasso(X, nlambda)
    gamma = nodewise_fit$gamma
    denominator = nodewise_fit$denominator
    rm(nodewise_fit)

    for (j in 1:p)
    {
        z_j = as.vector(X[, j] - X[, -j] %*% gamma[, j]) # vector, length = n
        z_j_norm = sqrt(sum(z_j^2))

        # debiased estimators & studentized statistics
        beta_debias[j] = beta_init[j] + sum(z_j * (Y - X %*% beta_init)) / denominator[j]
        statistics[j] = beta_debias[j] * denominator[j] / (sigma_hat * z_j_norm)
    }
    rm(X, Y, gamma, denominator, z_j, z_j_norm)

    result = list(statistics = statistics, beta_debias = beta_debias, beta_init = beta_init, sigma_hat = sigma_hat)
    return(result)
}
