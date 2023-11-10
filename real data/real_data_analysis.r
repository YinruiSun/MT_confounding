
# the data are available at https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/

library(ggplot2)
library(reshape2)
source("functions.r")

# data preprocess
# refer to https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/tutorial-4-r.html

data = list()
drug_names = list()

for (drug_class in c('NRTI', 'NNRTI', 'PI'))
{
    gene_file = paste(drug_class, 'GENO_PHENO.txt', sep='_')
    tsm_file = paste(drug_class, 'TSM.txt', sep='_')

    gene_df = read.delim(gene_file, na.string = c('NA', ''), stringsAsFactors = FALSE)
    tsm_df = read.delim(tsm_file, header = FALSE, stringsAsFactors = FALSE)
    names(tsm_df) = c('Position', 'Mutations')

    # Returns rows for which every column matches the given regular expression
    grepl_rows <- function(pattern, df)
    {
        cell_matches = apply(df, c(1, 2), function(x) grepl(pattern, x))
        apply(cell_matches, 1, all)
    }

    pos_start = which(names(gene_df) == 'P1')
    pos_cols = seq.int(pos_start, ncol(gene_df))
    valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[, pos_cols])
    gene_df = gene_df[valid_rows, ]

    # Flatten a matrix to a vector with names from concatenating row/column names
    flatten_matrix <- function(M, sep='.')
    {
        x <- c(M)
        names(x) <- c(outer(rownames(M), colnames(M),
                            function(...) paste(..., sep=sep)))
        x
    }

    # Construct preliminary design matrix
    muts = c(LETTERS, 'i', 'd')
    X = outer(muts, as.matrix(gene_df[, pos_cols]), Vectorize(grepl))
    X = aperm(X, c(2, 3, 1))
    dimnames(X)[[3]] <- muts
    X = t(apply(X, 1, flatten_matrix))
    mode(X) = 'numeric'

    # Remove any mutation/position pairs that never appear in the data
    X = X[, colSums(X) != 0]

    # Remove duplicate predictors
    X = X[, !duplicated(t(X))]

    # Extract response matrix
    Y = gene_df[, 4:(pos_start - 1)]

    data[[drug_class]] = list(Y = Y, X = X, TSM = tsm_df)
    drug = colnames(Y)
    drug_names[[drug_class]] = paste(drug, ' (', drug_class, ')', sep = '')
}

save(data, drug_names, file = 'HIV_data.RData')

load('HIV_data.RData')

# svd

k = 100 # top k singular values
svd_d = data.frame()
sizes = data.frame()

for (drug_class in c('NRTI', 'NNRTI'))
{
    for (i in 1:length(drug_names[[drug_class]]))
    {
        Y = data[[drug_class]]$Y[, i]
        X = data[[drug_class]]$X
        drug = drug_names[[drug_class]][i]

        # Remove patients with missing measurements
        missing = is.na(Y)
        X = X[!missing, ]

        ## Remove any mutation/position pairs that never appear in the data
        X = X[, colSums(X) != 0]

        # Remove duplicate predictors
        X = X[, !duplicated(t(X))]

        X = scale(X)

        # svd
        temp = RSpectra::svds(X, k = k, nu = 0, nv = 0)$d
        df = data.frame(Order = 1:k, SV = temp, drug = drug)
        svd_d = rbind(svd_d, df)

        # size
        sizes = rbind(sizes, data.frame(drug = drug, sample_size = dim(X)[1],
                                        dimension = dim(X)[2]))
    }
}

# plot

svd_d$drug = factor(svd_d$drug, ordered = T,
                    levels = unlist(drug_names[c('NRTI', 'NNRTI')]))
sizes$drug = factor(sizes$drug, ordered = T,
                    levels = unlist(drug_names[c('NRTI', 'NNRTI')]))
sizes$Text = paste('n = ', sizes$sample_size, ', p = ', sizes$dimension, sep = '')

p_svd = ggplot(svd_d, aes(x = Order, y = SV, group = drug)) +
    geom_point(size = 0.5) +
    labs(title = "", x = "", y = paste("Top", k, "Singular Values")) +
    facet_wrap(~ drug, ncol = 3) +
    theme(strip.text = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 10)) +
    geom_text(data = sizes, aes(label = Text), x = 50, y = 70, size = 5)

p_svd

# testing procedures

testing_procedures = function(Y, X, q, alpha)
{
    # Log-transform the drug resistance measurements
    Y = log(Y)

    # Remove patients with missing measurements
    missing = is.na(Y)
    Y = Y[!missing]
    X = X[!missing, ]

    # Remove duplicate predictors
    X = X[, !duplicated(t(X))]

    ## Remove any mutation/position pairs that never appear in the data
    X = X[, colSums(X) != 0]

    # standardization
    Y = scale(Y)
    X = scale(X)

    # decorrelating and debiasing method
    stat = test_stat_decorrelation(Y, X, initial = "trim", q = q)$statistics
    Decorrelate_Debias_T = multiple_testing(statistics = stat, alpha = alpha)$index_discovery
    Decorrelate_Debias_T = colnames(X)[Decorrelate_Debias_T]

    stat = test_stat_decorrelation(Y, X, initial = "decorrelation", q = q)$statistics
    Decorrelate_Debias_dc = multiple_testing(statistics = stat, alpha = alpha)$index_discovery
    Decorrelate_Debias_dc = colnames(X)[Decorrelate_Debias_dc]

    # doubly debiased method in Guo, Cevid and Buhlmann (2022) AOS
    stat = test_stat_doubly(Y, X)$statistics
    Doubly_Debias = multiple_testing(statistics = stat, alpha = alpha)$index_discovery
    Doubly_Debias = colnames(X)[Doubly_Debias]

    # standard debiased method in Zhang and Zhang (2014) JRSSB
    stat = test_stat_standard(Y, X)$statistics
    Standard_Debias = multiple_testing(statistics = stat, alpha = alpha)$index_discovery
    Standard_Debias = colnames(X)[Standard_Debias]

    # result
    discovery = list(Decorrelate_Debias_dc, Decorrelate_Debias_T,
                     Doubly_Debias, Standard_Debias)
    names(discovery) = c('Decorrelate & Debias-dc', 'Decorrelate & Debias-T',
                         'Doubly Debias', 'Standard Debias')
    return(list(discovery = discovery, size = dim(X)))
}

q = c(1, 1)
names(q) = c('NRTI', 'NNRTI')
alpha = 0.1

results = list()
for (drug_class in c('NRTI', 'NNRTI'))
{
    Y = data[[drug_class]]$Y
    X = data[[drug_class]]$X

    results[[drug_class]] = sapply(Y, function(y) testing_procedures(y, X, q[drug_class], alpha))
    colnames(results[[drug_class]]) = drug_names[[drug_class]]
}
# in the paper, we fix the randomness in cross validation steps for replication
# to replicate the results, set the following `foldid` for the argument `foldid` in `cv.glmnet` in initial estimation steps:
# n = dim(X)[1]
# num_group = 10
# foldid = rep(1:num_group, each = n %/% num_group)
# remainder = n %% num_group
# foldid = c(foldid, rep(num_group, remainder))

file = paste('HIV_results_alpha', alpha, '.RData', sep = '')
save(results, alpha, file = file)

load(file)

# result evaluation

get_position <- function(x)
{
    sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)
}

discoveries = data.frame()

for (drug_class in c('NRTI', 'NNRTI'))
{
    tsm = data[[drug_class]]$TSM
    res = results[[drug_class]]

    for (drug in drug_names[[drug_class]])
    {
        num_discovery = sapply(res['discovery', ][[drug]], function(discovery){
            positions = unique(get_position(discovery))
            total_discovery = length(positions)
            false_discovery = length(setdiff(positions, tsm$Position))
            true_discovery = total_discovery - false_discovery
            temp = c(true_discovery, false_discovery)
            names(temp) = c('true_discovery', 'false_discovery')
            return(temp)
        })

        discoveries_temp = data.frame(true_discovery = num_discovery[1, ],
                                      false_discovery = num_discovery[2, ],
                                      Method = colnames(num_discovery),
                                      drug = drug)
        discoveries = rbind(discoveries, discoveries_temp)
    }
}

# plot

discoveries_melt = melt(discoveries, id.vars = c('Method', 'drug'),
                        variable.name = "Discovery", value.name = "Num_Discovery")
discoveries_melt$drug = factor(discoveries_melt$drug, ordered = T,
                          levels = unlist(drug_names[c('NRTI', 'NNRTI')]))

p_bar = ggplot(discoveries_melt, aes(x = Method, y = Num_Discovery,
                                     fill = Discovery, group = drug)) +
    geom_bar(stat = "identity") +
    labs(title = "", x = "", y = "Number of Discoveries") +
    facet_wrap(~ drug, ncol = 3) +
    ylim(0, 50) +
    scale_fill_manual(values = c("true_discovery" = "#4169E1", "false_discovery" = "#F8766D"), 
                      labels = c("In TSM Set", "Not In TSM Set")) +
    theme(legend.position = "top",
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 13),
          strip.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
          axis.title.y = element_text(size = 15))