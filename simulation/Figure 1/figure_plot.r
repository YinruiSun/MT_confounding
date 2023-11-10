
library(ggplot2)
library(gridExtra)

# type = "identity"

load("decorrelate_debias_T_type_identity.RData")

df.fdr = data.frame(p = p_seq, fdr = fdrs[,1], Method = "Decorrelate & Debias-T", graph = "Identity")
df.power = data.frame(p = p_seq, power = powers[,1], Method = "Decorrelate & Debias-T", graph = "Identity")

load("decorrelate_debias_dc_type_identity.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Decorrelate & Debias-dc", graph = "Identity"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Decorrelate & Debias-dc", graph = "Identity"))

load("doubly_debias_type_identity.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Doubly Debias", graph = "Identity"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Doubly Debias", graph = "Identity"))

load("standard_debias_type_identity.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Standard Debias", graph = "Identity"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Standard Debias", graph = "Identity"))

# type = "ER"

load("decorrelate_debias_T_type_ER.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Decorrelate & Debias-T", graph = "Random Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Decorrelate & Debias-T", graph = "Random Graph"))

load("decorrelate_debias_dc_type_ER.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Decorrelate & Debias-dc", graph = "Random Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Decorrelate & Debias-dc", graph = "Random Graph"))

load("doubly_debias_type_ER.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Doubly Debias", graph = "Random Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Doubly Debias", graph = "Random Graph"))

load("standard_debias_type_ER.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Standard Debias", graph = "Random Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Standard Debias", graph = "Random Graph"))

# type = "band"

load("decorrelate_debias_T_type_band.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Decorrelate & Debias-T", graph = "Banded Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Decorrelate & Debias-T", graph = "Banded Graph"))

load("decorrelate_debias_dc_type_band.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Decorrelate & Debias-dc", graph = "Banded Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Decorrelate & Debias-dc", graph = "Banded Graph"))

load("doubly_debias_type_band.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Doubly Debias", graph = "Banded Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Doubly Debias", graph = "Banded Graph"))

load("standard_debias_type_band.RData")

df.fdr = rbind(df.fdr, data.frame(p = p_seq, fdr = fdrs[,1], Method = "Standard Debias", graph = "Banded Graph"))
df.power = rbind(df.power, data.frame(p = p_seq, power = powers[,1], Method = "Standard Debias", graph = "Banded Graph"))

# -----------------------------

df.fdr$graph = factor(df.fdr$graph, levels = c("Identity", "Random Graph", "Banded Graph"))
df.power$graph = factor(df.power$graph, levels = c("Identity", "Random Graph", "Banded Graph"))

p.fdr = ggplot(df.fdr, aes(p, fdr, color = Method, shape = Method, linetype = Method)) +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Method, graph))) +
    facet_grid(. ~ graph) +
    labs(x = "p", y = "FDR") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
    scale_shape_manual(values = c("Decorrelate & Debias-dc" = 16, "Decorrelate & Debias-T" = 17,
                                  "Doubly Debias" = 15, "Standard Debias" = 18)) +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))

p.power = ggplot(df.power, aes(p, power, color = Method, shape = Method, linetype = Method)) +
    geom_point(size = 2.5) +
    geom_line(aes(group = interaction(Method, graph))) +
    facet_grid(. ~ graph) +
    labs(x = "p", y = "Power") +
    scale_shape_manual(values = c("Decorrelate & Debias-dc" = 16, "Decorrelate & Debias-T" = 17,
                                  "Doubly Debias" = 15, "Standard Debias" = 18)) +
    theme(legend.position = "top",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))

# figure integration
grid.arrange(p.fdr, p.power, ncol = 1)