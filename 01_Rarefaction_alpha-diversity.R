###############################################################
# Script: 01_rarefaction_alpha_diversity.R
#
# Purpose:
# Assess sequencing depth, generate rarefaction curves,
# rarefy ASV counts, calculate alpha diversity metrics,
# and compare diversity among experimental groups.
#
# Input:
#   data/Abundance_allGroup.csv
#   data/Sample_metadata_all_group.csv
#
# Output:
#   results/Sequencing_Depth_Distribution.tiff
#   results/Rarefaction_Curves.tiff
#   results/Alpha_Diversity.tsv
#   results/Alpha_Diversity.tiff
#
# Author:Moirangthem Goutam Singh, Romi Wahengbam
# Date:9-7-2026
###############################################################

#Load packages
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggdist)

abundance_file <- "data/Abundance_allGroup.csv"
metadata_file  <- "data/Sample_metadata_all_group.csv"

output_dir <- "results"

chosen_depth <- 30000

set.seed(123)

# Read data
asv_table <- read.csv(abundance_file, row.names = 1, check.names = FALSE )

metadata <- read.csv(metadata_file, row.names = 1)

# Remove low-abundnace features
asv_table <- asv_df[rowSums(asv_df) > 0, ]

asv_matrix <- t(asv_table)

# Examine sequencing depth
lib_size <- rowSums(asv_mat)

summary(lib_size)

quantile(lib_size, c(0.05, 0.1, 0.25))

hist(lib_size,
     breaks = 30,
     col = "grey",
     main = "Sequencing depth per sample",
     xlab = "Reads per sample")

df_depth <- data.frame(lib_size = lib_size)
#Plot the sequencing dept
p_depth <- ggplot(df_depth, aes(x = lib_size, fill = after_stat(count))) +
  geom_histogram(bins = 30, color = "black") +
  
  scale_fill_gradientn(
    colours = c("#0D0887", "#7E03A8", "#CC4778", "#F89441", "#F0F921")
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    axis.title.x = element_text(size = 26),
    axis.text = element_text(size = 26),
    axis.title.y = element_text(size = 26)
  ) +
  labs(
    title = "Sequencing Depth Distribution",
    x = "Reads per sample",
    y = "Number of samples",
    fill = "Sample count"
  )
#Save the plot
ggsave(file.path(output_dir,"Sequencing_Depth_Distribution.tiff"),
       plot = p_depth,
       width = 10,
       height = 6, units = "in",
       dpi = 900, compression ="lzw")


# Generate rarefaction curves
set.seed(123)

rarecurve_out <- rarecurve(
  asv_matrix,
  step = 500,
  label = FALSE
)

step_size <- 500
sample_ids <- rownames(asv_matrix)

rare_df <- do.call(rbind, lapply(seq_along(rarecurve_out), function(i) {
  
  sp <- rarecurve_out[[i]]
  if (length(sp) == 0) return(NULL)
  
  data.frame(
    Sample = sample_ids[i],
    Reads = seq_len(length(sp)) * step_size,
    Observed_ASVs = as.numeric(sp)
  )
}))

#plot the rarefraction
p_rare <- ggplot(rare_df,
                 aes(x = Reads,
                     y = Observed_ASVs,
                     group = Sample)) +
  
  # Individual sample curves
  geom_line(color = "skyblue", alpha = 0.4) +
  
  # Overall smooth trend
  geom_smooth(aes(group = 1),
              color = "#D55E00",
              size = 1.2,
              se = FALSE) +
  
  # Threshold line
  geom_vline(xintercept = chosen_depth,
             linetype = "dashed",
             color = "red",
             size = 1) +
  
  # Annotation
  annotate("text",
           x = chosen_depth,
           y = max(rare_df$Observed_ASVs),
           label = paste("Depth =", chosen_depth),
           hjust = -0.1,
           size = 8) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 26),
    axis.text = element_text(size = 24),
    axis.title.y = element_text(size = 26)
  ) +
  
  labs(
    title = "Rarefaction Curves with Selected Depth",
    x = "Sequencing depth (reads)",
    y = "Observed ASVs"
  )

#Save the plot
ggsave(file.path(output_dir, "Rarefaction_Curves_cleanall_group.tiff"),
       p_rare,
       width = 8,
       height = 6,
       units = "in",
       dpi = 900, compression= "lzw")

#Rarefy the feature data
set.seed(123)

asv_rarefied <- rrarefy(
  asv_matrix,
  sample = chosen_depth
)

#Remove fail data
asv_rarefied <- asv_rarefied[rowSums(asv_rarefied) > 0, ]


#Compute alpha diversity
alpha_div <- data.frame(
  Sample = rownames(asv_rarefied),
  Observed = rowSums(asv_rarefied > 0),
  Shannon = diversity(asv_rarefied, index = "shannon"),
  Simpson = diversity(asv_rarefied, index = "simpson"),
  InvSimpson = diversity(asv_rarefied, index = "invsimpson"),
  Pielou = diversity(asv_rarefied, index = "shannon") /
    log(specnumber(asv_rarefied))
)

#Long format
alpha_long <- alpha_div %>%
  left_join(metadata %>% rownames_to_column("Sample"),
            by = "Sample") %>%   # adjust SampleID
  pivot_longer(
    cols = c(Observed, Shannon, Simpson, InvSimpson, Pielou),
    names_to = "Index",
    values_to = "Value"
  )


# Plot the alpha diversity

p1 <- ggplot(alpha_long,
             aes(x = Diet.x, y = Value, fill = Diet.x)) +
  stat_halfeye(
    adjust = 0.6,
    width = 0.8,
    justification = -0.2,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    width = 0.3,
    alpha = 0.5,
    size = 3,
    color = "black"
  ) +
  facet_wrap(~ Index, scales = "free_y", nrow = 2, ncol = 3) +
  scale_fill_manual(
    values = c( "Diet1" = "blue",  # Here we are using 4 different types of diet, can be customise based on the user's requirement
                "Diet2" = "green",
                "Diet3" = "orange",
                "Diet4" = "red"
    )
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.title = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "plain",size = 23,
      margin = margin(b = 10)
    ),
    panel.spacing.x = unit(1.2, "lines"),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 24),
    
    axis.title.x = element_text(size = 26, margin = margin(t = 8)),
    axis.title.y = element_text(size = 24, margin = margin(r = 8)),
    
    legend.text = element_text(size = 24)
  )


p2 <- p1 +
  stat_compare_means(
    aes(group = Diet.x),
    method = "wilcox.test",
    label = "p.signif",
    size = 10,
    label.x.npc = 0.5,   # middle
    label.y.npc = 0.98   # near top
  )

#Save the plot
ggsave(file.path(output_dir,"Alpha_diversity_index.tiff"),
  plot = p2,
  width = 15,
  height = 8,
  units = "in",
  dpi = 900, compression = "lzw"
)


#Statistics
indices <- c(
  "Observed",
  "Shannon",
  "Simpson",
  "InvSimpson",
  "Pielou"
)

lapply(indices, function(x){
  
  formula <- as.formula(
    paste(x,"~ Diet")
  )
  
  wilcox.test(formula,
              data = alpha_div)
  
})

#Save the results
write.csv(alpha_div, file.path(output_dir, "Alpha_Diversity.csv"), row.names = FALSE)
