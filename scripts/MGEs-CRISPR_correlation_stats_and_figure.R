library(dplyr)


##you need to inputs for this script, "genomes_group_matrices" and "crispr-cas_summary" file derived by the output from Defense finder
## genomes_group_matrices look like this

'''
genome	MGE_group_1	MGE_group_3	MGE_group_7	MGE_group_8
genomeA	0	0	0	0
genomeB	0	0	0	0
genomeC	1	0	0	0
genomeD	0	0	1	0
'''

##where per each genome there the number of of that MGE type
##I had one table per MGE, one for ICE, one for IME and one for plasmid, so you will see here I am merging them
##the groups where defined as in the paper (based on ANI)


# the crispr-cas_summary file should look like this which is basically a list of the genomes and a column with yes/no based on presence of CRISPR-Cas

'''
genome	crispr-cas	subtype
genomeA	yes	CAS_Class1-Subtype-I-F
genomeB	yes	CAS_Class1-Subtype-I-F
genomeC	no	no
'''

data_ICE <- read.csv("/path_to/group_matrixes/ICE_genome_group_matrix.csv", header = TRUE)
data_IME <- read.csv("/path_to/group_matrixes/IME_genome_group_matrix.csv", header = TRUE)
data_pl <- read.csv("/path_to/group_matrixes/plasmids_genome_group_matrix.csv", header = TRUE)

# Merge the dataframes
data0 <- merge.data.frame(data_pl,data_IME, by = "genome", all = TRUE)
data <- merge.data.frame(data0,data_ICE, by = "genome", all = TRUE)

data[is.na(data)] <- 0
data$row_sum <- rowSums(data[, !names(data) %in% "genome"])


totalMGEs <- data.frame(
  genome = data$genome,
  sum = data$row_sum,
  stringsAsFactors = FALSE
)

data_crispr <- read.csv("~/path_to/genomes/crispr_cas/crispr-cas_summary.csv")

data_MGE_CRISPR <- merge.data.frame(totalMGEs, data_crispr, by.x = "genome", by.y = "genome", all = TRUE)

data_MGE_CRISPR <- data_MGE_CRISPR %>%
  select(-subtype)


# Make sure there are no NAs
data_MGE_CRISPR$sum[is.na(data_MGE_CRISPR$sum)] <- 0
data_MGE_CRISPR$crispr.cas[is.na(data_MGE_CRISPR$crispr.cas)] <- "no"

# Calculate correlation
mwu_result <- wilcox.test(sum ~ crispr.cas, data = data_MGE_CRISPR, exact = FALSE)
# Calculate group means (for text)
means <- aggregate(sum ~ crispr.cas, data = data_MGE_CRISPR, mean)

#Print Results
print("--- Mann-Whitney U Test Results ---")
print(mwu_result)
print("--- Mean MGEs per Group ---")
print(means)


# Poisson GLM 
mge_model <- glm(sum ~ crispr.cas, 
                 data = data_MGE_CRISPR, 
                 family = poisson)

summary(mge_model)

##or

mge_model_adj <- glm(sum ~ crispr.cas, 
                     family = quasipoisson, 
                     data = cdata_MGE_CRISPR)
summary(mge_model_adj)

# Make a pretty plot
library(ggplot2)

pretty_plot<- ggplot(data_MGE_CRISPR, aes(x = crispr.cas, y = sum, fill = crispr.cas)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_y_continuous(
    limits = c(0, 9),          # Sets the total height of the axis
    breaks = c(0, 3, 6, 9)     # Sets the specific lines and numbers
  ) +
  labs(title = "MGE Abundance vs CRISPR Presence",
       x = "CRISPR-Cas System",
       y = "Number of MGEs per Genome") +
  theme_minimal() +
  scale_fill_manual(values = c("no" = "grey80", "yes" = "purple"))


pretty_plot

