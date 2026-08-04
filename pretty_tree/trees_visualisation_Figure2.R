##you probably don't need all these packages##

library(ape)
library(phangorn)
library(phyloseq)
library(dplyr)
library(ggtree)
library(gridExtra)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(tidyr)
library(tibble)


setwd("/path_to/pretty_tree")

#prepare the metadata
metadata <- read.csv("species_num_MGEs.csv",row.names = 1, header = TRUE)

#clean the rownames so match the ones in the tree
row.names(metadata)<- gsub("\\.", "_", rownames(metadata))

clean_labels <- function(labels) {
  # This pattern matches: (Segment1_Segment2_Segment3) _ (Anything else)
  # It keeps only the part in the parentheses.
  sub("^([^_]+_[^_]+_[^_]+)_.*", "\\1", labels)
}

rownames(metadata) <- clean_labels(rownames(metadata))
metadata_species <- metadata %>% select(Species)
metadata_MGEs <- metadata %>% select(Number.of.IMEs, Number.of.ICEs, Number.of.plasmids)

#crispr-cas metadata
metadata_crispr <- read.csv("3species_crispr_cas_summary.csv",row.names = 1, header = TRUE)
rownames(metadata_crispr) <- clean_labels(rownames(metadata_crispr))


#spacers metadata
metadata_spacers <- read.csv("spacernumber.csv",row.names = 1, header = TRUE)
row.names(metadata_spacers)<- gsub("\\.", "_", rownames(metadata_spacers))
rownames(metadata_spacers) <- clean_labels(rownames(metadata_spacers))

metadata_spacers_filtered <- metadata_spacers[
  rownames(metadata_spacers) %in% rownames(metadata_crispr),
  ,
  drop = FALSE
]
#import tree here
tree <- read.tree("3_species_tree.tree")

#check the tip labels and modify
tree$tip.label
tree$tip.label <- gsub("_genomic", "", tree$tip.label)
tree$tip.label <- gsub("\\.", "_", tree$tip.label)
tree$tip.label <-clean_labels(tree$tip.label)

#plot(tree)
rtree <- root(tree, "GCA_030168935_1")
rtree <- drop.tip(rtree, "GCA_030168935_1")
plot(rtree)
tip_colors <- setNames(metadata_species$Species[match(rtree$tip.label, rownames(metadata_species))], rtree$tip.label)


###COLORS!!!

color_mapping <-c("Xanthomonas campestris" = "#F9E79F", # Soft Maize Yellow
                  "Xanthomonas cissicola"  = "#AED6F1", # Soft Sky Blue [cite: 13]
                  "Xanthomonas oryzae"     = "#F5B7D1") # Pale Rose [cite: 14]

heatmap.colours_crispr <- c(
  #yes
  "#9B59B6",
  #1
  "black")

mge_colors <- c(
  "Number.of.IMEs"      = "#2C3E50", # Muted Slate (Replaces harsh Black)
  "Number.of.ICEs"      = "#E67E22", # Soft Peach (Replaces bright Orange)
  "Number.of.plasmids"  = "#16A085"  # Soft Mint/Teal (Replaces bright Red)
)

rn <- rownames(metadata_crispr)

heatmapData_crispr <- as.data.frame(sapply(metadata_crispr, as.character))
rownames(heatmapData_crispr) <- rn



#Create the base plot with bootstrap
p <- ggtree(rtree, layout = 'rectangular', size = 0.05) +
  geom_treescale(fontsize = 2) +
  geom_point(data = td_filter(!isTip & !is.na(as.numeric(label)) & as.numeric(label) > 0.80), 
             aes(color = as.numeric(label)), size = 0.1) +  
  scale_color_gradient(low = "gray90", high = "black", name = "Bootstrap")


p

##this one highlight the tip labels based on species

# We convert row names to a column called 'label' so it matches the tree data
metadata_to_join <- metadata_species %>% 
  mutate(label = rownames(.))

p$data <- p$data %>% 
  left_join(metadata_to_join, by = "label")

# Check if 'Species' actually exists now
if(!"Species" %in% colnames(p$data)) {
  stop("Species column still missing! Check your row name matching.")
}

# Add the species layer using the new scale
p_final <- p + 
  new_scale_color() + 
  geom_tippoint(aes(color = Species), size = 1) + 
  scale_color_manual(values = color_mapping)

p_final



##add crispr-cas column
p1 <- gheatmap(p_final, heatmapData_crispr,
               offset = 0,
               #name of the coloumn?
               colnames = TRUE, 
               #colnames with and angle?
               colnames_angle=30,
               #colname position
               colnames_position='top',
               width=0.05,
               color = NA
) +
  scale_fill_manual(values=heatmap.colours_crispr, na.value = "white")

p1


##need to add nor the number of MGEs
mge_long <- metadata_MGEs %>%
  rownames_to_column("label") %>% # Turn rownames into a 'label' column to match the tree
  pivot_longer(
    cols = c(Number.of.IMEs, Number.of.ICEs, Number.of.plasmids), 
    names_to = "MGE_Type", 
    values_to = "Count"
  )


# create a dummy row with a Count of 8 --> you will have to remove this manually from the figure with illustrator
# We attach it to the first tip label so it has a place to live on the Y-axis
ghost_row <- data.frame(
  label = "GCA_001928995_2", 
  MGE_Type = mge_long$MGE_Type[1], 
  Count = 10
)

# combine it with  real data
mge_long_with_ghost <- rbind(mge_long, ghost_row)

# pretty plot
p2 <- p1 + 
  new_scale_fill() + 
  geom_fruit(
    data = mge_long_with_ghost,
    geom = geom_bar,
    mapping = aes(y = label, x = Count, fill = MGE_Type),
    pwidth = 0.5, 
    orientation = "y", 
    stat = "identity",
    offset = 0.1,
    alpha = 0.9, # Slight transparency so ghost bar isn't jarring if it overlaps
    axis.params = list(
      axis = "x",
      at = c(0, 2, 4, 6, 8, 10),
      labels = c("0", "2", "4", "6", "8", "10"),
      text.size = 2
    ),
    grid.params = list(
      linetype = 2,
      color = "gray"
    )
  ) +
  # Use scale_fill_manual for your custom colors
  scale_fill_manual(
    values = mge_colors, 
    name = "MGE Categories"
    #    labels = c("IME", "ICE", "Plasmid") # This cleans up the legend text
  )
p2



