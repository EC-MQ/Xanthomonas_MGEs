#library(phyloseq)
library(dplyr)
library(ggtree)
library(ape)
library(gridExtra)
library(ggtreeExtra)
library(phangorn)
library(ggplot2)
library(ggnewscale)

setwd("/path_to/pretty_tree")

tree <- read.tree("ICE_core_with_other_genera.fasttree")
plot(tree)
rtree <- midpoint(tree)
plot(rtree)

rtree$tip.label

info <- read.csv("ICE_other_gen_metadata.csv", sep=',')
data <- info %>% select(ID_names,species)
rownames(data) <- data$ID_names
data <- data %>% select(-ID_names)



rn <- rownames(data)
heatmapData <- as.data.frame(sapply(data, as.character))
rownames(heatmapData) <- rn
unique(data$species)

#COLOURS!!!

tip_colors <- setNames(data$species[match(rtree$tip.label, rownames(data))], rtree$tip.label)
color_mapping <- c(
  
  " Achromobacter sp." =  "turquoise",
  "aeruginosa" =  "pink",
  "campestris" = "#DFFF00",
  "caspiana" = "purple4",
  "cissicola" = "#3498DB",
  
  "Cupriavidus sp." =  "brown4",
  "syringae" =  "green3",
  "viridiflava" = "darkgreen"
)


p <- ggtree(rtree, layout = "rectangular", size = 0.05) +
  geom_tiplab(align=TRUE, linetype="dashed", size=2, linesize=0.1,) +
  geom_treescale(fontsize = 2, linesize = 0.05) +
  
  # TIP COLORS (categorical)
  geom_tippoint(aes(color = tip_colors[label]), size = 1.8) +
  scale_color_manual(
    values = color_mapping,
    na.value = "transparent",
    name = "Tip group"
  ) +
  
  # allow a new color scale
  new_scale_color() +
  
  # NODE SUPPORT (numeric)
  geom_point(
    data = p$data %>%
      filter(!isTip) %>%
      mutate(support = as.numeric(label)) %>%
      filter(!is.na(support) & support > 0.80),
    aes(x = x, y = y, color = support),
    size = 1,
    shape = 18
  ) +
  scale_color_gradient(
    low = "gray90",
    high = "black",
    name = "Support"
  )
p

