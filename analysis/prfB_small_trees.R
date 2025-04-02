library(treeio)
library(ggtree)
library(phangorn)
library(tidyverse)
library(micro.gen.extra)

setwd("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Bioinformatics")

df_full = read.csv("FS_data_clean_8_1_24.csv") 

df = df_full %>%
  select(assembly, nuccore, in_frame_stop., internal_stop, phylum, assemblyStats.gcPercent) %>%
  rename(stop_identity = internal_stop, stop_presence = in_frame_stop., gc = assemblyStats.gcPercent)
rownames(df) = df$assembly

# Upload list of accesions, 5 genomes per phylum (with >10 genomes), were randomly selected for the tree. 
df_tree = read.table("assemblies_for_prfb_aa_tree.txt", col.names = "assembly") %>%
  left_join(df, by = join_by(assembly))
  

tree = read.newick("genomes_phylo.tre")
tree$tip.label = str_extract(tree$tip.label, "[^_]*_[^_]*")
tree_mid = midpoint(tree)

prfB_tree = read.newick("prfB_AA_small.tre")
prfB_tree_mid = midpoint(prfB_tree)



d1 = t1$data[t1$data$isTip,]  
d1$x[] = 1  
d2 = t2$data[t2$data$isTip,]  
d2$x[] = 2  

TTcon = rbind(d1, d2)  
# From... https://www.biostars.org/p/132909/


t1 = ggtree(tree_mid, branch.length = "none") %<+% df_tree + 
  theme_tree(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +
  #geom_tiplab(aes(label=phylum), align = TRUE, size = 1) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage")


t2 = ggtree(prfB_tree_mid, branch.length = "none") %<+% df_tree + 
  theme_tree(legend.position='none', plot.margin = unit(c(0,0,0,0),"cm")) +
  #geom_tiplab(aes(label=phylum), align = TRUE, size = 2) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  scale_x_reverse()

cols = sample(hue_pal()(19))

l1 = ggplot(TTcon, aes(x = x, y = y, color = phylum, group = label)) + 
  geom_line() +   
  theme_void() + 
  theme(legend.position="none", plot.margin = unit(c(1,0,1,0),"cm")) +
  scale_colour_manual(values=cols)

p = plot_grid(t1, l1 ,t2, nrow = 1, align = "hv")

ggsave("C:/Users/cassp/Box Sync/Feaga Lab/Cassidy Prince/prfB/Figures/coevolution_tree.png", p, width = 10, height = 6, dpi = 600, units = "in")


l2 = ggplot(TTcon, aes(x = x, y = y, color = phylum, group = label)) + 
  geom_line() +   
  theme_void() + 
  theme(plot.margin = unit(c(1,0,1,0),"cm")) +
  scale_colour_manual(values=cols)
l2






ggtree(tree_mid, branch.length = "none") %<+% df_tree + 
  theme_tree(plot.margin = unit(c(0,0,0,0),"cm")) +
  geom_tiplab(aes(label=phylum), align = TRUE, size = 3) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage")



###########


ggtree(prfB_tree_mid) %<+% df_tree + 
  xlim(NA, 2) + 
  geom_tiplab(align = TRUE, size = 3) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage")


ggsave("prfB_AA_small.png", test, width = 10, height = 6.4, dpi = 600, units = "in")


ggtree(prfB_tree_mid) %<+% df_tree + 
  xlim(NA, 8) + 
  scale_x_reverse() +
  geom_tiplab(align = TRUE, size = 2) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 1.5, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage")


# Tree with internal node numbers 
ggtree(prfB_tree_mid) %<+% df_tree + 
  xlim(NA, 8) + 
  scale_x_reverse() +
  geom_tiplab(aes(label=phylum), align = TRUE, size = 2) + 
  geom_text(aes(label=node), size = 3) 

clade_nodes = c(166, 171, 174, 179, 184, 185, 115, 110, 107, 102, 98, 160, 123, 132, 128, 154, 151, 156, 145, 141)




p = ggtree(prfB_tree_mid) %<+% df_tree + 
  xlim(NA, 8) + 
  scale_x_reverse() +
  geom_tiplab(aes(label=phylum), align = TRUE, size = 2) 

collapse_tree(p, clade_nodes)


#Maybe tree_subset instead...