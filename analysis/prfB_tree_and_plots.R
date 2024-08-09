# Load packages and set working directory.

library(treeio)
library(ggtree)
library(ggprism)
library(ggbreak)
library(ggnewscale)
library(ggpubr)
library(tidyverse)
library(castor)
library(phangorn)
library(aplot)
library(grDevices)
library(scales)
library(ggpmisc)

setwd("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Bioinformatics")

### --- DATA PREP --- ###

# Import clean data and select important columns.
df_full = read.csv("FS_data_clean_8_1_24.csv") 

df = df_full %>%
  select(assembly, nuccore, in_frame_stop., internal_stop, phylum, assemblyStats.gcPercent) %>%
  rename(stop_identity = internal_stop, stop_presence = in_frame_stop., gc = assemblyStats.gcPercent)
rownames(df) = df$assembly

# Four Actinobacteria in the dataset claim to have an internal stop codon, but based on manual inspection their prfB genes are misannotated. Fixing the data to reflect this.
actinos = df %>%
  filter(phylum == "Actinobacteriota") %>%
  mutate(stop_presence = "no", stop_identity = "no stop")

df = df %>%
  filter(phylum != "Actinobacteriota") %>%
  bind_rows(actinos) %>%
  mutate()

df$stop_identity[df$stop_identity != "TAG" & df$stop_identity != "TAA" & df$stop_identity != "TGA"] = "no stop"

# Make a separate dataframe that has only phyla with more than 10 representatives.
df_filt = df %>%
  group_by(phylum) %>%
  filter(n() > 10)

df1 = data.frame(df$stop_presence)
rownames(df1) = rownames(df)
df2 = data.frame(df$stop_identity)
rownames(df2) = rownames(df)
df3 = data.frame(df$phylum)
rownames(df3) = rownames(df)


### --- TREE VISUALIZATION --- ###

# Import 16S rRNA maximum likelihood tree.
tree = read.newick("C:\\Users\\cassp\\OneDrive\\Documents\\Feaga Lab\\prfB\\16S_fasttree.tre")

# Upload NCBI RefSeq Assembly database accessions (beginning with "GCF_") and corresponding NCBI Nucleotide (nuccore) database accessions (beginning with "NC_" or "NZ_". 
# Rename tree tips to be the Assembly accessions instead on Nucleotide accessions.
df_GCF_nc = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//GCF_nuccore_reps_clean.csv")
tree = get_subtree_with_tips(tree, only_tips = df_GCF_nc$nuccore)$subtree
tree_tip_GCF = left_join(data.frame(tree$tip.label), df_GCF_nc, by = join_by("tree.tip.label" == "nuccore"), multiple = "any")
tree$tip.label = tree_tip_GCF$assembly

# Get subtree for only genomes with prfB data.
subtree = get_subtree_with_tips(tree, only_tips = rownames(df1))$subtree
tree_mid = midpoint(subtree)

p = ggtree(tree_mid, layout='circular', size=0.2) 

p1 = gheatmap(p, df1, offset=-0.2, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("yes" = "#961415", "no" = "gray80"), labels = c("no frameshift", "frameshift"), na.value = "white") + 
  theme(text=element_text(size=18)) + 
  new_scale_fill()

p2 = gheatmap(p, df2, offset=-0.2, width=0.1, font.size=1, colnames = FALSE, color=NA) + 
  scale_fill_manual(values=c("no stop" = "gray80", "TGA" = "#961415", "TAA" = "#520e15", "TAG" = "white"), name="Internal stop codon \nidentity or modification", na.value = "white") +
  theme(text=element_text(size=15)) +
  new_scale_fill()

p3 = gheatmap(p1, df3, offset=0.3, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=15)) +
  scale_fill_discrete(name = "Phylum/group", na.value = "white")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\16S_tree_phy_8_1_24.png", p3, units = "in", width = 17, height = 13, dpi = 600)

### --- TREE WITH BARCHART --- ###

# Select only phyla with more than 10 genomes.

n_vals = data.frame(table(df_filt$phylum))

# Identify teh proportion of genomes that have the frameshift in each phylum.
table_all = as.data.frame(prop.table(table(df_filt$phylum, df_filt$stop_presence), margin = 1)*100)
table_yes = table_all[table_all$Var2 == "yes",]

# Randomly sample one assembly per phylum for the collapsed tree. This was performed once and the assembly accessions were saved to the file "random_genomes_new.csv" for reuse. 
random = df_filt %>%
  rownames_to_column('rowname') %>%
  group_by(phylum) %>%
  sample_n(1) %>%
  column_to_rownames('rowname')

random = random %>%
  inner_join(table_yes, by = c('phylum' = 'Var1'))

#write.csv(random, "random_genomes_new.csv")

random = read.csv("random_genomes_new.csv")
random = column_to_rownames(random, 'X')

# Filter tree for the random representative genomes and midpoint root.
subtree_bar = get_subtree_with_tips(tree, only_tips = rownames(random))$subtree
tree_mid = midpoint(subtree_bar)

p = ggtree(tree_mid, size = 0.8) %<+% random + xlim(NA, 9) + geom_tiplab(aes(label=phylum), align = TRUE, size = 6)

p1 = ggplot(random, aes(assembly, Freq)) + 
  geom_col(color="black", fill = "gray80", width = 1, linewidth = 0.7) +
  theme_prism() +
  xlab("") +
  ylab("% genomes with frameshift") +
  theme(axis.text.y=element_blank(), text=element_text(size=17)) +
  scale_y_continuous(expand= c(0,0), limits = c(0,110)) +
  geom_text(position = position_dodge2(preserve = 'single',width = 0.9), hjust = -.3, aes(label = n_vals$Freq, size = 17)) +
  coord_flip()

p2 = p1 %>% insert_left(p, width = 2)


png(filename = "C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\collapsed_tree_8_6_24.png", units = "in", width = 16, height = 9, res = 600)

p2

dev.off()

### --- OTHER PLOTS --- ###

# GC content violin plots
df_gc = df %>%
  mutate(FS = recode(stop_identity, TGA = "frameshift", TAA = "frameshift", 'no stop' = "no frameshift"))

summ = df_gc %>%
  group_by(FS) %>%
  summarize(n = paste0("n = ", n()), gc = 81)

vplot = ggplot(df_gc, aes(x = factor(FS, level=c('no frameshift', 'frameshift')), y = gc, fill = FS))+
  geom_violin(trim=TRUE, width = 0.7)+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  theme_prism()+
  scale_y_continuous(expand= c(0,0), limits = c(20, 84))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#961415"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("GC content (%)") + 
  theme(legend.position = "none") +
  geom_text(data = summ, aes(label = n)) +
  stat_compare_means(method = "t.test", label.y = 83, label.x = 1.5, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\GC_violin_8_1_24.png", vplot, width = 4.5, height = 6, dpi = 600, units = "in") 

# Supplemental Figure *** gc no Actinobacteria.
df_gc_no_act = df_gc %>%
  filter(phylum != "Actinobacteriota")

summ_no_act = df_gc_no_act %>%
  group_by(FS) %>%
  summarize(n = paste0("n = ", n()), gc = 80)

vplot_no_act = ggplot(df_gc_no_act, aes(x = factor(FS, level=c('no frameshift', 'frameshift')), y = gc, fill = FS))+
  geom_violin(trim=TRUE, width = 0.7)+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  theme_prism()+
  scale_y_continuous(expand= c(0,0), limits = c(20, 89))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#961415"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("GC content (%)") + 
  theme(legend.position = "none") +
  geom_text(data = summ_no_act, aes(label = n)) +
  stat_compare_means(method = "t.test", label.y = 85, label.x = 1.5, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\GC_violin_no_act_8_6_24.png", vplot_no_act, width = 4.5, height = 6, dpi = 600, units = "in")

# What is the mean GC content for FS vs no FS in both GC datasets?
df_gc %>%
  group_by(FS) %>%
  summarize(mean = mean(as.numeric(gc)))

df_gc_no_act %>%
  group_by(FS) %>%
  summarize(mean = mean(as.numeric(gc)))

# Figure **** Premature stop codon usage barchart
df_stop = df %>% 
  filter(stop_identity != "no stop") %>%
  group_by(stop_identity) %>%
  summarize(n = n()) %>%
  add_row(stop_identity = "TAG", n = 0)

legend_ord = levels(with(df_stop, reorder(stop_identity, -n)))
yticks = c(0, 100, 200, 2000, 4000, 6000, 8000)

stop_plot = ggplot(df_stop, aes(x = factor(stop_identity, level=c('no stop', 'TGA', 'TAA', 'TAG')), y = n, fill = stop_identity)) +
  geom_col(color = "black") + 
  theme_prism() +
  scale_fill_manual(breaks = legend_ord, values = c("TGA" = "#961415", "TAA" = "#520e15", "TAG" = "#EFAFB3")) +
  scale_y_continuous(expand= c(0,0), limits = c(0, 8500), breaks = yticks) +
  scale_y_break(c(210, 2000), scales = 4, expand= c(0,0), space = 0.5) +
  geom_text(aes(y=n+1, label=n), vjust= -0.5, color="black", size=5) +
  geom_text(aes(y=n+1, label=paste0(round(100*n/sum(n), 1), "%")), vjust= 1.6, color="white", size=5) +
  xlab("Premature stop codon identity") +
  ylab("Number of genomes") +
  theme(text = element_text(size = 18))

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\stop_identity_8_1_24.png", stop_plot, width =  6, height = 5, units = "in")

# TGA codon usage violin plot

df_cds = data.frame(read.csv("cds_new_stops.csv"))

df_term_stops = data.frame(table(df_cds$genome_ID[df_cds$terminal_stop == "TAA" | df_cds$terminal_stop == "TAG" | df_cds$terminal_stop == "TGA"], df_cds$terminal_stop[df_cds$terminal_stop == "TAA" | df_cds$terminal_stop == "TAG" | df_cds$terminal_stop == "TGA"]))%>% 
  inner_join(df_GCF_nc, join_by(Var1 == nuccore)) %>%
  inner_join(df, by = 'assembly') %>%
  mutate(stop_presence = recode(stop_presence, no = "no frameshift", yes = "frameshift")) %>%
  group_by(assembly, Var2) %>% 
  mutate(Sum=sum(Freq)) %>% 
  ungroup() %>%
  group_by(assembly) %>%
  mutate(total_stops = sum(Freq)) %>%
  ungroup() %>%
  mutate(Prop = 100*Sum/total_stops) %>%
  distinct(assembly, Var2, .keep_all = TRUE)
  
props_TGA = df_term_stops %>%
  filter(Var2 == "TGA")

summ = props_TGA %>%
  group_by(stop_presence) %>%
  summarize(n = paste0("n = ", n()), Prop = 110)

plot = ggplot(data = props_TGA, aes(x = factor(stop_presence, level=c('no frameshift', 'frameshift')), y = Prop, fill = stop_presence))+
  geom_violin(trim=TRUE, width = 0.7)+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  theme_prism()+
  scale_y_continuous(expand= c(0,0), limits = c(0, 115), breaks = c(0, 25, 50, 75, 100))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#961415"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("TGA codon usage (%)") + 
  theme(legend.position = "none") +
  geom_text(data = summ, aes(label = n)) +
  stat_compare_means(method = "t.test", label.y = 105, label.x = 1.4, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\TGA_usage_violin_8_5_24.png", plot, width = 4.5, height = 6, dpi = 600, units = "in")

# Supplemental Figure **** TGA usage Remove Actinobacteriota and rerun.
df_term_stops_no_act = df_term_stops %>%
  filter(phylum != "Actinobacteriota")

props_TGA_no_act = df_term_stops_no_act %>%
  filter(Var2 == "TGA")

summ = props_TGA_no_act %>%
  group_by(stop_presence) %>%
  summarize(n = paste0("n = ", n()), Prop = 100)

plot = ggplot(data = props_TGA_no_act, aes(x = factor(stop_presence, level=c('no frameshift', 'frameshift')), y = Prop, fill = stop_presence))+
  geom_violin(trim=TRUE, width = 0.7)+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  theme_prism()+
  scale_y_continuous(expand= c(0,0), limits = c(0, 110), breaks = c(0, 25, 50, 75, 100))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#961415"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("TGA codon usage (%)") + 
  theme(legend.position = "none") +
  geom_text(data = summ, aes(label = n)) +
  stat_compare_means(method = "t.test", label.y = 105, label.x = 1.4, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\TGA_usage_violin_noact_8_5_24.png", plot, width = 4.5, height = 6, dpi = 600, units = "in")


# Supplemental Figure ### GC content correlation with stop codon usage

df_term_stops$Var2 = factor(df_term_stops$Var2, levels=c("TAA", "TGA", "TAG"))
plot = ggplot(df_term_stops, aes(x = gc, y = Prop, color = Var2)) +
  geom_point(size = 3, alpha = 0.7)+
  theme_prism()+
  scale_color_manual(values=c("TAA" = "#520E15", "TGA" = "#961415", "TAG" = "#CB757C"))+
  scale_y_continuous(expand= c(0,0), limits = c(0, 100)) +
  theme(text = element_text(size = 17)) +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2")), label.x = "center") +
  xlab("GC content (%)")+
  ylab("% of stop codons")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\GC_stop_usage_8_6_24.png", plot, width = 7.5, height = 5.5, dpi = 600, units = "in")

