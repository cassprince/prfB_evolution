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

p2 = gheatmap(p1, df3, offset=0.3, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=15)) +
  scale_fill_discrete(name = "Phylum/group", na.value = "white")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\16S_tree_phy_10_9_24.png", p2, units = "in", width = 17, height = 13, dpi = 600)

# Figure S1: Stop identity on tree

sup_df1 = data.frame(df$stop_identity)
rownames(sup_df1) = rownames(df)

sup_p1 = gheatmap(p, sup_df1, offset=-0.2, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("no stop" = "gray80", "TGA" = "#961415", "TAA" = "#520e15", "TAG" = "white"), name="Internal stop codon \nidentity or modification", na.value = "white") + 
  theme(text=element_text(size=18)) + 
  new_scale_fill()

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\FigS1_8_26_24.png", sup_p1, units = "in", width = 17, height = 13, dpi = 600)

### --- FIGURE 3A AND FIGURE S2: TREES WITH HEATMAPS --- ###

# Select only phyla with more than 10 genomes.

n_vals = data.frame(table(df_filt$phylum))

# Identify the proportion of genomes that have the frameshift in each phylum.
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

random = read.csv("random_genomes_new.csv") %>%
  select(-X)
rownames(random) = random$assembly
random$n = n_vals$Freq

# Filter tree for the random representative genomes and midpoint root.
subtree_bar = get_subtree_with_tips(tree, only_tips = random$assembly)$subtree
tree_mid = midpoint(subtree_bar)

# Make tree with all phyla (Figure S2).
p_all = ggtree(tree_mid, size = 0.8) %<+% random + 
  xlim(NA, 15) + 
  geom_tiplab(aes(label=phylum), align = TRUE, size = 6) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage", limits = c(0,100)) +
  new_scale_fill()


all = gheatmap(p_all, random %>% select(Freq), offset = 10, width=0.4, font.size=3.5, colnames = FALSE, color=NA, colnames_angle = 90, colnames_offset_y = -0.2) + 
  scale_fill_viridis_c(option="F", direction = -1, name="Percent\nwith frameshift") 

all

text_all = ggplot(random, aes(assembly, y = 0)) + 
  geom_text(hjust = 0, aes(label = n, size = 17)) +
  coord_flip() +
  theme_void() +
  theme(legend.position="none")
  
all_n = text %>% insert_left(all, width = 10)


ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\full_phyla_heatmap_3_7_25.png", all_n, width = 11, height = 6.4, dpi = 600, units = "in")

# Deep root only (Figure 3A).
deep_phyla = c("Spirochaetota", "Deinococcota", "Fusobacteriota", "Synergistota", "Thermotogota")

random_deep = random %>%
  filter(phylum %in% deep_phyla)

subtree_deep = get_subtree_with_tips(tree, only_tips = random_deep$assembly)$subtree
tree_mid_deep = midpoint(subtree_deep)


p_deep = ggtree(tree_mid_deep, size = 0.8) %<+% random_deep + 
  xlim(NA, 15) + 
  geom_tiplab(aes(label=phylum), align = FALSE, size = 6) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage", limits = c(0,100)) +
  new_scale_fill()


deep = gheatmap(p_deep, random_deep %>% select(Freq), offset = 7, width=0.4, font.size=3.5, colnames = FALSE, color=NA, colnames_angle = 90, colnames_offset_y = -0.2) + 
  scale_fill_viridis_c(option="F", direction = -1, name="Percent\nwith frameshift") +
  theme(text = element_text(size = 12))

deep

text_deep = ggplot(random_deep, aes(assembly, y = 0)) + 
  geom_text(hjust = 0, aes(label = n, size = 17)) +
  coord_flip() +
  theme_void() +
  theme(legend.position="none")

deep_n = text_deep %>% insert_left(deep, width = 7)

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\deep_phyla_heatmap_3_7_25.png", deep_n, width = 7, height = 3, dpi = 600, units = "in")

### --- OTHER PLOTS --- ###

# Figure 5A: GC % violin plot
df_gc = df %>%
  mutate(FS = recode(stop_identity, TGA = "frameshift", TAA = "frameshift", 'no stop' = "no frameshift"))

summ = df_gc %>%
  group_by(FS) %>%
  summarize(n = paste0("n = ", n()), gc = 81)

vplot = ggplot(df_gc, aes(x = factor(FS, level=c('no frameshift', 'frameshift')), y = gc, fill = FS))+
  geom_violin(trim=TRUE, width = 0.7)+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  theme_classic()+
  scale_y_continuous(expand= c(0,0), limits = c(20, 89))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#BA272F"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("GC content (%)") + 
  theme(legend.position = "none",
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  geom_text(data = summ, aes(label = n), , size = 5) #+
  #stat_compare_means(method = "t.test", label.y = 83, label.x = 1.5, label = "p.format")
vplot

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\GC_violin_3_17_25.png", vplot, width = 4, height = 4, dpi = 600, units = "in") 

# Figure 5C: GC % no Actinobacteriota.
df_gc_no_act = df_gc %>%
  filter(phylum != "Actinobacteriota")

summ_no_act = df_gc_no_act %>%
  group_by(FS) %>%
  summarize(n = paste0("n = ", n()), gc = 80)

vplot_no_act = ggplot(df_gc_no_act, aes(x = factor(FS, level=c('no frameshift', 'frameshift')), y = gc, fill = FS))+
  geom_violin(trim=TRUE, width = 0.7)+
  geom_boxplot(width=0.1, outlier.shape = NA)+
  theme_classic()+
  scale_y_continuous(expand= c(0,0), limits = c(20, 89))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#BA272F"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("GC content (%)") + 
  theme(legend.position = "none") +
  geom_text(data = summ_no_act, aes(label = n), size = 5) #+
  #stat_compare_means(method = "t.test", label.y = 85, label.x = 1.5, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\GC_violin_no_act_3_17_25.png", vplot_no_act, width = 4, height = 4, dpi = 600, units = "in")

# What is the mean GC content for FS vs no FS in both GC datasets?
df_gc %>%
  group_by(FS) %>%
  summarize(mean = mean(as.numeric(gc)))

df_gc_no_act %>%
  group_by(FS) %>%
  summarize(mean = mean(as.numeric(gc)))

# Figure 2B: Premature stop codon usage barchart
df_stop = df %>% 
  filter(stop_identity != "no stop") %>%
  group_by(stop_identity) %>%
  summarize(n = n()) %>%
  add_row(stop_identity = "TAG", n = 0)

legend_ord = levels(with(df_stop, reorder(stop_identity, -n)))
yticks = c(0, 100, 200, 2000, 4000, 6000, 8000)

stop_plot = ggplot(df_stop, aes(x = factor(stop_identity, level=c('no stop', 'TGA', 'TAA', 'TAG')), y = n, fill = stop_identity)) +
  geom_col(color = "black") + 
  theme_classic() +
  scale_fill_manual(breaks = legend_ord, values = c("TGA" = "#961415", "TAA" = "#520e15", "TAG" = "#EFAFB3")) +
  scale_y_continuous(expand= c(0,0), limits = c(0, 8500), breaks = yticks) +
  scale_y_break(c(210, 2000), scales = 3, expand= c(0,0), space = 0.2) +
  geom_text(aes(y=n+1, label=n), vjust= -0.5, color="black", size=3.8) +
  geom_text(aes(y=n+1, label=paste0(round(100*n/sum(n), 1), "%")), vjust= 1.6, color="white", size=3.8) +
  xlab("Premature stop \ncodon identity") +
  ylab("Number of genomes") +
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 15), 
        legend.position="none",
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))
stop_plot


ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\stop_identity_3_13_25.png", stop_plot, width =  3, height = 3.5, units = "in", dpi = 600)

# Figure 5B: TGA codon usage violin plot

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
  theme_classic()+
  scale_y_continuous(expand= c(0,0), limits = c(0, 115), breaks = c(0, 25, 50, 75, 100))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#BA272F"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("TGA codon usage (%)") + 
  theme(legend.position = "none",
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  geom_text(data = summ, aes(label = n), size = 5)# +
  #stat_compare_means(method = "t.test", label.y = 105, label.x = 1.4, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\TGA_usage_violin_3_17_25.png", plot, width = 4, height = 4, dpi = 600, units = "in")

# Figure 5D: TGA usage no Actinobacteriota.
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
  theme_classic()+
  scale_y_continuous(expand= c(0,0), limits = c(0, 110), breaks = c(0, 25, 50, 75, 100))+
  scale_fill_manual(values = c("no frameshift" = "gray80", "frameshift" = "#BA272F"), breaks=c("no frameshift", "frameshift")) +
  theme(text = element_text(size = 17))+
  xlab("")+
  ylab("TGA codon usage (%)") + 
  theme(legend.position = "none",
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  geom_text(data = summ, aes(label = n), size = 5) #+
  #stat_compare_means(method = "t.test", label.y = 105, label.x = 1.4, label = "p.format")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\TGA_usage_violin_noact_3_17_25.png", plot, width = 4, height = 4, dpi = 600, units = "in")


# Figure S3 GC content correlation with stop codon usage

df_term_stops$Var2 = factor(df_term_stops$Var2, levels=c("TAA", "TGA", "TAG"))
plot = ggplot(df_term_stops, aes(x = gc, y = Prop, color = Var2)) +
  geom_point(size = 3, alpha = 0.7)+
  theme_classic()+
  scale_color_manual(values=c("TAA" = "#520E15", "TGA" = "#961415", "TAG" = "#CB757C"), name = "Stop codon \nidentity")+
  scale_y_continuous(expand= c(0,0), limits = c(0, 100)) +
  theme(text = element_text(size = 20), 
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black")) +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2")), label.x = "center") +
  xlab("GC content (%)")+
  ylab("Stop codon usage (%)")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\GC_stop_usage_4_3_25.png", plot, width = 8, height = 5.7, dpi = 600, units = "in")


### --- TABLES --- ###

# Table S1
table_S1 = df_full %>%
  unite("taxonomy", domain:species, sep = ";") %>%
  rename(stop_identity = internal_stop, FS_presence = in_frame_stop., gc = assemblyStats.gcPercent) %>%
  mutate(FS_presence = recode(FS_presence, no = "no frameshift", yes = "frameshift")) %>%
  select(assembly, FS_presence, stop_identity, gc, taxonomy)

write.csv(table_S1, "C://Users//cassp//Cornell University//Heather Feaga - Cassidy prfB manuscript//Table_S1.csv",row.names = FALSE)

# Table S2

table_S2 = props_TGA %>%
  rename(proportion_TGA_stops = Prop, FS_presence = stop_presence, total_TGA_stops = Sum) %>%
  mutate(FS_presence = recode(FS_presence, no = "no frameshift", yes = "frameshift")) %>%
  select(assembly, FS_presence, total_TGA_stops, total_stops, proportion_TGA_stops)

write.csv(table_S2, "C://Users//cassp//Cornell University//Heather Feaga - Cassidy prfB manuscript//Table_S2.csv", row.names = FALSE)

# Summary information about Alphaproteobacteria

df_alphas = df_full %>%
  filter(phylum == "Alphaproteobacteria") 

df_alphas %>%
  group_by(in_frame_stop.) %>%
  summarize(mean_size = mean(assemblyStats.totalSequenceLength), mean_gc = mean(assemblyStats.gcPercent))

table(df_alphas$order, df_alphas$in_frame_stop.)

########## prfB DNA tree
tree_prfB = read.newick("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Bioinformatics\\prfB_3.tre")
tree_prfB$tip.label = gsub("'", "", tree_prfB$tip.label)
tree_prfB$tip.label = gsub("\\:.*", "", tree_prfB$tip.label)


tree_prfB = get_subtree_with_tips(tree_prfB, only_tips = df_GCF_nc$nuccore)$subtree
tree_prfB_tip_GCF = left_join(data.frame(tree_prfB$tip.label), df_GCF_nc, by = join_by("tree_prfB.tip.label" == "nuccore"), multiple = "any")

tree_prfB$tip.label = tree_prfB_tip_GCF$assembly

# Get subtree for only genomes with prfB data.
subtree_prfB = get_subtree_with_tips(tree_prfB, only_tips = rownames(df1))$subtree

tree_prfB_mid_2 = midpoint(subtree_prfB)

p = ggtree(tree_prfB_mid_2, layout='circular', size=0.2) 

p1 = gheatmap(p, df1, offset=-0.2, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("yes" = "#961415", "no" = "gray80"), labels = c("no frameshift", "frameshift"), na.value = "white") + 
  theme(text=element_text(size=18)) + 
  new_scale_fill()

p2 = gheatmap(p1, df3, offset=0.4, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=15)) +
  scale_fill_discrete(name = "Phylum/group", na.value = "white")

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\prfB_tree_phy_3_10_11_24.png", p2, units = "in", width = 17, height = 13, dpi = 600)

# Mycoplasmatota stops

df_GCF_nc = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//GCF_nuccore_reps_clean.csv")
df_GCF_nc$nuccore = gsub("NZ_JAZHCZ02", "NZ_JAZHCY01", df_GCF_nc$nuccore)

df_cds_myco = data.frame(read.csv("cds_myco_stops.csv"))
df_cds_myco$genome_ID = gsub("NZ_JAZHCZ02", "NZ_JAZHCY01", df_cds_myco$genome_ID)

df_term_stops_myco = data.frame(table(df_cds_myco$genome_ID[df_cds_myco$terminal_stop == "TAA" | df_cds_myco$terminal_stop == "TAG" | df_cds_myco$terminal_stop == "TGA"], df_cds_myco$terminal_stop[df_cds_myco$terminal_stop == "TAA" | df_cds_myco$terminal_stop == "TAG" | df_cds_myco$terminal_stop == "TGA"]))%>% 
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

props_TGA = df_term_stops_myco %>%
  filter(Var2 == "TGA")


### --- FIGURE 4: MYCOPLASMA TREE --- ###

df_myco = read.csv("myco_data_clean_10_10_24.csv")

df_myco_trna = read.table("myco_trna.csv", skip = 3, sep = "\t") 
colnames(df_myco_trna) = c("nuccore", "trna_number", "start", "end", "aa_type", "anticodon", "intron_start", "intron_end", "inf_score", "note")
df_myco_trna = df_myco_trna %>% 
  mutate(nuccore, nuccore = str_trim(nuccore))%>%
  left_join(df_GCF_nc, by = join_by("nuccore")) %>%
  select(-X)

df_myco_final = data.frame(table(df_myco_trna$assembly, df_myco_trna$aa_type)) %>%
  filter(Var2 == "Sup") %>%
  right_join(df_myco, by = join_by("Var1" == "assembly")) %>%
  select(-X, -Var2) %>%
  rename(assembly = Var1, sup_count = Freq) %>%
  mutate(prfB = assembly %in% df$assembly) %>%
  left_join(select(df, assembly, stop_presence)) %>%
  left_join(props_TGA %>% select(assembly, Prop))

df_prfb = data.frame(df_myco_final$prfB)
rownames(df_prfb) = df_myco_final$organism.organismName
df_FS = data.frame(df_myco_final$stop_presence)
rownames(df_FS) = df_myco_final$organism.organismName
df_sup = data.frame(df_myco_final$sup_count)
df_sup$df_myco_final.sup_count = as.character(df_sup$df_myco_final.sup_count)
rownames(df_sup) = df_myco_final$organism.organismName
df_gc = data.frame(df_myco_final$assemblyStats.gcPercent)
rownames(df_gc) = df_myco_final$organism.organismName
df_prop = data.frame(df_myco_final$Prop)
rownames(df_prop) = df_myco_final$organism.organismName
df_fam = data.frame(df_myco_final$family)
rownames(df_fam) = df_myco_final$organism.organismName


tree_myco = read.newick("myco_phylophlan.tre")
tree_myco$tip.label = paste0(str_split_i(tree_myco$tip.label, "_", 1), "_", str_split_i(tree_myco$tip.label, "_", 2))



tree_myco_tip_GCF = data.frame(tree_myco$tip.label) %>%
  left_join(df_myco_final %>% select(assembly, organism.organismName), by = join_by("tree_myco.tip.label" == "assembly"))

tree_myco$tip.label = tree_myco_tip_GCF$organism.organismName
tree_myco$node.label = as.numeric(tree_myco$node.label)*100

keep_tips = distinct(data.frame(tree_myco$tip.label))

tree_myco = get_subtree_with_tips(tree_myco, only_tips = keep_tips$tree_myco.tip.label)$subtree

tree_myco_mid = midpoint(tree_myco)

p = ggtree(tree_myco_mid, layout='rectangular', size=0.2) +
  geom_treescale(x = 0.06, y = 100, linesize=0.5, fontsize=3) +
  geom_tiplab(size = 2.5) +
  geom_nodepoint(aes(fill = as.numeric(label)), size = 1, shape = 21) +
  xlim(NA, 2.5) +
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  new_scale_fill()

p1 = gheatmap(p, df_prfb, offset=1.4, width=0.05, colnames = FALSE, color=NA) +
  scale_fill_manual(values = c("TRUE" = "#961415", "FALSE" = "gray80"), name = "prfB present?", na.value = "white") + 
  new_scale_fill()

p2 = gheatmap(p1, df_sup, offset=1.45, width=0.05, colnames = FALSE, color=NA) +
  scale_fill_manual(values = c("0" = "gray80", "1" = "#961415", "2" = "#961415"), name = "number of supressor tRNAs", na.value = "white") + 
  new_scale_fill()

p3 = gheatmap(p2, df_FS, offset=1.5, width=0.05, colnames = FALSE, color=NA) +
  scale_fill_manual(values = c("no" = "gray80", "yes" = "white"), name = "FS?", na.value = "white") + 
  new_scale_fill()

p4 = gheatmap(p3, df_gc, offset=1.55, width=0.05, colnames = FALSE, color=NA) +
  scale_fill_gradient(low = "#edc0c4", high = "#961415", name = "GC content", na.value = "white") + 
  new_scale_fill()

p5 = gheatmap(p4, df_prop, offset=1.60, width=0.05, colnames = FALSE, color=NA) +
  scale_fill_gradient(low = "#eddbc0", high = "#964f14", name = "TGA proportion", na.value = "white") + 
  new_scale_fill()

phylo = gheatmap(p3, df_fam, offset=0.2, width=0.05, colnames = FALSE, color=NA) +
  scale_fill_discrete(name = "family", na.value = "white") + 
  new_scale_fill()

p5

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\myco_tree_10_21_24.png", p5, units = "in", width = 8.5, height = 10.5, dpi = 600)
