library(tidyverse)
library(ggprism)
library(ggbreak)
library(readxl)
library(ggpubr)
library(rstatix)
library(ggsignif)

df_raw = read_excel("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Western Blot\\FS efficiencies_5_14_24.xlsx", sheet = "raw")
df_an_ms = read_excel("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Western Blot\\FS efficiencies_5_14_24.xlsx", sheet = "analyzed") %>%
  filter(species == "Ms")
df_an_bs = read_excel("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Western Blot\\FS efficiencies_8_8_24.xlsx", sheet = "analyzed")

df_an = rbind(df_an_ms, df_an_bs)

df_an %>%
  group_by(strain) %>%
  summarize(mean = mean(fs_efficiency), sd = sd(fs_efficiency)) 

level_order = c('CP127', 'CP271', 'CP253', 'CP267') 

df_an_filt = df_an %>%
  filter(strain != "CP252") %>% 
  filter(strain != "CP109")

p = ggplot(df_an_filt, aes(x = factor(strain, level = level_order), y = fs_efficiency)) +
  geom_bar(stat = "summary", fill = "gray80", color = "black") +
  geom_jitter(height = 0.5, width = 0.2) +
  scale_y_continuous(expand= c(0,0), limits = c(0, 65)) +
  labs(x = "", y = "Frameshifting efficiency (%)") +
  theme_prism() +
  theme(text = element_text(size = 18)) +
  geom_errorbar(stat="summary",  
                fun.ymin=function(x) {mean(x)-sd(x)/sqrt(length(x))},  
                fun.ymax=function(x) {mean(x)+sd(x)/sqrt(length(x))},
                width=0.2)+
  geom_signif(comparisons = list(c("CP127", "CP253")), 
              test = "t.test", test.args = list(var.equal = FALSE), textsize = 5)

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures\\fs_efficiency_bars_8_8_24.png", p, width = 5, height = 3.75, dpi = 600, units = "in")
