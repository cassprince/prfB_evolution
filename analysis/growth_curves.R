# Load packages and data.

library(tidyverse)
library(readxl)
library(ggprism)

OD600_30 = read_excel("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Plate Reader\\CRP_24H_30dC_Aug_06_2024.xlsx", sheet = "OD600", col_types = "numeric")

OD600_37 = read_excel("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Data\\Plate Reader\\CRP_24H_37dC_Aug_06_2024.xlsx", sheet = "OD600", col_types = "numeric")

# Select wells for strains at 30 dC and pivot.

df_30 = data.frame(cbind(Time = OD600_30$Time,
                      CP201.1_0.05 = rowMeans(select(OD600_30, B2,B3,B4)),
                      CP201.2_0.05 = rowMeans(select(OD600_30, C2,C3,C4)),
                      CP201.3_0.05 = rowMeans(select(OD600_30, D2,D3,D4)),
                      CP202.1_0.05 = rowMeans(select(OD600_30, B5,B6,B7)),
                      CP202.2_0.05 = rowMeans(select(OD600_30, C5,C6,C7)),
                      CP202.3_0.05 = rowMeans(select(OD600_30, C5,C6,C7)),
                      CP201.1_0.005 = rowMeans(select(OD600_30, E2,E3,E4)),
                      CP201.2_0.005 = rowMeans(select(OD600_30, F2,F3,F4)),
                      CP201.3_0.005 = rowMeans(select(OD600_30, G2,G3,G4)),
                      CP202.1_0.005 = rowMeans(select(OD600_30, E5,E6,E7)),
                      CP202.2_0.005 = rowMeans(select(OD600_30, F5,F6,F7)),
                      CP202.3_0.005 = rowMeans(select(OD600_30, G5,G6,G7))
))

df_30_pivot = df_30 %>%
  pivot_longer(cols = colnames(df_30[-1]), names_to= "name", values_to = "value") 

# Make new columns for info on strain, dilution, and replicate.
df_30_names = df_30_pivot %>%
  mutate(strain = sub("\\..*", "", df_30_pivot$name)) %>%
  mutate(rep = gsub(".*\\.(.+)_.*", "\\1", df_30_pivot$name)) %>%
  mutate(dilution = str_extract(name, "[^_]+$")) %>%
  mutate(new_name = paste0(as.character(strain), "_", as.character(dilution)))

# Summarize data (mean and standard deviation) and plot data from 0.005 starting OD.
summ_30 = df_30_names %>%
  group_by(new_name, Time) %>%
  mutate(mean = mean(value), sd = sd(value)) %>%
  filter(dilution == "0.005")

plot_30 = ggplot(summ_30, aes(x = Time, y = mean, color = strain)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.4,
                position=position_dodge(0.05)) +
  theme_prism() +
  labs(x = "Time (h)", y = "log10(OD600)") +
  scale_x_continuous(limits = c(0,17)) +
  scale_y_continuous(trans='log10', n.breaks = 8) +
  theme(text=element_text(size = 20)) + 
  scale_color_manual(values = c("CP201" = "gray50", "CP202" = "#961415"))

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures//prfB_curve_30dc_8_7_24.png", plot_30, dpi = 600, width = 7, height = 5, units = "in")

# Select wells for strains at 37 dC and pivot.

df_37 = data.frame(cbind(Time = OD600_37$Time,
                         CP201.1_0.05 = rowMeans(select(OD600_37, B2,B3,B4)),
                         CP201.2_0.05 = rowMeans(select(OD600_37, C2,C3,C4)),
                         CP201.3_0.05 = rowMeans(select(OD600_37, D2,D3,D4)),
                         CP202.1_0.05 = rowMeans(select(OD600_37, B5,B6,B7)),
                         CP202.2_0.05 = rowMeans(select(OD600_37, C5,C6,C7)),
                         CP202.3_0.05 = rowMeans(select(OD600_37, C5,C6,C7)),
                         CP201.1_0.005 = rowMeans(select(OD600_37, E2,E3,E4)),
                         CP201.2_0.005 = rowMeans(select(OD600_37, F2,F3,F4)),
                         CP201.3_0.005 = rowMeans(select(OD600_37, G2,G3,G4)),
                         CP202.1_0.005 = rowMeans(select(OD600_37, E5,E6,E7)),
                         CP202.2_0.005 = rowMeans(select(OD600_37, F5,F6,F7)),
                         CP202.3_0.005 = rowMeans(select(OD600_37, G5,G6,G7))
))

df_37_pivot = df_37 %>%
  pivot_longer(cols = colnames(df_37[-1]), names_to= "name", values_to = "value") 

# Make new columns for info on strain, dilution, and replicate.
df_37_names = df_37_pivot %>%
  mutate(strain = sub("\\..*", "", df_37_pivot$name)) %>%
  mutate(rep = gsub(".*\\.(.+)_.*", "\\1", df_37_pivot$name)) %>%
  mutate(dilution = str_extract(name, "[^_]+$")) %>%
  mutate(new_name = paste0(as.character(strain), "_", as.character(dilution)))

# Summarize data (mean and standard deviation) and plot data from 0.005 starting OD.
summ_37 = df_37_names %>%
  group_by(new_name, Time) %>%
  mutate(mean = mean(value), sd = sd(value)) %>%
  filter(dilution == "0.005")

plot_37 = ggplot(summ_37, aes(x = Time, y = mean, color = strain)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  theme_prism() +
  labs(x = "Time (h)", y = "log10(OD600)") +
  scale_x_continuous(limits = c(0,10)) +
  scale_y_continuous(trans='log10', n.breaks = 8) +
  theme(text=element_text(size = 20)) + 
  scale_color_manual(values = c("CP201" = "gray50", "CP202" = "#961415"))

ggsave("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\prfB\\Figures//prfB_curve_37dc_8_7_24.png", plot_37, dpi = 600, width = 7, height = 5, units = "in")
