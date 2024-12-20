# Develop ESP indicator for juvenile snow crab energetic condition:
#annual mean % DWt of hepatopancreas

# Erin Fedewa

#2025 improvements: Use model-based bayesian output for annual mean to control
#for julian date and crab size? 

# load ----
library(tidyverse)
library(ggridges)
library(akgfmaps)
library(patchwork)
library(viridis)

#######################################

sc_condition <- read.csv("./output/condition_master.csv")

#data wrangling 
sc_condition %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         lme != "NBS",
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year) %>%
  summarise(avg_condition = mean(Perc_DWT, na.rm=T)) -> cond

#Write output
missing <- data.frame(year = 2020)

cond %>%
  bind_rows(missing) %>%
  arrange(year) %>%
  write.csv(file="./output/2024_esp_opilio_condition.csv")

#####################################################
#Plots 
sc_condition %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  group_by(year, lme) %>%
  summarise(sd = sd(Perc_DWT, na.rm = TRUE),
            cond = mean(Perc_DWT, na.rm=T)) %>%
  bind_rows(missing) %>%
  arrange(year) %>%
  filter(lme != "NA") %>%
  mutate(lme = recode(lme, EBS = "Eastern Bering Sea",
                      NBS = "Northern Bering Sea")) -> plot

#plot annual mean as bar plot 
plot %>% 
  ggplot(aes(as.factor(year), cond, fill=year)) +
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin = cond-sd, ymax = cond+sd), width = 0.2, color="gray45") +
  facet_wrap(~lme) +
  theme_bw() +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 11, axis_text_size =10) +
  labs(y = "Snow Crab Condition (% DWT)", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(yintercept=22.6, color="red") + #adding prelim threshold from starvation lab experiment
  geom_rect(aes(xmin=0, xmax=Inf, ymin=22.6, ymax=22.6 + 2.9), fill="red", alpha = 0.05) + 
  geom_rect(aes(xmin=0, xmax=Inf, ymin=22.6-2.9, ymax=22.6), fill="red", alpha = 0.05) +
  geom_vline(data = subset(plot, lme == "Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(plot, lme == "Eastern Bering Sea"), aes(x = .7, y=44, label = "Mid-collapse"),
            size = 2.5, color = "#D55E00") +
  geom_text(data = subset(plot, lme == "Eastern Bering Sea"), aes(x = 3, y=44, label = "Post-collapse"),
            size = 2.5, color = "#084594") 
ggsave("./figures/condition_starvation.png", height = 6, width = 7, units = "in", dpi = 300)

#plot annual mean as dot plot
plot %>%
  ggplot(aes(as.factor(year), cond)) +
  geom_point(size=2) + 
  labs(y = "Snow Crab Condition", x = "") +
  theme(axis.text=element_text(size=12)) +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = mean(cond, na.rm=TRUE)), linetype = 5)+
  geom_hline(aes(yintercept = quantile(cond, .10, na.rm=TRUE)), linetype = 3)+
  geom_hline(aes(yintercept = quantile(cond, .90, na.rm=TRUE)), linetype = 3)+
  theme_bw()

#size comp of crab sampled
sc_condition %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         lme != "NBS",
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1) %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  filter(lme != "NA", sex != "NA") %>% #one crab collected outside the sampling design
  group_by(year) %>%
  ggplot() +
  geom_density(aes(x=cw, fill=Sex), position = "stack", binwidth = 2) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_grid(lme~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")

############################################################################
#Improvements for 2025- exploring a model controlling for seasonality and 
  #crab size 

