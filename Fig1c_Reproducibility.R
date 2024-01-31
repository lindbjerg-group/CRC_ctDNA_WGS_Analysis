#Author: Amanda Frydendahl

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(psych)



#### Fig. 1c - reproducibility ####
# Set working directory
mywd <- "..." #Set your working dir
setwd(mywd)

# Read data from CSV file
dta <- read.csv2(file = file.path(mywd, "reproducibility.csv")) %>% 
  select(!PtID)

# Top panel - Create a tiled heatmap for ctDNA detection
dta_long <- dta %>% 
  select(LabA.detected, LabB.detected) %>% 
  mutate(ptID = row_number()) %>% 
  pivot_longer(cols = !ptID, 
               names_to = "lab", 
               values_to = "detection")

ggplot(dta_long, 
       aes(x = ptID, y = lab, fill = as.factor(detection), 
           col = "black", label = as.factor(detection))) +
  geom_tile() + 
  scale_color_manual("Detection", values = c("black")) +
  scale_fill_manual(values = c("#92BDDE", "#9F3C3A"), 
                    labels = c("ctDNA-negative", "ctDNA-positive")) +
  scale_y_discrete(limits = c("LabB.detected", "LabA.detected"), 
                   labels = c("LabB", "LabA")) +
  ylab("") +
  xlab("") +
  guides(fill = guide_legend(title = ""), color = "none") +  # Replace FALSE with "none"
  theme_classic() + 
  theme(axis.line = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank())


# Bottom panel - Scatter plot with estimated tumor fractions
ggplot(dta[dta$class != "none detected",], 
       aes(x = LabA_tumor_fraction, y = LabB_tumor_fraction, fill = class)) +
  geom_point(size = 5, shape = 21) + 
  geom_smooth(data = dta[dta$class == "both detected",],
              aes(x = LabA_tumor_fraction, y = LabB_tumor_fraction),
              col = "black",
              method = "lm", 
              show_guide = FALSE) +
  scale_fill_manual(values = c("#9F3C3A", "#949494"), 
                    labels = c("in both laboratories", "in one laboratory"),
                    guide = guide_legend(title = "ctDNA positive")) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  labs(x = "Tumor Fraction - Lab \"A\" ", y = "Tumor Fraction - Lab \"A\"") +
  theme_bw()


# Pearson correlation test
cor.test(dta[dta$class == "both detected",]$LabA_tumor_fraction,
         dta[dta$class == "both detected",]$LabB_tumor_fraction, method = "pearson")

# Cohen's Kappa estimate
cohen.kappa(x = cbind(dta$LabA.detected, dta$LabB.detected))

# McNemar test
dta_tab <- table(dta$LabA.detected, dta$LabB.detected)
mcnemar.test(dta_tab)





