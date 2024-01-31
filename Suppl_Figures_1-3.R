#Author: Amanda Frydendahl

#load packages
library(ggplot2)
library(dplyr)
library(ggpubr)
library(survival)
library(survminer)


# Set working directory
mywd <- "..." #Set your working dir
setwd(mywd)

# Read data from CSV file
dta <- read.csv2(file = file.path(mywd, "ctDNA_results.csv"))



#### Supplementary Figure 1 - swimmer plot of all patients ####

# Reorder PtID based on last_scan and Relapse
dta <- dta %>%
  arrange(Relapse, last_scan) %>%
  mutate(PtID = factor(PtID, levels = unique(PtID)))


ggplot(dta, aes(x = last_scan/(365.25/12), y = PtID)) +
  #add follow-up time (last scan)
  geom_linerange(aes(xmin = 0, xmax = ifelse(last_scan/(365.25/12) < 40, 
                                             last_scan/(365.25/12), 40))) +
  #add adjuvant chemotherapy
  geom_segment(aes(x = ACT_start_days/(365.25/12), 
                   y = as.factor(PtID),
                   xend = ACT_end_days/(365.25/12), 
                   yend = as.factor(PtID),
                   alpha = 0.8), color = "#76A6BB", size = 1.5) + 
  #add plasma samples
  geom_point(aes(x = sample_timepoint/(365.25/12), 
                 y = PtID, 
                 fill = as.factor(ctDNA_status), 
                 stroke = 0.4), 
             shape = 21, size = 2) +
  #add intervention start day
  geom_point(aes(x = Intervention_start_days/(365.25/12), 
                 y = PtID, 
                 color = "intervention"),
             size = 1, shape = 4, stroke = 1) +
  #add recurrence
  geom_point(aes(x = time_to_recurrence_days/(365.25/12), 
                 y = PtID, 
                 color = "recurrence"),
             shape = "I", size = 2, stroke = 3) +
  geom_vline(xintercept = 0) +
  scale_fill_manual("", values = c("black", "white"), 
                    labels = c("ctDNA-positive", "ctDNA-negative")) +
  scale_color_manual("", values = c("#FAB95D", "#AC1917"), 
                     labels = c("Intervention start", "Recurrence")) +
  scale_alpha_continuous("", 
                         labels = c("Adjuvant chemotherapy")) +
  labs(y = "", x = "Time since operation (months)") + 
  xlim(0,40) +
  theme_bw() + 
  theme(axis.line = element_line(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), 
        legend.position = "top", 
        axis.text.y = element_text(size = 8))
# change keys for Intervention start and recurrence manually to fit true annotation


#### Supplementary Figure 2a - KM plot for "late" post-OP samples ####
#Extract late post-OP samples (drawn later than 2 weeks post-OP, but before ACT)
postOP_late <- dta[dta$time_category == "post-op" & dta$sample_timepoint > 14,] %>%
  mutate(RFS_days = ifelse(Relapse == 1, time_to_recurrence_days, last_scan)) %>% 
  mutate(ctDNA_status = ifelse(ctDNA_status == "ND", 0, 1))


#Plot RFS as KM plot
ggsurvplot(survfit(Surv(RFS_days/(365.25/12), Relapse) ~ ctDNA_status, 
                   data = postOP_late, 
                   conf.type=c("log-log")),
           pval=T, 
           conf.int = T,
           risk.table=T, ncensor.plot = F,
           title= paste0("ctDNA post-OP"),
           palette=c("#748BAA", "#952422"),
           xlim = c(0, 36),
           legend= "top", legend.title="", legend.labs=c("ctDNA-negative","ctDNA-positive"),
           ylab="Recurrence-free survival",
           risk.table.height=0.22, xlab="Time since operation (months)", 
           break.time.by=4, font.y = c(12, "bold", "black"))

#cox regression
c <- coxph(Surv(time = RFS_days/(365.25/12), event = Relapse)~ctDNA_status, data = postOP_late)
summary(c)


#### Supplementary Figure 2b - Detection rate for "late" post-OP samples ####
# Plot post-OP detection rate
postOP_late <- dta[dta$time_category == "post-op" & dta$sample_timepoint > 14,]

postOP_late %>% 
  group_by(ctDNA_status, Relapse) %>% 
  summarize(n=n()) %>% 
  mutate(freq=n/sum(n)) %>% 
  ggplot(., aes(x = as.factor(ctDNA_status), y = freq, 
                fill = factor(Relapse, levels = c("0", "1")))) + 
  geom_bar(stat = "identity", col = "black") + 
  geom_text(position = position_fill(vjust = 0.5),  
            aes(label = paste0("n = ",n, "\n(", round(freq,2)*100, "%)"))) +
  scale_fill_manual("", values = c("white", "#6D6D6D"), 
                    labels = c("No Recurrence", "Recurrence")) +
  scale_x_discrete(labels = c(paste0("ctDNA-positive\nn = ", nrow(postOP_late[postOP_late$ctDNA_status == "D",])), 
                              paste0("ctDNA-negative\nn = ", nrow(postOP_late[postOP_late$ctDNA_status == "ND",]))
  )) +
  labs(x = "", 
       y = "Detection rate") +
  theme_classic()



#### Supplementary Figure 3a - False positives, swimmer plot ####

#Identify false positive, true positive and negative samples and extract results from those patients
#Only include samples drawn within the follow-up time
dta_FP <- dta %>%
  filter(Relapse == 0 & sample_timepoint < last_scan) %>%
  mutate(False_Positive = case_when(
    is.na(ACT_end_days) ~ sample_timepoint > 0 & ctDNA_status == "D",
    sample_timepoint > ACT_end_days & ctDNA_status == "D" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(False_Positive = case_when(
    False_Positive ~ "False positive",
    ctDNA_status == "D" ~ "True positive",
    TRUE ~ "Negative"
  ))

#Patients with false positive samples
pts <- dta_FP[dta_FP$False_Positive == "False positive",]$PtID %>% unique()


#sort dta_subset based on time to recurrence in descending order

dta_subset <- dta_FP[dta_FP$PtID %in% pts,]
dta_subset <- dta_subset[order(-dta_subset$last_scan), ]

ggplot(dta_subset, aes(x = last_scan/(365.25/12), y = reorder(PtID, last_scan))) +
  geom_linerange(aes(xmin = 0, xmax = ifelse(last_scan/(365.25/12) < 40, 
                                             last_scan/(365.25/12), 40))) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = ACT_start_days/(365.25/12), 
                   y = as.factor(PtID),
                   xend = ACT_end_days/(365.25/12), 
                   yend = as.factor(PtID),
                   alpha = 0.2), color = "#76A6BB", size = 3) + 
  geom_point(aes(x = sample_timepoint/(365.25/12), 
                 y = PtID, 
                 fill = as.factor(False_Positive), 
                 stroke = 0.4), shape = 21, size = 5) +
  scale_fill_manual("", values = c("#AC1919", "white", "black"), 
                    labels = c("False Positive", "Negative", "True Positive")) +
  scale_alpha_continuous("", 
                         labels = c("ACT")) +
  scale_x_continuous(lim = c(NA,40), breaks = seq(0, max(dta_subset$last_scan/(365.25/12)), by = 4)) +
  labs(y = "", x = "Time since operation (months)") + 
  theme_pubr() + 
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank())

# Count the number of false positive surveillance samples
num_false_positive_surveillance <- nrow(dta_FP[dta_FP$False_Positive == "False positive",])

#### Supplementary Figure 3b - False positives, tumor fraction ####

# Extract true positive samples (Pre-OP samples)
dta_TP_preOP <- dta %>%
  filter(time_category == "pre-op" & ctDNA_status == "D") %>%
  mutate(False_Positive = "True Positive PreOP")

# Count the number of true positive pre-op samples
num_true_positive_preop <- nrow(dta_TP_preOP)

# Extract true positive samples (surveillance samples)
dta_TP_surv <- dta %>%
  filter(Relapse == 1 & sample_timepoint < time_to_recurrence_days) %>%
  mutate(False_Positive = case_when(
    is.na(ACT_end_days) ~ sample_timepoint > 0 & ctDNA_status == "D",
    sample_timepoint > ACT_end_days & ctDNA_status == "D" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(False_Positive = case_when(
    False_Positive ~ "True positive surveillance",
    ctDNA_status == "D" ~ "False positive",
    TRUE ~ "Negative"
  )) %>%
  filter(False_Positive == "True positive surveillance")

# Count the number of true positive surveillance samples
num_true_positive_surveillance <- nrow(dta_TP_surv)

#merge dataframe with true positive and false positive samples
dta_total <- rbind(dta_FP[dta_FP$False_Positive == "False positive",], 
      dta_TP_preOP, 
      dta_TP_surv)


#Plot
#Define which groups you want to compare
my_comparisons <- list( c("False positive", "True positive surveillance"), 
                        c("False positive", "True Positive PreOP"), 
                        c("True positive surveillance", "True Positive PreOP") )

ggplot(dta_total, aes(x = False_Positive, y = as.numeric(TumorFraction))) + 
  geom_boxplot(width = 0.3, outlier.shape = NA, size = 0.7) +
  geom_jitter(aes(col = False_Positive), width = 0.15, size = 2) + 
  scale_x_discrete(limits = c("False positive", "True positive surveillance", "True Positive PreOP"), 
                   labels = c(paste0("Surveillance samples\nfrom patients without\nrecurrence\nn = ", 
                                     num_false_positive_surveillance), 
                              paste0("Surveillance samples\nfrom recurrence\npatients\nn = ",
                                     num_true_positive_surveillance), 
                              paste0("Pre-OP\nsamples\nn = ", 
                                     num_true_positive_preop))) +
  scale_color_manual("", values = c("False positive" = "#AC1919", "True Positive PreOP" = "black", 
                                    "True positive surveillance" = "black"), 
                     labels = c("False positive", "True positive", "True positive")) +
  scale_y_continuous(trans = 'log10', 
                     breaks = 10^(-6:1)) +
  annotation_logticks(
    base = 10,
    sides = "l") +
  labs(title = "", x = "", y = "Estimated tumor fractions") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox") +
  theme_bw() + 
  theme(axis.text.x  = element_text(size = 8), 
        legend.position = c(0.8, 0.15),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
#remove one of the keys for "True positive" manually
