#Author: Amanda Frydendahl

#load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(survival)
library(survminer)


# Set working directory
mywd <- "..." #Set your working dir
setwd(mywd)

# Read data from CSV file
dta <- read.csv2(file = file.path(mywd, "ctDNA_results.csv"))


#### Pre-OP ####


#### Fig. 3a - Pre-OP detection rate ####
# Filter pre-OP samples and handle zero TumorFraction values
preOP <- dta[dta$time_category == "pre-op",]
preOP <- preOP %>% mutate(
  TumorFraction  = as.numeric(ifelse(TumorFraction  == 0, NA, TumorFraction))
)

#plot pre-OP detection rate
preOP %>% 
  group_by(ctDNA_status) %>% 
  summarize(n=n()) %>% 
  mutate(freq=n/sum(n)) %>% 
  ggplot(., aes(x = 1, y = freq, fill = factor(ctDNA_status, levels = c("ND", "D")))) + 
  geom_bar(stat = "identity", col = "black") + 
  geom_text(position = position_fill(vjust = 0.5),  
            aes(label = paste0("n = ",n, "\n(", round(freq,2)*100, "%)"))) +
  scale_fill_manual("", values = c("#748BAA", "#952422"), 
                    labels = c("ctDNA-negative", "ctDNA-positive")) +
  labs(x = paste0("n = ", nrow(preOP)), 
       y = "Detection rate") +
  theme_classic() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


#### Fig. 3b - Pre-OP tumor fractions ####
# Plot pre-OP tumor fractions for ctDNA positive samples
ggplot(preOP[preOP$ctDNA_status == "D",],
       aes(x = 1, y = TumorFraction)) + 
  geom_violin(width = 0.5, fill = "#952422", alpha = 0.5) +
  geom_boxplot(width= 0.2, fill = "#952422", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.08) +
  labs(x = paste0("n = ", nrow(preOP[!is.na(preOP$TumorFraction),])), 
       y = "Tumor fraction") +
  scale_y_log10() +
  annotation_logticks(
    base = 10,
    sides = "l") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.margin = margin(r = 350))


#### Fig. 3c - Total mutational burden and pre-OP ctDNA status ####
# Read data with total mutational burden and MSI status from CSV file
TMB <- read.csv2(file = file.path(mywd, "TMB_MSIstatus.csv"))

# Merge TMB data and ctDNA results
TMB_ctDNA <- merge(preOP, TMB, by = "PtID")

# Find median TMB of pre-OP ctDNA positive and negative patients
TMB_ctDNA[TMB_ctDNA$time_category == "pre-op" & TMB_ctDNA$ctDNA_status == "D",]$TMB %>% summary()
TMB_ctDNA[TMB_ctDNA$time_category == "pre-op" & TMB_ctDNA$ctDNA_status == "ND",]$TMB %>% summary()


# For synchronous tumors, use max TMB when plotting
TMB_ctDNA <- TMB_ctDNA %>% 
  select(PtID, time_category, ctDNA_status, TMB, MSI_cosmic) %>% 
  group_by(PtID) %>%
  mutate(TMB = ifelse(n() > 1, max(TMB, na.rm = TRUE), TMB)) %>% 
  unique()

# Plot TMB, stratified by pre-OP ctDNA status and colored by MSI status
TMB_ctDNA[TMB_ctDNA$time_category == "pre-op",] %>% 
  ggplot(aes(x = ctDNA_status , y = TMB)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, aes(col = MSI_cosmic), size = 2) + 
  #scale_y_continuous(limits = c(0, 50000)) +
  stat_compare_means(label.x = 1.3) +
  scale_color_manual("", values = c("#952422", "black")) +
  scale_x_discrete(label = c(paste0("ctDNA-positive\nn= ",
                                    nrow(TMB_ctDNA[TMB_ctDNA$time_category == "pre-op" & 
                                                TMB_ctDNA$ctDNA_status == "D",])), 
                             paste0("ctDNA-negative\nn= ",
                                    nrow(TMB_ctDNA[TMB_ctDNA$time_category == "pre-op" & 
                                                TMB_ctDNA$ctDNA_status == "ND",])))) +
  scale_y_continuous(trans = 'log10',
                     limits = c(500, 500000),
                     breaks = 10^(1:6),
                     labels = comma) +
  annotation_logticks(
    base = 10,
    sides = "l") + 
  labs(title = "", y = "Total mutational burden", x = "") +
  theme_bw()


#### Post-OP ####

#### Fig 3d - Post-OP detection rate #### 

# Filter post-OP samples and handle zero TumorFraction values
postOP <- dta[dta$time_category == "post-op",]
postOP <- postOP %>% mutate(
  TumorFraction  = as.numeric(ifelse(TumorFraction  == 0, NA, TumorFraction))
)

# Plot post-OP detection rate
postOP %>% 
  group_by(ctDNA_status, Relapse) %>% 
  summarize(n = n(), .groups = "drop_last") %>%  # Use .groups = "drop_last" to drop the last grouping variable
  mutate(freq = n / sum(n)) %>% 
  ggplot(., aes(x = ctDNA_status, y = freq, fill = factor(Relapse, levels = c("0", "1")))) + 
  geom_bar(stat = "identity", col = "black") + 
  geom_text(position = position_fill(vjust = 0.5), 
            aes(label = paste0("n = ", n, "\n(", round(freq, 2) * 100, "%)"))) +
  scale_fill_manual("", values = c("white", "#6D6D6D"), labels = c("No Recurrence", "Recurrence")) +
  scale_x_discrete(labels = c(paste0("ctDNA-positive\nn = ", sum(postOP$ctDNA_status == "D")), 
                              paste0("ctDNA-negative\nn = ", sum(postOP$ctDNA_status == "ND")))) +
  labs(x = "", y = "Detection rate") +
  theme_classic() + 
  theme(legend.position = "top")


#### Fig 3E - Post-OP KM plot ####

#Calculate recurrence-free survival and make ctDNA_status binary (0,1)
#For recurrence patients: time to recurrence
#For non recurrence patients: last scan
postOP <- postOP %>%
  mutate(RFS_days = ifelse(Relapse == 1, time_to_recurrence_days, last_scan)) %>% 
  mutate(ctDNA_status = ifelse(ctDNA_status == "ND", 0, 1))


#Plot RFS as KM plot
ggsurvplot(survfit(Surv(RFS_days/(365.25/12), Relapse) ~ ctDNA_status, 
                   data = postOP, 
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
c <- coxph(Surv(time = RFS_days/(365.25/12), event = Relapse)~ctDNA_status, data = postOP)
summary(c)



#### Fig 3f - Post-OP false negative rate ####

# Add postOP_category and call_category columns to the postOP dataframe
postOP <- postOP %>%
  mutate(postOP_category = ifelse(sample_timepoint > 14, "Post_day14", "Pre_day14")) %>%
  mutate(call_category = case_when(
    ctDNA_status == 1 & Relapse == 1 ~ "TruePositive",
    ctDNA_status == 1 & Relapse == 0 ~ "FalsePositive",
    ctDNA_status == 0 & Relapse == 1 ~ "FalseNegative",
    ctDNA_status == 0 & Relapse == 0 ~ "TrueNegative",
    TRUE ~ NA_character_
  ))

# Create a summary dataframe (dta) for further analysis
dta <- postOP %>%
  filter(call_category %in% c("FalseNegative", "TruePositive")) %>%
  group_by(postOP_category, call_category) %>%
  summarize(n = n(), .groups = "drop_last") %>%  # Drop the last grouping variable
  mutate(freq = n / sum(n) * 100)

# Calculate total number of samples drawn within 14 days post-OP
total_n_Pre_day14 <- dta %>%
  filter(postOP_category == "Pre_day14") %>%
  summarize(total_n = sum(n)) %>%
  pull(total_n) %>%
  as.integer()

# Calculate total number of samples drawn later than 14 days post-OP
total_n_Post_day14 <- dta %>%
  filter(postOP_category == "Post_day14") %>%
  summarize(total_n = sum(n)) %>%
  pull(total_n) %>%
  as.integer()

# Plot the false negative and true positive rate as piecharts
ggplot(dta, aes(x = "", y = freq,
                fill = factor(call_category, levels = c("TruePositive", "FalseNegative")))) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual("", values = c("#93BEDF", "#26547C"),
                    labels = c("True positive rate\n(sensitivity)", "False negative rate")) +
  geom_text(aes(y = freq, label = paste0(round(freq, 0), " %")),
            position = position_stack(vjust = 0.5)) +
  facet_wrap(~factor(postOP_category, levels = c("Pre_day14", "Post_day14")),
             labeller = as_labeller(c("Pre_day14" = paste0("Expected positive samples drawn\nwithin 14 days post-OP, n = ",
                                                           total_n_Pre_day14),
                                      "Post_day14" = paste0("Expected positive samples drawn\nlater than 14 days post-OP, n = ",
                                                            total_n_Post_day14))),
             dir = "v") +
  theme_void() +
  theme(legend.position = "bottom")
