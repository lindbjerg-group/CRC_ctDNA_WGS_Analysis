#Author: Amanda Frydendahl

#load packages
library(ggplot2)
library(dplyr)
library(scales)
library(survival)
library(survminer)


# Set working directory
mywd <- "..." #Set your working dir
setwd(mywd)

# Read data from CSV file
dta <- read.csv2(file = file.path(mywd, "ctDNA_results.csv"))


#### Fig. 4a - ctDNA before and after ACT ####
# Extract patients with both post-OP and post-ACT samples
# Categorize patients based on post-OP to post-ACT ctDNA status
dta_ <- dta %>% 
  filter(Sample_drawn_pre_relapse == 1) %>%
  group_by(PtID) %>%
  filter(all(c("post-op", "post-act") %in% time_category)) %>%
  filter(time_category %in% c("post-op", "post-act")) %>%
  mutate(label = case_when(
    all(ctDNA_status == "ND") ~ "neg-neg",
    all(ctDNA_status == "D") ~ "pos-pos",
    any(ctDNA_status == "ND") & any(ctDNA_status == "D") & ctDNA_status[time_category == "post-op"] == "D" ~ "pos-neg",
    any(ctDNA_status == "ND") & any(ctDNA_status == "D") & ctDNA_status[time_category == "post-op"] == "ND" ~ "neg-pos",
    TRUE ~ NA_character_
  )) %>%
  ungroup() %>% 
  mutate(
    TumorFraction  = as.numeric(TumorFraction)) %>%
  as.data.frame()


#Right panel: Plot post-OP and post-ACT tumor fractions for patients where both samples were collected
p1 <- ggplot(dta_, aes(x = time_category, y = TumorFraction, group = PtID)) + 
  geom_line(aes(col = as.factor(label), linetype = ifelse(label %in% c("pos-neg", "neg-neg"), "dashed", "solid")), 
            linewidth = 1) +  # Replace size with linewidth
  geom_point(aes(col = as.factor(label)), size = 3) + 
  scale_color_manual(values = c("#919292", "#773D3D", "#A1B8DA", "#BF6061"), 
                     name = "") +
  scale_linetype_manual(values = c("dashed", "solid"),
                        labels = c("No recurrence", "Recurrence"),
                        name = "") +
  scale_x_discrete("", limits = c("post-op", "post-act"), 
                   labels = c("Post-OP", "Post-ACT"), 
                   expand = c(0.2,0.2)) +
  scale_y_continuous("Tumor fraction in plasma", trans = 'log10',
                     breaks = c(0.1, 0.01, 0.001, 0),
                     labels = c(expression(10^-1), expression(10^-2), expression(10^-3), "0")) +
  annotation_logticks(
    base = 10,
    sides = "l") +
  theme_bw() + 
  theme(legend.position = "top")


#Left panel: Bar plot showing number of patients in each post-OP to post-ACT category
# Count the number of each combination of Relapse and label
dta_count <- dta_ %>% 
  select(PtID, Relapse, label) %>% 
  unique() %>% 
  group_by(Relapse, label) %>% 
  summarize(n = n(), .groups = 'drop')

#Make plot
p2 <- ggplot(dta_count, aes(x = label, y = n, fill = as.factor(label), col = "black", 
             linetype = as.factor(Relapse))) + 
  geom_bar(stat = "identity") + 
  scale_linetype_manual("", values = c("dashed", "solid")) +
  geom_text(aes(label = paste0("n = ", n)), position=position_stack(vjust=0.5)) + 
  labs(title = "") + 
  scale_fill_manual("", values = c("#919292", "#773D3D", "#A1B8DA", "#BF6061")) + 
  scale_color_manual("", values = c("black")) +
  scale_x_discrete(limits = c("neg-pos", "pos-pos", 
                              "pos-neg", "neg-neg")) +
  labs(x = "ctDNA and recurrence status", y = "# patients") +
  guides(color = "none") +
  theme_bw() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

#plot right and left panel side-by-side
ggarrange(p1, p2, common.legend = T)




#### Fig. 4b - Swimmer plot of recurrence patients with neg-neg ctDNA status ####

#Extract samples from the 5 recurrence patients that were ctDNA negative post-OP and post-ACT
dta_subset <- dta[dta$Sample_drawn_pre_relapse == 1,] %>%
  group_by(PtID) %>%
  filter(all(c("post-op", "post-act") %in% time_category) &
           all(ctDNA_status[time_category == "post-op"] == "ND") &
           all(ctDNA_status[time_category == "post-act"] == "ND") &
           all(Relapse == 1) &
           !is.na(ACT_start_days)) %>%
  ungroup() %>% 
  as.data.frame()

pts <- dta_subset$PtID %>% unique()

#sort dta_subset based on time to recurrence in descending order
dta_subset <- dta_subset[order(-dta_subset$time_to_recurrence_days), ]

ggplot(dta_subset, aes(x = time_to_recurrence_days/(365.25/12), y = reorder(PtID, time_to_recurrence_days))) +
  geom_linerange(aes(xmin = 0, xmax = time_to_recurrence_days/(365.25/12))) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = ACT_start_days/(365.25/12), 
                   y = as.factor(PtID),
                   xend = ACT_end_days/(365.25/12), 
                   yend = as.factor(PtID),
                   alpha = 0.2), color = "#76A6BB", size = 3) + 
  geom_point(aes(x = sample_timepoint/(365.25/12), 
                 y = PtID, 
                 fill = as.factor(ctDNA_status), 
                 stroke = 0.4), shape = 21, size = 5) +
  geom_point(aes(x = time_to_recurrence_days/(365.25/12), 
                 y = PtID, 
                 color = "red"), 
             shape = "|", size = 7) +
  scale_fill_manual("", values = c("black", "white"), 
                    labels = c("ctDNA-positive", "ctDNA-negative")) +
  scale_color_manual("", values = c("red"), 
                    labels = c("Recurrence")) +
  scale_alpha_continuous("", 
                         labels = c("ACT")) +
  labs(y = "", x = "Time since operation (months)") + 
  theme_pubr() + 
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank())
#Post-OP and post-OP plasma samples illustrated by red dashed line added manually


#### Fig. 4c - Post-ACT detection rate ####

# Extract post-ACT samples (drawn before relapse) and handle zero TumorFraction values
postACT <- dta[dta$time_category == "post-act" & dta$Sample_drawn_pre_relapse == 1,]
postACT <- postACT %>% mutate(
  TumorFraction  = as.numeric(ifelse(TumorFraction  == 0, NA, TumorFraction))
)

# Plot post-OP detection rate
postACT %>% 
  group_by(ctDNA_status, Relapse) %>% 
  summarize(n=n()) %>% 
  mutate(freq=n/sum(n)) %>% 
  ggplot(., aes(x = ctDNA_status, y = freq, 
                fill = factor(Relapse, levels = c("0", "1")))) + 
  geom_bar(stat = "identity", col = "black") + 
  geom_text(position = position_fill(vjust = 0.5),  
            aes(label = paste0("n = ",n, "\n(", round(freq,2)*100, "%)"))) +
  scale_fill_manual("", values = c("white", "#6D6D6D"), 
                    labels = c("No Recurrence", "Recurrence")) +
  scale_x_discrete(labels = c(paste0("ctDNA-positive\nn = ", nrow(postACT[postACT$ctDNA_status == "D",])), 
                              paste0("ctDNA-negative\nn = ", nrow(postACT[postACT$ctDNA_status == "ND",]))
  )) +
  labs(x = "", 
       y = "Detection rate") +
  theme_classic() + 
  theme(legend.position = "top")


#### Fig 4D - Post-ACT KM plot ####

#Calculate recurrence-free survival (RFS) from end of ACT and make ctDNA_status binary (0,1)
#For recurrence patients: time to recurrence
#For non recurrence patients: last scan
postACT <- postACT %>%
  mutate(RFS_days = ifelse(Relapse == 1, time_to_recurrence_days - ACT_end_days, 
                           last_scan - ACT_end_days)) %>%
  mutate(ctDNA_status = ifelse(ctDNA_status == "ND", 0, 1))


#Plot RFS as KM plot
ggsurvplot(survfit(Surv(RFS_days/(365.25/12), Relapse) ~ ctDNA_status, 
                   data = postACT, 
                   conf.type=c("log-log")),
           pval=T, 
           conf.int = T,
           risk.table=T, ncensor.plot = F,
           title= paste0("ctDNA post-ACT"),
           palette=c("#748BAA", "#952422"),
           xlim = c(0, 36),
           legend= "top", legend.title="", legend.labs=c("ctDNA-negative","ctDNA-positive"),
           ylab="Recurrence-free survival",
           risk.table.height=0.22, xlab="Time since operation (months)", 
           break.time.by=4, font.y = c(12, "bold", "black"))

#cox regression
c <- coxph(Surv(time = RFS_days/(365.25/12), event = Relapse)~ctDNA_status, data = postACT)
summary(c)


#### Fig. 4e - Cumulative incidence of ctDNA detection after end of treatment ####
#calculate sample timepoint since end of treatment (EOT)
#This is either from the post-OP sample or the post-ACT sample
dta_EOT <- dta %>% 
  mutate(sample_timepoint_EOT = ifelse(
    is.na(ACT_end_days), 
    sample_timepoint, 
    sample_timepoint - ACT_end_days)) %>% 
  filter(sample_timepoint_EOT > 0)

#summarize data
dta_EOT_fitdata <- dta_EOT %>%
  group_by(PtID) %>%
  summarise(
    Relapse = ifelse(any(Relapse == 1), 1, 0),
    ctDNA_D = ifelse(any(ctDNA_status == "D"), 1, 0),
    min_ctDNA_time = ifelse(any(ctDNA_status == "D"), min(sample_timepoint_EOT[ctDNA_status == "D"]), max(sample_timepoint_EOT))
  )

fit1 <- survfit(Surv(min_ctDNA_time/(365.25/12), ctDNA_D) ~ Relapse, 
                data = dta_EOT_fitdata)

#Make plot
ggsurvplot(fit1,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette=c("#71A2B6","#982649"),
           pval=T,
           pval.coord = c(10, 0.3),
           xlim = c(0, 36),
           ylab="Cumulative incidence of ctDNA detection",
           xlab="Time since end of treatment (months)",
           legend="top", legend.title="", legend.labs=c("No recurrence","Recurrence"),
           fun = "event",
           risk.table = T,
           break.time.by=2, font.y = c(12, "bold", "black"))


#### Cox regression using ctDNA as a time-dependent variable (surveillance analysis) ####
# Create RFS_days column, filter, and select columns
dta_EOT_ <- dta_EOT %>%
  mutate(RFS_days = ifelse(Relapse == 1, time_to_recurrence_days, last_scan)) %>%
  filter(sample_timepoint < RFS_days) %>%
  select(PtID, SampleID, sample_timepoint, ctDNA_status, Relapse, RFS_days)

# Count number of samples per PtID
n_samples <- dta_EOT_ %>%
  count(PtID) %>%
  filter(n > 1) %>%
  as.data.frame()

# Merge counts with the main dataset
dta_EOT_ <- merge(dta_EOT_, n_samples, by = "PtID")

# Only include patients with more than 1 sample
dta_EOT_ <- dta_EOT_[dta_EOT_$n > 1, ]

# Prepare data for Cox proportional hazards regression
mdat <- dta_EOT_ %>%
  select(PtID, Relapse, RFS_days) %>%
  unique()

# Merge data for time-dependent Cox regression
newdat <- tmerge(data1 = mdat, data2 = mdat, id = PtID, tstop = RFS_days)
newdat <- tmerge(data1 = newdat, data2 = dta_EOT_, id = PtID,
                 Recurrence = event(RFS_days, Relapse), 
                 ctDNA = tdc(sample_timepoint, ctDNA_status))

# Perform Cox proportional hazards regression
c <- coxph(data = newdat,
           Surv(time = tstart, time2 = tstop, event = Relapse) ~ factor(ctDNA, levels = c("ND", "D")))

# Summarize Cox regression results
summary(c)



#### Fig. 4f - Lead time analysis ####
#Extract recurrence patients with serially collected plasma samples
#Exclude patients where there is more than 6 months (183 days) between last plasma sample before recurrence and recurrence 
dta_ <- dta %>%
  filter(sample_timepoint > 0 & Relapse == 1) %>%
  group_by(PtID) %>%
  filter(min(Sample_to_relapse_days, na.rm = TRUE) <= 183)

#first positive plasma sample (mFIRST)
mFIRST <- dta_[dta_$ctDNA_status == "D",] %>% group_by(PtID) %>% 
  mutate(
    TTR_mFIRST = min(sample_timepoint)) %>% as.data.frame() %>% select(PtID, TTR_mFIRST) %>% unique()


#first positive plasma sample (mEOT)
tmp <- dta_ %>%
  filter((is.na(ACT_start_days) & (time_category == "post-op" | time_category == "surveillance")) |
           (!is.na(ACT_start_days) & (time_category == "post-act" | time_category == "surveillance")))

mEOT <- tmp[tmp$ctDNA_status == "D",] %>% group_by(PtID) %>% 
  mutate(
    TTR_mEOT = min(sample_timepoint)) %>% as.data.frame() %>% select(PtID, TTR_mEOT) %>% unique()

#merge mFIRST and mEOT dataframes to main data (dta_)
dta_ <- merge(dta_, mFIRST, by = "PtID", all.x = T)
dta_ <- merge(dta_, mEOT, by = "PtID", all.x = T)


#add overall lead time (OA_lead_time) to data
dta_$OA_lead_time <- dta_$time_to_recurrence_days - dta_$TTR_mFIRST

#add lead time from EOT to recurrence (EOT_lead_time)
dta_$EOT_lead_time <- dta_$time_to_recurrence_days - dta_$TTR_mEOT

#Find plot order by adding TTR_order as the longest lead time - either EOT_lead_time or OA_lead_time
dta_ <- dta_ %>% mutate(
  TTR_leadtime = 0, #radiological relapse is defined as time = 0
  TTR_order = ifelse(is.na(OA_lead_time) & is.na(EOT_lead_time), NA, ifelse(
    OA_lead_time > EOT_lead_time, OA_lead_time, EOT_lead_time)))



#Remove unnecessary columns
dta_ <- dta_ %>% select(PtID, time_to_recurrence_days, TTR_mFIRST, TTR_mEOT, 
                        OA_lead_time, EOT_lead_time, TTR_leadtime, TTR_order) %>% 
  unique()

#order PtID labels
labels <- dta_ %>% arrange(-TTR_order) %>% select(PtID) %>% unique()

#Make plot
dta_ %>%
  arrange(-TTR_order) %>%
  mutate(PtID = factor(PtID, levels = PtID)) %>%
  ggplot() +
  geom_linerange(aes(xmin = -OA_lead_time, xmax = TTR_leadtime, y = as.factor(PtID)), size = 0.7) +
  geom_linerange(aes(xmin = -EOT_lead_time, xmax = TTR_leadtime, y = as.factor(PtID)), size = 0.7) +
  geom_point(aes(x = -OA_lead_time, y = as.factor(PtID), color = "First ctDNA+"), size = 6) +
  geom_point(aes(x = -EOT_lead_time, y = as.factor(PtID), color = "First ctDNA+ post EOT"), size = 6) +
  geom_point(aes(x = TTR_leadtime, y = as.factor(PtID), color = "Recurrence"), size = 6) +
  scale_x_continuous(labels = seq(-36, 36, by = 4),
                     breaks = seq(-1100, 1100, by = 30.25 * 4)) +
  scale_color_manual("", values = c("#BFCDE0", "#235789", "#BF6060"),
                     labels = c("First ctDNA detection", "First ctDNA detection after EOT", "Radiological recurrence")) +
  labs(title = "", y = "Patient #", x = "Time to radiological recurrence (months)", color = "Event") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#calculate leadtime from first ctDNA detection to recurrence
summary(dta_[!is.na(dta_$EOT_lead_time),]$OA_lead_time/(365.25/12), na.rm=T)

#test for diff in TTR to TTR_mFIRST
wilcox.test(dta_$time_to_recurrence_days,
            dta_[!is.na(dta_$TTR_mFIRST),]$TTR_mFIRST)


#calculate leadtime from first ctDNA detection after EOT to recurrence
summary(dta_[!is.na(dta_$EOT_lead_time),]$EOT_lead_time/(365.25/12), na.rm=T)

#test for diff in TTR to TTR_mEOT 
wilcox.test(dta_$time_to_recurrence_days,
            dta_[!is.na(dta_$TTR_mEOT),]$TTR_mEOT)

