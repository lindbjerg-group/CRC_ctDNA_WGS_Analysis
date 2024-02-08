library(ggplot2)
library(ggpubr)
library(argparse)
library(cowplot)
library(stringr)
library(ragg)

COLORS6 <- c("#2EBAED", "#000000", "#DE1C14","#D4D2D2", "#ADCC54", "#F0D0CE")

INDEL_CATEGORIES <- tibble::tibble(
  "muttype" = c(
    rep("C_deletion", 6), rep("T_deletion", 6), rep("C_insertion", 6),
    rep("T_insertion", 6), rep("2bp_deletion", 6), rep("3bp_deletion", 6),
    rep("4bp_deletion", 6), rep("5+bp_deletion", 6), rep("2bp_insertion", 6),
    rep("3bp_insertion", 6), rep("4bp_insertion", 6), rep("5+bp_insertion", 6),
    rep("2bp_deletion_with_microhomology", 1), rep("3bp_deletion_with_microhomology", 2),
    rep("4bp_deletion_with_microhomology", 3), rep("5+bp_deletion_with_microhomology", 5)
  ),
  "muttype_sub" = c(
    rep(c(seq_len(5), "6+"), 2),
    rep(c(0:4, "5+"), 2),
    rep(c(seq_len(5), "6+"), 4),
    rep(c(0:4, "5+"), 4), 1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"
  )
)
INDEL_CATEGORIES = paste0(INDEL_CATEGORIES$muttype, "_", INDEL_CATEGORIES$muttype_sub)

plot_96_profile = function (mut_matrix, colors = NA, ymax = NA, condensed = FALSE, relative = TRUE) 
{
  freq <- full_context <- substitution <- context <- NULL
  if (is.na(colors)) {
    colors <- COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  
  if(relative) {
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  } else {
    norm_mut_matrix = mut_matrix
  }
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
    dplyr::mutate(substitution = stringr::str_replace(full_context, 
                                                      "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
                                                                                                                 "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", 
                        values_to = "freq") %>% dplyr::mutate(sample = factor(sample, 
                                                                              levels = unique(sample)))
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution, 
                                width = width)) + geom_bar(stat = "identity", colour = "black", 
                                                           size = 0.2) + scale_fill_manual(values = colors) + facet_grid(sample ~ 
                                                                                                                           substitution, scales = "free_y") + ylab("Relative contribution") + 
    guides(fill = "none") + theme_bw() + theme(axis.title.y = element_text(size = 12, 
                                                                          vjust = 1), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
                                              axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5), 
                                              strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9), 
                                              panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, 
                                                                                                           "lines"))
  if(!is.na(ymax)) {
    plot = plot + coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, ymax, 0.1))
  }
  
  return(plot)
}

#' Calculate relative contributions
#'
relativeContributions = function(contributions) {
  contributions_long_rel = contributions
  for(sample in unique(contributions_long_rel$Samples)) {
    rows_in = contributions_long_rel$Samples == sample
    contributions_long_rel$value[rows_in] = contributions_long_rel$value[rows_in] / sum(contributions_long_rel$value[rows_in]) 
  }
  return(contributions_long_rel)
}

#' Group signatures
#'
groupSigs = function(contribution_long, groupSignatures, groupName) {
  cnames = as.character(unique(contribution_long$variable))
  contribution_long$variable= as.character(contribution_long$variable)
  if(sum(contribution_long$variable %in% groupSignatures) > 0) {
    contribution_long$variable[contribution_long$variable  %in% groupSignatures] = groupName
    
    temp = aggregate(value ~ Samples + variable, data = contribution_long[contribution_long$variable == groupName,], sum)
    
    contribution_long = rbind(contribution_long[contribution_long$variable != groupName,], temp)
    
    cnames = cnames[!(cnames %in% groupSignatures)]
    cnames = c(cnames, groupName)
    contribution_long$variable  = factor(contribution_long$variable , levels = cnames)
  }
  return(contribution_long)
}

cos_sim <- function(x, y) {
  res <- x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res <- as.numeric(res)
  return(res)
}

indel_names_COSMIC_to_MutPatterns = function(contexts) {
  
  indel_components = str_split_fixed(contexts, ":", 4)
  indel_components[,2] = str_replace(indel_components[,2], "Del", "deletion")
  indel_components[,2] = str_replace(indel_components[,2], "Ins", "insertion")
  indel_components[indel_components[,1] == 5,1] = "5+"
  indel_components[,1] = paste0(indel_components[,1], "bp")
  indel_components[indel_components[,4] == 5,4] = "5+"
  
  new_names = paste0(indel_components[,1], "_", indel_components[,3], "_", indel_components[,2], "_", indel_components[,4])
  new_names = str_replace(new_names, "1bp_", "")
  new_names = str_replace(new_names, "_R", "")
  new_names = str_replace(new_names, "_M_deletion_", "_deletion_with_microhomology_")
  
  new_names = str_replace(new_names, "deletion_5\\+", "deletion_6+")
  new_names = str_replace(new_names, "deletion_4", "deletion_5")
  new_names = str_replace(new_names, "deletion_3", "deletion_4")
  new_names = str_replace(new_names, "deletion_2", "deletion_3")
  new_names = str_replace(new_names, "deletion_1", "deletion_2")
  new_names = str_replace(new_names, "deletion_0", "deletion_1")
  
  #new_names[order(new_names)][new_names[order(new_names)] == mut_pat[order(mut_pat)]]
  #new_names[order(new_names)][new_names[order(new_names)] != mut_pat[order(mut_pat)]]
  #mut_pat[order(mut_pat)][new_names[order(new_names)] != mut_pat[order(mut_pat)]]
  
  return(new_names)
}

#makesignaturegroups <- function(novelsignatures=""){
  signature_groups <- list()
  signature_groups[["Ageing"]] <- c("SBS1")
  signature_groups[["SBS5"]] <- c("SBS5")
  signature_groups[["Slippage"]] <- c("ID1","ID2")
  signature_groups[["APOBEC"]] <- c("SBS2", "SBS13")
  signature_groups[["HRD"]] <- c("SBS3")
  signature_groups[["NHEJ"]] <- c("ID6","ID8")
  signature_groups[["MMRd"]] <- c("SBS6","SBS15","SBS20", "SBS21","SBS26","SBS44","ID7", "DBS7", "DBS10")
  signature_groups[["Smoking"]] <- c("SBS4", "ID3","DBS2")
  signature_groups[["Smoking(Bladder)"]] <- c("SBS92")
  signature_groups[["Tobacco"]] <- c("SBS29")
  signature_groups[["UV_SBS7"]] <- c("SBS7a","SBS7b","SBS7c","SBS7d","ID13", "DBS1", "DBS1")
  signature_groups[["UV_SBS38"]] <- c("SBS38")
  signature_groups[["POLE"]] <- c("SBS10a","SBS10b", "DBS3")
  signature_groups[["POLE+MMRd"]] <- c("SBS14")
  signature_groups[["E.coli"]] <- c("SBS88","ID18")
  signature_groups[["TMZ"]] <- c("SBS11")
  signature_groups[["Platinum"]] <- c("SBS31", "SBS35", "DBS5")
  signature_groups[["Other"]] <- " "
  signature_groups[["ROS"]] <- c("SBS18")
  signature_groups[["MUTYH"]] <- c("SBS36")
  signature_groups[["PolD1"]] <- c("SBS10c","SBS10d")
  #signature_groups[["Novel signature"]] <- novelsignatures
  signature_groups[["Mutagen exposure"]] <- c("SBS22", "SBS24","SBS25",  "SBS32", "SBS42", "SBS86", "SBS87", "SBS90")
  signature_groups[["Somatic hypermutation"]] <- c("SBS84","SBS85", "SBS9")
  signature_groups[["TOP2A"]] <- c("ID17")
  #signature_groups[["BER(NTHL1)"]] <- c("SBS30")
  signature_groups[["SBS30"]] <- c("SBS30")
  signature_groups[["SBS17"]] <- c("SBS17a", "SBS17b")
  signature_groups[["SBS93"]] <- c("SBS93")
  signature_groups[["SBS57"]] <- c("SBS57")
  signature_groups[["SBS37"]] <- c("SBS37")
#  signatureGroups <- rep(names(signature_groups),sapply(signature_groups,length))
#  names(signatureGroups) <- unlist(signature_groups)
#  signatureGroups
#}
  
#makemutsigcolors <- function(){
  mut.sig.colors <- c("Ageing"='#3B99D4',
                      "SBS5"='#b4dbed',
                      "Slippage"='#73b3d1',
                      "APOBEC"='#8ED14B',
                      "HRD"='#F0B849',
                      "NHEJ"='#F06B49',
                      "MMRd"='#D92B45',
                      "Smoking"='#B221E8',
                      "Smoking(Bladder)"='#ba6dd6',
                      "Tobacco"='#5a157a',
                      "UV_SBS7"='#FFEC00',
                      "UV_SBS38"='#e6d400',
                      "POLE"='#19413E',
                      "POLE+MMRd"='#28857F',
                      "E.coli"='#9de3cb',
                      "TMZ"='#076B82',
                      "Platinum"='#7D7C7C',
                      "Other"='#E8E8E8',
                      "ROS"='#aba7c7',
                      "MUTYH"='#A7C7C5',
                      "PolD1"='#1bbfb5',
                      "Novel signature"='#f2b8f0',
                      "Mutagen exposure"='#aef0a1',
                      "Somatic hypermutation"='#f054b9',
                      "TOP2A"='#decfad',
                      #"BER(NTHL1)"='#bda46c',
                      "SBS30"='#bda46c',
                      "SBS17"='#786846',
                      "SBS93"="#F4A460",
                      "SBS57"="#C71585",
                      "SBS37"="#1E90FF")
#  mut.sig.colors
#}

colors_real_artefact = c("#990011FF",
                         "#FCF6F5FF")
names(colors_real_artefact) = c("Artefact", "Real")

  #to ensure same color for same signatures in all plots
  colors_list = list()
  colors_list[["SBS1"]] = "#3B99D4"
  colors_list[["SBS2"]] = "#b0d14b"
  colors_list[["SBS3"]] = "#F06B49"
  colors_list[["SBS4"]] = "#B221E8"
  colors_list[["SBS5"]] = "#8c564b"
  colors_list[["SBS6"]] = "#d62728"
  colors_list[["SBS7a"]] = "#FFEC00"
  colors_list[["SBS7b"]] = "#e8ff00"
  colors_list[["SBS8"]] = "#2ca02c"
  colors_list[["SBS9"]] = "#17becf"
  colors_list[["SBS10a"]] = "#19413E"
  colors_list[["SBS10b"]] = "#276661"
  colors_list[["SBS11"]] = "#076B82"
  colors_list[["SBS13"]] = "#8ED14B"
  colors_list[["SBS14"]] = "#8A2BE2"
  colors_list[["SBS15"]] = "#D2691E"
  colors_list[["SBS17a"]] = "#006400"
  colors_list[["SBS17b"]] = "#FF1493"
  colors_list[["SBS18"]] = "#1E90FF"
  colors_list[["SBS20"]] = "#C71585"
  colors_list[["SBS21"]] = "#FA8072"
  colors_list[["SBS26"]] = "#FF00FF"
  colors_list[["SBS29"]] = "#50136D"  
  colors_list[["SBS31"]] = "#7D7C7C"
  colors_list[["SBS35"]] = "#a3a2a2"
  colors_list[["SBS40"]] = "#4169E1"
  colors_list[["SBS41"]] = "#DA70D6"
  colors_list[["SBS44"]] = "#4B0082"
  colors_list[["SBS88"]] = "#44BD94"
  colors_list[["SBS92"]] = "#BE7CD6"
  colors_list[["SBS93"]] = "#C7A7BC"
  colors_list[["SBS96F"]] = "#FF69B4"
  colors_list[["SBS30"]] = "#483D8B"
  colors_list[["FFPE1"]] = "#6B8E23"
  colors_list[["FFPE2"]] = "#00FFFF"
  colors_list[["Residuals"]] = "#DB7093"
  colors_list[["MSI"]] = "#D92B45"
  

#additional colors if there are signatures that are not in the upper list
additional_colors = c(#"#e377c2",
  #"#ff7f0e",
  #         "#9467bd",
  #         "#bcbd22",
  #         "#8c564b",
  #         "#d62728",
  #         "#2ca02c",
  #         "#17becf",
  #         "#FF4500",
  #         "#8A2BE2",
  #         "#D2691E",
  #         "#006400",
  #         "#FF1493",
  #         "#1E90FF",
  #         "#C71585",
  #         "#FA8072",
  #         "#FF00FF",
  #         "#F4A460",
  #         "#228B22",
  #         "#4169E1",
  #         "#DA70D6",
  #         "#4B0082",
  #         "#8FBC8F",
  #         "#0000FF",
  #         "#DB7093",
  #         "#483D8B",
  #         "#6B8E23",
  #         "#00FFFF",
  #         "#DB7093",
           "#228B22",
           "#0000FF",
           "#F4A460",
           "#bcbd22",
           "#9467bd",
           "#FF4500",
           "#ff7f0e",
           "#e377c2",
           "#663399",
           "#00FF00")

parser <- ArgumentParser(description='Plot signature activities obtained with SigProfiler.')
parser$add_argument('-i', '--infile',
                    dest='infile',
                    required=TRUE,
                    help="The input file with activities.",
                    type="character"
)
parser$add_argument('-o', '--outfile',
                    dest='outfile',
                    required=TRUE,
                    help='The output file.',
                    type="character"
)
parser$add_argument('-g', '--group_sigs',
                    dest='group_sigs',
                    help='A list of signatures to group.',
                    type="character",
                    nargs='+'

)
parser$add_argument('-n', '--group_name',
                    dest='group_name',
                    help='Name of signatures to group.',
                    type="character",
                    default="Group"
)
parser$add_argument('-f', '--ffpe',
                    dest='ffpe',
                    help='The list of ffpe IDs. If given samples will be sorted by FF, FFPE',
                    type="character"
)
parser$add_argument('-e', '--ff',
                    dest='ff',
                    help='The list of ff IDs. If given samples will be sorted by FF, FFPE',
                    type="character"
)
parser$add_argument('-m', '--msi_th',
                    dest='msi_th',
                    help='The MSI threshold, i.e., samples above this threshold (n mutations) will be grouped together on the plot.',
                    type="character"
)
parser$add_argument('-x', '--exclusive',
                    dest='exclusive',
                    help='If different than NULL, the group of all signatures other than the one specified in group_sigs will be created and named according to name given here. E.g., "Real".',
                    type="character"
)
parser$add_argument('-w', '--width',
                    dest='width',
                    default = 0.3,
                    help='the width of one sample bar, in cm.',
                    type="numeric"
)
parser$add_argument('-r', '--residuals',
                    dest='residuals',
                    choices=c("Yes", "No"),
                    default="No",
                    help='Whether to plot residuals, if yes, 96 matrix and signatures have to be provided.',
                    type="character"
)
parser$add_argument('-t', '--matrix',
                    dest='matrix',
                    help='A path to the file with 96 matrix.Rows are contexts, cols are samples, first column contains context names.',
                    type="character"
)
parser$add_argument('-s', '--signatures',
                    dest='signatures',
                    help='A path to the file with signatures. Rows are contexts, cols are signatures, first column contains context names.',
                    type="character"
)
parser$add_argument('-d', '--red_labels',
                    dest='red_labels',
                    help='A list of examples that should be marked with red text in x axis labels.',
                    type="character",
                    nargs='+'
)
parser$add_argument('-p', '--exclude_samples',
                    dest='exclude_samples',
                    help='A list of examples that should be excluded from plot.',
                    type="character",
                    nargs='+'
)
parser$add_argument('-c', '--exclude_signatures',
                    dest='exclude_signatures',
                    help='A list of signatures that should be excluded from plot.',
                    type="character",
                    nargs='+'
)
parser$add_argument('-l', '--group_all',
                    dest='group_all',
                    choices=c("Yes", "No"),
                    default="No",
                    help='Should all signatures be grouped according to predefined groups',
                    type="character"
)
parser$add_argument('-u', '--plot_reconstructed',
                    dest='plot_reconstructed',
                    choices=c("Yes", "No"),
                    default="Yes",
                    help='Should reconstructed profiles be plotted',
                    type="character"
)
parser$add_argument('-a', '--order_alphabetical',
                    dest='order_alphabetical',
                    choices=c("Yes", "No"),
                    default="No",
                    help='Should the samples be ordered alphabetically. If No, they will be ordered by the number of mutations.',
                    type="character"
)

#parser$print_help()

args <- parser$parse_args(
c(
"-i", "/data/cohorts/Results/Signatures/MOMA3-Bladder/INDEL/Fitting_MPOS20/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", 
"-o", "/data/cohorts/Results/Signatures/MOMA3-Bladder/INDEL/Fitting_MPOS20/Residuals_test.png", 
"-r", "Yes",
"-t", "/data/cohorts/Results/Signatures/MOMA3-Bladder/INDEL/Fitting_MPOS20/features_ID83.txt",
"-s", "/data/cohorts/Results/Signatures/MOMA3-Bladder/INDEL/Fitting_MPOS20/custom_signatures.txt"
  )
)

args <- parser$parse_args(
  c(
    "-i", "/data/jurica/CGP/INDELs/MOMA2-Colon-mutect-svaba-PASS-intersect/Fitting_Colon/Assignment_Solution/Activities/Assignment_Solution_Activities.txt", 
    "-o", "/data/jurica/CGP/SNVs/MOMA2-Colon_new/Fitting/Assignment_Solution_Activities_Plot_R_MSIgrouped.png", 
    "-f", "/data/jurica/CGP/ffpe.ids.txt", 
    "-e", "/data/jurica/CGP/ff.ids.txt", 
    "-g", "SBS6","SBS14","SBS15","SBS20","SBS21","SBS26","SBS44", 
    "-n", "MSI", "-m", "50000")
)

args <- parser$parse_args()


activities = read.csv(args$infile, sep = "\t", check.names = F)

if(colnames(activities)[1] != "Samples") { #the file is our custom file, and not coming from sigprofiler
  
  samples = colnames(activities)[-1]
  signatures = activities[,1]
  
  activities = as.data.frame(t(activities[,-1, drop = FALSE]))
  colnames(activities) = signatures
  activities = cbind(Samples = samples, activities)
  
}

nPerSig = colSums(activities[,-1])
sigs = names(nPerSig[nPerSig > 0])

if(!is.null(args$exclude_signatures)) {
  sigs = sigs[!(sigs %in% args$exclude_signatures)]
}

activities = activities[,c("Samples", sigs)]

activities$Samples = str_replace(activities$Samples, "_notshared", "")

#
if(!is.null(args$exclude_samples)) {
  samples_out = c()
  for(s in args$exclude_samples) {
    sample_out = which(grepl(paste0(s,"-"), activities$Samples))
    if(length(sample_out) > 0) {
      samples_out = c(samples_out, sample_out)
    }
    sample_out = which(grepl(paste0("-", s), activities$Samples))
    if(length(sample_out) > 0) {
      samples_out = c(samples_out, sample_out)
    }
  }
  activities = activities[-samples_out,]
}

if(args$residuals == "Yes") {
  cosine_similarities = data.frame(Samples = activities$Samples, cosine_similarity = 0)
  
  #activities = activities[order(-rowSums(activities[,-1])),]
  
  signatures = read.csv(args$signatures, sep = "\t", check.names = F)
  signatures = signatures[, c(colnames(signatures)[1],colnames(activities)[-1])]
  context_matrix = read.csv(args$matrix, sep = "\t", check.names = F)
  colnames(context_matrix)[1] <- "MutationType"
  colnames(signatures)[1] <- "MutationType"
  
  #order to make sure equal ordering of contexts
  #signatures = signatures[order(signatures$MutationType),]
  #context_matrix = context_matrix[order(context_matrix$MutationType),]
  reconstructed = round(as.matrix(signatures[,-1]) %*% as.matrix(t(activities[,-1])), 1)
  colnames(reconstructed) = activities$Samples
  reconstructed = as.data.frame(reconstructed)
  activities$Residuals = 0
  
  context_type = "DBS"
  if(grepl("\\[", signatures$MutationType[1])) {
    context_type = "SBS"
  }
  if(grepl(":", signatures$MutationType[1])) {
    signatures$MutationType = indel_names_COSMIC_to_MutPatterns(signatures$MutationType)
    context_matrix$MutationType = indel_names_COSMIC_to_MutPatterns(context_matrix$MutationType)
    signatures = signatures[match(INDEL_CATEGORIES, signatures$MutationType),]
    context_matrix = context_matrix[match(INDEL_CATEGORIES, context_matrix$MutationType),]
    context_type = "ID"
  }
  
  plots = list()
  
  for(i in 1:nrow(activities)) {
    total_mut = sum(activities[i,-1]) #sums also residuals, but at this point they are 0
    orig_row = context_matrix[,as.character(activities$Samples[i])]
    recon_row = reconstructed[,as.character(activities$Samples[i])]
    cosine_similarities$cosine_similarity[i] = cos_sim(orig_row, recon_row)
    
    if(args$plot_reconstructed == "Yes") {
      mat_96 = data.frame(original = orig_row, reconstructed = recon_row)
      row.names(mat_96) = signatures$MutationType
      
      if(context_type == "SBS") {
        p = plot_96_profile(mat_96, condensed = T, relative = F)  
      }
      if(context_type == "DBS") {
        row.names(mat_96) = str_replace(row.names(mat_96), pattern = ">", "_") 
        p = MutationalPatterns::plot_dbs_contexts(mat_96, same_y = F, condensed = T)
      }
      if(context_type == "ID") {
        p = MutationalPatterns::plot_indel_contexts(mat_96, condensed = T) 
        p = p + theme(legend.key.size = unit(0.3, "cm"), legend.position = "top", legend.title=element_blank()) + guides(fill=guide_legend(nrow=3,byrow=TRUE)) 
      }
      p = p + ggtitle(paste0(activities$Samples[i], ", cosine similarity: ", round(cosine_similarities$cosine_similarity[i], 3))) + 
        theme(strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
        theme(text=element_text(size=7)) + ylab("Number of mutations") + xlab("")
      
      plots[[i]] = p
    }
    
    res = orig_row - recon_row
    res = sum(abs(res))
    #rescale activities so activities + residuals = total mutations
    activities[i,-1] = activities[i,-1] * (total_mut - res) / total_mut
    activities$Residuals[i] = res
  }
  
  if(args$plot_reconstructed == "Yes") {
    p_all = plot_grid(plotlist = plots, ncol = 1)
    ggsave(paste0(dirname(args$outfile), "/original_vs_reconstructed_profiles.pdf"), p_all, width = 20, height = 9*length(plots), units = "cm", limitsize = F)
  }
}

if(!is.null(args$ff) & !is.null(args$ffpe)) {
  ffpe = read.csv(args$ffpe, header = F)
  ff = read.csv(args$ff, header = F)
  
  activities$Samples[activities$Samples %in% ffpe$V1] = paste0(activities$Samples[activities$Samples %in% ffpe$V1], "_FFPE")
  activities$Samples[activities$Samples %in% ff$V1] = paste0(activities$Samples[activities$Samples %in% ff$V1], "_FF")
  
  activities.long = reshape2::melt(activities, id.vars=c("Samples"))
  
  ##FF and FFPE grouped and ordered by n muts
  cohort_order = rowSums(activities[,-1])
  names(cohort_order) = activities$Samples
  order_cohort_FF = cohort_order[!grepl("FFPE", names(cohort_order))]
  order_cohort_FF = order_cohort_FF[order(-order_cohort_FF)]
  order_cohort_FFPE = cohort_order[grepl("FFPE", names(cohort_order))]
  order_cohort_FFPE = order_cohort_FFPE[order(-order_cohort_FFPE)]
  
  order = c(names(order_cohort_FF), names(order_cohort_FFPE))
  
  activities.long$Samples = factor(activities.long$Samples, levels = order)
  
} else {
  cohort_order = rowSums(activities[,-1])
  names(cohort_order) = activities$Samples
  if(args$order_alphabetical == "Yes") {
    cohort_order = cohort_order[order(names(cohort_order))]
  } else {
    cohort_order = cohort_order[order(-cohort_order)]
  }
  
  activities.long = reshape2::melt(activities, id.vars=c("Samples"))
  
  activities.long$Samples = factor(activities.long$Samples, levels = names(cohort_order))
}

if(!is.null(args$msi_th)) {
  #HACK for MOMA-Colon-All INDEL
  #act = read.csv("/data/cohorts/Results/Signatures/MOMA-Colon-All/SNV/Fitting_afterNS_donorsExcluded_SBS36_SBS3_SBS5_APOBEC/activities_fitting.txt", sep = "\t", check.names = F)
  #cnames = colnames(act)
  #cnames[cnames %in% ff$V1] = paste0(cnames[cnames %in% ff$V1], "_FF")
  #cnames[cnames %in% ffpe$V1] = paste0(cnames[cnames %in% ffpe$V1], "_FFPE")
  #colnames(act) = cnames
  #act.n = colSums(act[,-1])
  #MSI = names(act.n[act.n > 50000])
  nPerSample = aggregate(value ~ Samples , data = activities.long, sum)
  MSI = nPerSample$Samples[nPerSample$value > as.numeric(args$msi_th)]
}

if(args$group_all == "Yes") {
  for(group in names(signature_groups)) {
    activities.long =  groupSigs(activities.long, signature_groups[[group]], group)
  }
}

if(!is.null(args$group_sigs)) {
  activities.long =  groupSigs(activities.long, args$group_sigs, args$group_name)
}



if(!is.null(args$exclusive)) {
  #colors = c(colors_real_artefact, colors)
  activities.long =  groupSigs(activities.long, sigs[!(sigs %in% args$group_sigs)], args$exclusive)
  colors = colors_real_artefact
} else {
  sigs = unique(activities.long$variable)
  #take subset of relevant colors
  colors = rep("", length(sigs))
  names(colors) = sigs
  
  colors[names(colors) %in% names(mut.sig.colors)] = mut.sig.colors[names(colors)[names(colors) %in% names(mut.sig.colors)]]
  
  #colors_temp = unlist(colors_list[which(names(colors_list) %in% sigs)])
  #colors[names(colors) %in% ] = colors_temp
  
  #if(length(colors_temp) < length(sigs)) { #if some signatures are not in the list
  if(sum(names(colors) %in% names(mut.sig.colors)) < length(sigs)) {
    #if(length(colors_temp) == 0) {
    if(sum(names(colors) %in% names(mut.sig.colors)) == 0) {
      colors = unlist(colors_list)[1:length(sigs)]
      names(colors) = sigs
    } else {
      colors[colors == ""] = additional_colors[1:sum(colors == "")]
    } 
  }
}
 
if(!is.null(args$msi_th)) {
  activities.long$status = "MSS"
  activities.long$status[activities.long$Samples %in% MSI] = "MSI"
}

#if(!is.null(args$msi_th)) {
#  activities.long$status = paste0(activities.long$status, "_", activities.long$FF)
#  activities.long.rel$status = paste0(activities.long.rel$status, "_", activities.long.rel$FF)
#}


if(!is.null(args$ff) & !is.null(args$ffpe)) {
  activities.long$FF = "FF"
  activities.long$FF[grepl("_FFPE", activities.long$Samples)] = "FFPE"
}

p = ggbarplot(activities.long, x = "Samples", y = "value", fill = "variable", palette = colors) + 
  theme(strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  theme(text=element_text(size=7)) + 
  ylab("Numer of mutations") +
  rotate_x_text() + xlab("")+ 
  theme(legend.key.size = unit(0.5, "cm"), legend.margin=margin(t = -0.1, b = -0.1, unit='cm')) +
  labs(fill = "Signatures")

##
##
if(FALSE) {
  p.snv = ggbarplot(activities.long.svn, x = "Samples", y = "value", fill = "variable", palette = colors) + 
    theme(strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
    theme(text=element_text(size=7)) + 
    ylab("Number of mutations") +
    rotate_x_text() + xlab("")+ 
    theme(legend.key.size = unit(0.5, "cm"), legend.margin=margin(t = -0.1, b = -0.1, unit='cm')) +
    labs(fill = "Signatures")
  
  p.rel.snv = ggbarplot(activities.long.rel.snv, x = "Samples", y = "value", fill = "variable", palette = colors) + 
    theme(strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
    theme(text=element_text(size=7)) + 
    ylab("Relative contribution") +
    rotate_x_text() + xlab("") + 
    theme(legend.key.size = unit(0.5, "cm"), legend.margin=margin(t = -0.1, b = -0.1, unit='cm')) +
    labs(fill = "Signatures")
  
  p.snv.indel = plot_grid(p.snv, p, nrow = 2)
  p.rel.snv.indel = plot_grid(p.rel.snv, p.rel, nrow = 2)
  
  p.both = plot_grid(p.snv.indel, p.rel.snv.indel, nrow = 2)
  
  
  width = nrow(activities) * 0.28
  #height = 30
  height = 17
  
  args$outfile = "/data/jurica/CGP/SNVs/MOMA2-Colon_new/MOMA2-Colon_fitting_activities_snv_indel.png"
  
  if(grepl(".pdf", args$outfile)) {
    ggsave(args$outfile, p.both, width = width, height = height, units = "cm", limitsize = F, device = "pdf")
  } else {
    ggsave(args$outfile, p.both, width = width, height = height, units = "cm", limitsize = F, dpi = 300, device = "png")
  }
}
##
##

activities.long.rel = relativeContributions(activities.long)

p.rel = ggbarplot(activities.long.rel, x = "Samples", y = "value", fill = "variable", palette = colors) + 
  theme(strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  theme(text=element_text(size=7)) + 
  ylab("Relative contribution") +
  rotate_x_text() + xlab("") + 
  theme(legend.key.size = unit(0.5, "cm"), legend.margin=margin(t = -0.1, b = -0.1, unit='cm')) +
  labs(fill = "Signatures")

if(!is.null(args$msi_th)) {
  p = p + facet_grid(. ~ status, scales = "free", space = "free_x")
  p.rel = p.rel + facet_grid(. ~ status, scales = "free", space = "free_x")
}

#mark samples in x axis labels with red
if(!is.null(args$red_labels)) {
  #to color x axis labels 
  labels = ggplot_build(p.rel)$layout$panel_params[[1]]$x$breaks
  
  excluded = args$red_labels #c("6062","5471","5148","4855","5683","5701")
  colours = rep("black", times = length(labels))
  for(e in excluded) {
    colours[grepl(e, labels)] = "red"
  }
  
  p = p + theme(axis.text.x = element_text(colour = colours))
  p.rel = p.rel + theme(axis.text.x = element_text(colour = colours))
}

if(FALSE) {
  p.FF = p
  p.rel.FF = p.rel
  
  p.FFPE = p
  p.rel.FFPE = p.rel
  
  p.FF = p.FF + ggtitle("FF")
  p.rel.FF = p.rel.FF + ggtitle("FF")
  p.FFPE = p.FFPE + ggtitle("FFPE")
  p.rel.FFPE = p.rel.FFPE + ggtitle("FFPE")
  
  nFF = nrow(ff)
  nFFPE = nrow(ffpe)
  
  p.both.FF = plot_grid(p.FF, p.FFPE, ncol = 2, rel_widths = c(nFF / (nFF+nFFPE), nFFPE / (nFF+nFFPE)))
}  

#p = p + coord_flip()
#p.rel = p.rel + coord_flip()

if(args$residuals == "Yes") {
  cosine_similarities$Samples = factor(cosine_similarities$Samples, levels = levels(activities.long$Samples))
  p.sim = ggbarplot(cosine_similarities, x = "Samples", y = "cosine_similarity", fill = colors[ncol(activities)-1]) + 
    theme(strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
    theme(text=element_text(size=7)) + 
    ylab("Cosine similarity (original vs. reconstructed)") +
    rotate_x_text() + xlab("")
  
  p.both = plot_grid(p, p.rel, p.sim, nrow = 3)
  height = 36
} else {
  p.both = plot_grid(p, p.rel, nrow = 2)
  height = 24
}

width = nrow(activities) * args$width

if(width < 10) {
  width = 10
}

#width = 18
#height = 10

if(grepl(".png", args$outfile)) {
  ggsave(args$outfile, p.both, width = width, height = height, units = "cm", limitsize = F, dpi = 300, device = "png")
  ggsave(str_replace(args$outfile, ".png", ".pdf"), p.both, width = width, height = height, units = "cm", limitsize = F, device = "pdf")
} 
if(grepl(".pdf", args$outfile)) {
  ggsave(str_replace(args$outfile, ".pdf", ".png"), p.both, width = width, height = height, units = "cm", limitsize = F, dpi = 300, device = "png")
  ggsave(args$outfile, p.both, width = width, height = height, units = "cm", limitsize = F, device = "pdf")
}
