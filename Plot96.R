library(optparse)
library(data.table)
library(ggplot2)
library(cowplot)
library(stringr)
library(grid)
library(gridExtra)


option_list = list(
  make_option(c("-t", "--threads"), type="numeric", default=8,
              help="threads to use", metavar="numeric"),
  make_option(c("-i", "--input"), type="character", default="/data/NoiseSupression/MASTERTABLE.csv",
              help="Table with input data. Currently expect the data to have a column named NC_3, containing the trinucleotide context 
                for each SNV.", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Either output directory (ending with a /), or a basename for plot(s), or an (almost) full path consisting of the output 
                dir and basename. Defaults to current working directory.", metavar="character"),
  make_option(c("-f", "--filter"), type="character", default=NULL,
              help="Condition for initial filtering - decides what data is included in the analysis. Passed as string to data.table, but OMIT
                BACKTICKS when writing on command line.", metavar="character"),
  make_option(c("-A", "--splitrows"), type="character", default=NULL,
              help="First level of splitting. Determines how an individual plot is split to rows. Can be a name of column containing a 
                categorical variable (factor), or a codition like the one in 'filter'.", metavar="character"),
  make_option(c("-B", "--splitplotsrow"), type="character", default=NULL,
              help="Second level of splitting. Determines how data is split to individual plots, represented in rows in the resulting pdf. 
                Can be a name of column containing a categorical variable (factor), or a codition like the one in 'filter'.", metavar="character"),
  make_option(c("-C", "--splitplotscol"), type="character", default=NULL,
              help="Third level of splitting. Determines how data is split to individual plots, represented in columns in the resulting pdf. 
                Can be a name of column containing a categorical variable (factor), or a codition like the one in 'filter'. If you don't specify 
                anything plots will be laid out in 2 columns.", metavar="character"),
  make_option(c("-D", "--splitfiles"), type="character", default=NULL,
              help="Fourth level of splitting. Determines how plots are organised into individual pdfs. Can be a name of column containing a 
                categorical variable (factor), or a codition like the one in 'filter'.", metavar="character"),
  make_option(c("-s", "--sort"), type="character", default="no",
              help="Options: 'byname'(alphabetically), 'bycount'(plots with highest numbers of SNVs first), 'byperc'(look at percentages which are 
                not 100, and sort plots descendingly), 'no'(leave order as in input data).", metavar="character"),
  make_option(c("-p", "--showperc"), type="character", default="yes",
              help="Percentage of SNVs in each horizontal facet (plot row) is calculated in relation to facet with highest number, and shown on the plot 
                by default. If you set this to 'no', percentage won't be shown. Useful when e.g. comparing entire cohort with public data such as PCAWG.",
                metavar="character"),
  make_option(c("-z", "--fontsize"), type="numeric", default=NULL,
              help="Font size for plot row titles (what it says in grey box to the right).", metavar="numeric"),
  make_option(c( "--figureformat"), type="character", default="pdf",
              help="Format to save in. Default is pdf.", metavar="character"),
  make_option(c("-e", "--dropempty"), type="character", default="no",
              help="Whether to drop levels when they're empty ('yes') or show even the empty ones ('no'; default). Note that this will not influence levels
                inside one facet of 96 plot, i.e. all contexts will be kept no matter whether empty or not. Use this option when plotting e.g. detections
                from all plasmas in multiple samples.", metavar="character"),
  make_option(c("-x", "--hideall"), type="character", default="no",
              help="By default, scripts adds 'all SNVs' as a row in the plot. If this is set to 'yes', that part is skipped. Primarily intended for use 
                with custom pre-prepared tables (e.g. already melted ones, containing a column with sequential filters), or when plotting e.g. tumor SNVs
                vs detections in plasma.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt_parser@description <- "Calculates and plots 96 matrices per sample, split by specified condition. Tip: to get 'cohort:subjectID.sampleID' as a splitting level without passing it through paste(), use shorthand 'fullID'. 'longID' for 'subjectID.sampleID'."
opt_parser@usage <- "Rscript /data/scripts/NoiseSupression/plots/Plot96.R -i /path/to/input -o /path/to/outputname -f 'column1 == TRUE' -A 'column2 == value2' -B 'column3+column4 > value3' -C 2 -D 'column5' -s 'no' -x 'no'"
opts = parse_args(opt_parser);
input <- opts$input
output <- opts$output
initialfilter <- opts$filter
rowsplit <- opts$splitrows
plotsplitrow <- opts$splitplotsrow
plotsplitcol <- opts$splitplotscol
filesplit <- opts$splitfiles
orderplots <- opts$sort
showperc <- opts$showperc
labfontsize <- opts$fontsize
dropempty <- opts$dropempty
dontshowall <- opts$hideall
corenum <- opts$threads
setDTthreads(corenum)



####################################
### get data and evaluate conditions
####################################

### parse conditions

# extract column names and prepare command-line input for passing to data.table:
cond_raw <- list(initialfilter=initialfilter, 
                 rowsplit=rowsplit, 
                 plotsplitrow=plotsplitrow, plotsplitcol=plotsplitcol, 
                 filesplit=filesplit)                 # all conditions
if(grepl("\\.tsv$", input)) {
    allcolnames <- unlist(str_split(readLines(con=input, n=1), "\t"))
} else {
    allcolnames <- unlist(str_split(readLines(con=input, n=1), ","))
}

parseCondition <- function(condition, allcolnames) {
  
  if(is.null(condition)) {
    list(NA, NA)
  } else {
    # locate all matches:
    locations <- as.data.table(do.call(rbind, str_locate_all(condition, allcolnames)))
    # keep the longest, and in case of a tie, the one which appears first:
    # (reason being that sometimes a column name can be included in another columns' name, or be a part of the contents of the column):
    locations <- locations[!is.na(start), ][, len := end-start][order(-len, start)]
    locations <- locations[, multmatches := ifelse((duplicated(start) | duplicated(end)), "yes", "no")][multmatches=="no", ]
    # extract required column name from match and add it to the string of those which are needed:
    cols <- str_sub(condition, start = locations$start, end = locations$end)
    
    # add backticks to appropriate places (to enable parsing by data.table later):
    mayneedticks <- locations[which(str_detect(cols, "[^a-zA-Z0-9]")), ]
    # additional check: we don't want to put ticks around something which is part of a string already (e.g. 'risky' character in column 
    # contents referenced in condition). This checks if there are single quotation marks immediately before or after the phrase.
    # Note that this will fail if the risky phrase is the same as one of the column names and in the middle of the string.
    checkforticks <- str_sub(condition, start = mayneedticks$start-1, end = mayneedticks$end+1)
    definiteticks <- mayneedticks[!which(str_detect(checkforticks, "'"))]
    if(length(definiteticks > 0)) {
      tickpositions <- sort(c(definiteticks$start, definiteticks$end+1))
      for(i in 1:length(tickpositions)) {
        condition <- stringi::stri_sub_replace(condition, from=tickpositions[i], to=tickpositions[i]-1, value = "`")
        tickpositions <- tickpositions+1
      }
    }
    # return corrected condition and needed column(s):
    list(condition, cols)
  }
}
cond_corrected <- lapply(cond_raw, parseCondition, allcolnames)
cond <- lapply(cond_corrected, function(x) x[[1]])
requiredcols <- unique(unlist(lapply(cond_corrected, function(x) x[[2]])), use.names = FALSE)
requiredcols <- requiredcols[!is.na(requiredcols)]

# special option to include full sample info (cohort, subject, sample) but avoid an ugly title*:
# (which we would get if we passed e.g. str_c() to command line arguments)
if("longID" %in% unlist(cond)) {
  requiredcols <- c(requiredcols, "subjectID", "sampleID")
} else if("fullID" %in% unlist(cond)) {
  requiredcols <- c(requiredcols, "cohort", "subjectID", "sampleID")
}


### evaluate conditions

# read in data:
#selectcols <- unique(c("CHROM", "POS", "REF", "ALT", "NC_3", requiredcols))
selectcols <- unique(c("REF", "ALT", "NC_3", requiredcols))

dt <- fread(input, select = selectcols)

# create long or full ID columns if requested, and adjust cond and cond_corrected accordingly:
if("longID" %in% unlist(cond)) {
  dt[, longID := str_c(subjectID, ".", sampleID)]
  whichcond <- which(unlist(cond)== "longID")
  cond[[whichcond]] <- "longID"
  cond_corrected[[whichcond]][[2]] <- "longID"
}
if("fullID" %in% unlist(cond)) {
  dt[, fullID := str_c(cohort, ":", subjectID, ".", sampleID)]
  whichcond <- which(unlist(cond)== "fullID")
  cond[[whichcond]] <- "fullID"
  cond_corrected[[whichcond]][[2]] <- "fullID"
}

# filter and split data:
eval(parse(text=str_c("dt <- dt[", cond$initialfilter, ", ]")))
for(i in 2:length(cond)) {
  if(is.na(cond[[i]])) {
    eval(parse(text=str_c("dt[, ", names(cond)[i], " := 'all SNVs']")))
  } else {
    eval(parse(text=str_c("dt[, ", names(cond)[i], " := ", cond[[i]], "]")))
  }
}

# initial filtering is done so removing that from cond, for easier handling in rest of the script:
cond$initialfilter <- NULL

# keep NAs created by conditions as special category:
# (this is not the same as drop empty! we always want to see if there's NA for a filter or condition)
dt[, (names(cond)) := lapply(.SD, as.character), .SDcols = names(cond)]
dt[, (names(cond)) := lapply(.SD, function(x) ifelse(is.na(x), "information not available", x)), .SDcols = names(cond)]




##############################
### make nice labels and names
##############################

### extract context and substitution type

# extract substitution type and context from trinucleotide context:
dt[, context_ref := ifelse(REF %in% c("C", "T"), REF, chartr(old="GA", new="CT", REF))
   ][, context_alt := ifelse(REF %in% c("C", "T"), ALT, chartr(old="ACTG", new="TGAC", ALT))
     ][, substitution := str_c(context_ref, ">", context_alt)
       ][, context := str_c(str_sub(NC_3, 1, 1), ".", str_sub(NC_3, 3, 3))]


### get counts and percentages

tb <- dt[, .N, by=c(names(cond), "substitution", "context")]   # N in individual bar
tb <- tb[!str_detect(context, "[^ATCG\\.]")]     # remove anything that has a non-ATCG character (rarely happens)

if(dontshowall == "no") {   # adds 'all SNVs' category
  
  # add "all SNVs" category for plotting:
  tb[, N_96 := sum(N), by=c(names(cond)[-1], "substitution", "context")]
  tb_all <- unique(copy(tb)[, c("rowsplit", "N") := NULL])
  setnames(tb_all, "N_96", "N")
  tb_all[, rowsplit := "all SNVs"]
  tb[, N_96 := NULL]
  tb <- rbind(tb_all, tb, use.names=TRUE, fill=FALSE)
}

tb[, N_row := sum(N), by=c(names(cond))]   # how many SNVs in horizontal facet of plot
# percentages for categories will be expressed in relation to N_base = N_row with highest number of SNVs:
tb[, N_base := max(N_row), by=c(names(cond)[-1])]
tb[, perc_row := round(N_row / N_base *100, 1)]   


### nicer labels for TRUE/FALSE categories

# rename levels if they are currenty named TRUE/FALSE:
colstochange <- lapply(unique(tb[, names(cond), with=FALSE]), function(x) all(x %in% c("TRUE", "FALSE", "information not available")))
colstochange <- names(cond)[unlist(colstochange)]
tb[, (names(cond)) := lapply(.SD, function(x) ifelse(x=="TRUE", "After filtering", 
                                              ifelse(x=="FALSE", "Filtered out", x))), .SDcols = names(cond)]


### create output file name

# part of file name reflecting filtering:
partoffilename1 <- unique(unlist(lapply(cond_corrected[1], function(x) x[[2]])), use.names = FALSE)
partoffilename1 <- partoffilename1[!is.na(partoffilename1)]
partoffilename1 <- paste0("_FILTER_", paste0(partoffilename1, collapse="_"))
partoffilename1 <- str_remove_all(partoffilename1, "[^a-zA-Z0-9_\\-]")   # remove risky characters
if(partoffilename1=="_FILTER_") partoffilename1 <- ""    # remove if empty

# part of file name reflecting faceting:
partoffilename2 <- unique(unlist(lapply(cond_corrected[-c(1,5)], function(x) x[[2]])), use.names = FALSE)
partoffilename2 <- partoffilename2[!is.na(partoffilename2)]
partoffilename2 <- paste0("_BY_", paste0(partoffilename2, collapse="_"))
partoffilename2 <- str_remove_all(partoffilename2, "[^a-zA-Z0-9_\\-]")
if(partoffilename2=="_BY_") partoffilename2 <- ""

partoffilename <- paste0(partoffilename1, partoffilename2)

# function which sanitizes conditions so they're safe to use as filenames:
sanitizeFilename <- function(x) {
  x <- str_replace_all(x, pattern = "<={0,1}", replacement = "lessthan")
  x <- str_replace_all(x, pattern = ">={0,1}", replacement = "morethan")
  x <- str_replace_all(x, pattern = "==", replacement = "equals")
  x <- str_replace_all(x, pattern = "&", replacement = "AND")
  x <- str_replace_all(x, pattern = "\\|", replacement = "OR")
  x <- str_remove_all(x, "[^a-zA-Z0-9_\\-]")
  return(x)
}




#######################
### levels and ordering
#######################

### plot level

# ensure all plots show all 96 substitutions, whether they're empty or not:
contextlevels <- paste0(rep(c("A", "C", "G", "T"), each=4), ".", c("A", "C", "G", "T"))
subslevels <- str_c(rep(c("C", "T"), each=3), ">", c("A", "G", "T", "A", "C", "G"))
tb[, context := factor(context, levels = contextlevels)
   ][, substitution := factor(substitution, levels = subslevels)]

# correctly order horizontal facets ('rows') in a plot:
# (order of rows is kept the same as in original dataset, except for 'allSNVs' (which will come first) and a few common applications where it makes sense)
# (all plots will have either the same order of rows, or (in case of dropempty==TRUE) be ordered by the same general principle)
rowsplit_categories <- unique(tb$rowsplit)
rowsplit_order <- unique(c("all SNVs", "Filtered out", "After filtering", "Shared", "Not shared", rowsplit_categories))
rowsplit_order_final <- rowsplit_order[rowsplit_order %in% rowsplit_categories]
tb[, rowsplit := factor(rowsplit, levels = unique(rowsplit_order_final))]


### pdf level

# get plots into order in which they will appear in pdf(s):
# (splitting into individual pdfs later in the script will not influence order established here)
tbh <- unique(copy(tb)[, c("substitution", "context") := NULL])       # helper table

if(orderplots == "byname") {
  plotsplitrow_order <- sort(unique(tbh$plotsplitrow))
} else if(orderplots == "bycount") {
  plotsplitrow_order <- tbh[order(N_base, decreasing = TRUE)][, .SD[1], by=c(names(cond)[-1])]$plotsplitrow
} else if(orderplots == "byperc") {
  plotsplitrow_order <- tbh[perc_row != 100][order(N_row, decreasing = TRUE)][, .SD[1], by=c(names(cond)[-1])]$plotsplitrow
} else {   # else leave as in input
  plotsplitrow_order <- unique(tbh$plotsplitrow)
}
tb[, plotsplitrow := factor(plotsplitrow, levels = unique(plotsplitrow_order))]

# plotsplitcol has to be same in all pdfs - sorting in the same way as rowsplit
plotsplitcol_categories <- unique(tb$plotsplitcol)
plotsplitcol_order <- unique(c("all SNVs", "Filtered out", "After filtering", "Shared", "Not shared", plotsplitcol_categories))
plotsplitcol_order_final <- plotsplitcol_order[plotsplitcol_order %in% plotsplitcol_categories]
tb[, plotsplitcol := factor(plotsplitcol, levels = unique(plotsplitcol_order_final))]

tb <- tb[order(plotsplitcol, plotsplitrow, rowsplit)]




####################
### plot single plot
####################

# define some universal plotting parameters:
plotcolors <- c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
names(plotcolors) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

# function to produce an individual plot:
singlePlot <- function(x, showperc=showperc, dropempty=dropempty) {
  
  # make title (which shows row, column, and file splitting, if they exist):
  title_filesplit <- ifelse(is.null(filesplit), "", str_c(filesplit, ": "))
  title_splitcol <- ifelse(is.null(plotsplitcol), "", str_c(plotsplitcol, ": "))
  title_splitrow <- ifelse(is.null(plotsplitrow), "", str_c(plotsplitrow, ": "))
  
  plottitle <- c(paste0(title_filesplit, x$filesplit[1], "\n"),
                 paste0(title_splitcol, x$plotsplitcol[1], "\n"),
                 paste0(title_splitrow, x$plotsplitrow[1], "\n"))
  plottitle <- paste0(str_remove(plottitle, "all SNVs\n"), collapse = "")
  
  # make subtitle - show filters (initial and current):
  title_initialfilter <- ifelse(is.null(initialfilter), "none", initialfilter)
  title_rowsplit <- ifelse(is.na(cond$rowsplit), "none", cond$rowsplit)
  plotsubtitle <- str_c("Dataset filtered before analysis: ", title_initialfilter,
                        "\nCurrent filter: ", title_rowsplit)
  
  # on demand drop empty horizontal facets i.e. plot rows (achieved by re-leveling rowsplit variable):
  if(dropempty=="yes") {
    x[, rowsplit := factor(as.character(rowsplit), levels = unique(as.character(x$rowsplit)))]
  }
  nrrowsplot <- length(unique(x$rowsplit))
  
  # create labels and set their size:
  if(showperc=="yes") {
    labs <- unique(str_c(x$rowsplit, "\n", "N=", x$N_row, " (", x$perc_row, "%)"))
  } else {
    labs <- unique(str_c(x$rowsplit, "\n", "N=", x$N_row))
  }
  names(labs) <- unique(str_c(x$rowsplit))
  maxlabchar <- max(nchar(labs))
  if(!is.null(labfontsize)) {
    fontsize_y <- labfontsize
  } else {
    fontsize_y <- ifelse(maxlabchar > 40, 7, ifelse(maxlabchar > 30, 8, ifelse(maxlabchar > 20, 10, 12)))
  }
  
  # plot:
  ggplot(x, aes(x = context, y = N, fill = substitution, width = 1)) + 
    geom_bar(stat = "identity", position = "identity", colour = "black", linewidth = 0.2) + 
    scale_fill_manual(values = plotcolors) + 
    facet_grid(rowsplit ~ substitution, 
               scales = "free_y", 
               drop = FALSE,
               labeller = labeller(.rows=labs)) + 
    labs(y = "Number of mutations") + 
    guides(fill = "none") +
    theme_bw() +
    scale_x_discrete(drop=FALSE) +
    theme(plot.title = ggtext::element_textbox_simple(),
          axis.title.y = element_text(size = 12, vjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = fontsize_y),
          panel.grid.major.x = element_blank(),
          panel.spacing.x = unit(0, "lines"), 
          aspect.ratio = 1) +
    ggtitle(plottitle, plotsubtitle)
}




###############
### tidy pdf(s)
###############

# all pdfs will have equal number of columns (which is 2 if not specified otherwise):
nrcolspdf <- length(unique(tb$plotsplitcol))
if(is.na(cond$plotsplitcol) & !is.na(cond$plotsplitrow) & dropempty=="no") nrcolspdf <- 2
if(is.na(cond$plotsplitcol) & !is.na(cond$plotsplitrow) & dropempty=="yes") nrcolspdf <- 1   # just one column if we're dropping empty (otherwise we get ton of blank spaces)
onlyoneplot <- length(unique(tb$plotsplitcol))==1 & length(unique(tb$plotsplitrow))==1
if(isTRUE(onlyoneplot)) nrcolspdf <- 1   # or in special case where there's only 1 plot, there's of course only 1 pdf column

# make and export all necessary plots:
bypdfdata <- split(tb, by="filesplit")
lapply(bypdfdata, function(onepdfdata) {

  # some per-pdf dimensions:
  nrrowsperplot <- onepdfdata[, uniqueN(rowsplit), by=plotsplitrow]$V1   # numbers of horizontal facets in all plots in this pdf
  nrrowspdf <- max(onepdfdata[, uniqueN(plotsplitrow), by=plotsplitcol]$V1, na.rm = TRUE)   # we want all columns in pdf to have the same number of rows
  
  # all facets are 4 cm in height, plus additional 2 cm space per plot for margins and titles:
  # NOTE: if there are columns defined, we take the space which biggest plot occupies as 'basic' unit when calculating how big the pdf will be in the end.
  # This has no visible effect when all levels are kept (dropempty='no'), but if we're dropping empty levels it will cause some big margins and blank
  # spaces. But this is necessary if we want to keep columns parallel.
  
  # However, if the columns are not defined, default is col=2 when all plots are the same size (keeping all levels), BUT single column if empty levels are
  # dropped. That way we can easily calculate minimal plot length for the second case and avoid gaps.
  # (in the future, it would be good to introduce minimal height per row - that way we could column-ize differently-sized plots without getting too much 
  # blanks - possible usage in conjunction with rel_heights in plot_grid)
  if(nrcolspdf > 1) {
    heightfactor <- max(nrrowsperplot, na.rm = TRUE) * 4 + 2
    pdfheight <- nrrowspdf * heightfactor
    rel_heights <- rep(1, nrrowspdf)
  } else {
    pdfheight <- sum(nrrowsperplot, na.rm = TRUE) * 4 + length(nrrowsperplot) * 3
    rel_heights <- (nrrowsperplot*4) + 2
  }

  # create plots:
  byplotdata <- split(onepdfdata, by=c("plotsplitcol", "plotsplitrow"))
  onepdf <- lapply(byplotdata, function(x) singlePlot(x, showperc=showperc, dropempty=dropempty))
  
  # create layout (add blank spaces to ensure proper "faceting" in pdfs):
  if((nrrowspdf==1)|(nrcolspdf==1)) {   # no spaces needed
    onepdf_withspaces <- onepdf
  } else {
    spaces <- onepdfdata[, uniqueN(plotsplitrow), by=plotsplitcol][order(plotsplitcol)]$V1
    blankspaces <- nrrowspdf - spaces
    plotspaces <- c(0, cumsum(spaces))
    onepdf_withspaces <- list()
    for(i in 1:(length(plotspaces)-1)) {
      onepdf_withspaces <- c(onepdf_withspaces, 
                             onepdf[(plotspaces[i]+1):plotspaces[i+1]],
                             replicate(blankspaces[i], NULL))
    }
  }
  
  # create filename:
  if(is.null(output)) output <- str_c(getwd(), "/")
  #filesplit_full <- str_c(filesplit, unique(onepdfdata$filesplit), plotsplitcol, plotsplitrow, rowsplit, initialfilter)   # longer filename option - for Paz
  filesplit_full <- str_c(filesplit, unique(onepdfdata$filesplit))                                                         # file-specific part of filename
  sanitize_filename <- sanitizeFilename(filesplit_full)
  
  if(str_ends(output, "/")) {
    filename <- paste0("96plot_", sanitize_filename, partoffilename)
    full_filename <- str_c(output, filename, ".",opts$figureformat)
  } else if(str_ends(output, str_c(".",opts$figureformat))) {   # no need to sanitize if you provided full filename for output
    full_filename <- output
  } else {
    full_filename <- str_c(output, "_96matrix", sanitize_filename, ".",opts$figureformat)
  }

  # adjust appearance of output if there is splitting to multiple plots in pdf but no particular division for columns:
  # (if division to columns is not defined then 2 cols is default UNLESS we're dripping empty levels - in that case, anything except 1 column results 
  # in ton of empty spaces if plots vary greatly in size)
  onlysinglerows <- all(nrrowsperplot==1)
  if(is.null(plotsplitcol) & !is.null(plotsplitrow) & dropempty=="no" & onlyoneplot==FALSE) {
    pdfheight <- ceiling(pdfheight/2)  
    if(onlysinglerows) pdfheight <- pdfheight + 6    # for pdfs where all plots have a single row, ratio of space needed for plot title and actual plot is such that we end up with lack of space when there is a small number of plots; this will fix the problem and won't influence much when there is a lot of plots
    plot_grid(plotlist = onepdf_withspaces, ncol = nrcolspdf, byrow = TRUE, rel_heights = rel_heights, align = "h", axis="t")
  } else {
    plot_grid(plotlist = onepdf_withspaces, ncol = nrcolspdf, byrow = FALSE, rel_heights = rel_heights, align = "h", axis="t")
  }
  
  # export:
  ggsave(filename = full_filename, width = nrcolspdf * 22, height = pdfheight, units = "cm", device = opts$figureformat, limitsize = FALSE)
})

