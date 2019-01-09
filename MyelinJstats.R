# Statistical analysis and graph production for percentage (%) neurite density and % 
# myelination from MyelinJ image analysis. Each experimental condition has its own
# folder, within which are comma separated values (csv) files for each experiment. Each
# cvs file contains % neurite density and % myelination for each image.
# %%%%% statistics %%%%%%
# ggpubr is used for statistical analysis where the mean % neurite density and % myelination
# is calculated for each experimental repeat and compared using unpaired T test. If multiple
# comparisons are performed then the false discovery rate (FDR) is used to adjust for
# multiple comparisons. This will compare all experimental conditionditions to each other. The user
# also has the option of comparing all experimental conditionsitions to control only by putting
# (control) next to the name of the control (using the GUI in ImageJ). csv files are saved
# denoting all comparisons and exact p values. 
# %%%%% graphs %%%%%%
# ggpubr is used for producing graphs. Two graphs are produced for % neurite density and
# % myelination. The first graph illustrates the variation between images for each experiment.
# The second is in the style of a graph for publication where the repeats for each 
# experimental conditionsitions are conpared as a bar graph, with the individual points also
# plotted. Any statistical comparison with p < 0.05 will be displayed. In order to to this
# in ggbpubr a function for removing non-significant comparisons had to be made (removeNS). 
# A further two graphs will be produced if the user has defined a control, only significant
# comparisons with the control will be displayed on the graph. 



library(ggpubr)

# get argument passed to command line (location to set as working directory)
args <- commandArgs(trailingOnly = TRUE)
# if file name has any spaces the filepath is split at each space, so it is pasted back into
# one string so it can be used as a filepath. 
location <- paste(args, collapse = " ")
setwd(location)

subdirectories <- c()
neurite.average <- c()
myelination.average <- c()
myelination.raw <- c()
experiment.name <- c()
neurite.raw <- c()
allExperiment.name <- c()
subdirectories2  <- c()
control <- c()
remove = list()


getMax <- function(values)
{
  # get the maximum value from a list and add 10%
  # @param values numeric list
  # @return result numeric value for largest number in value + 10% 
  result <- max(values) + max(values)/10 
  return(result)
}


removeNS <- function(statistics.analysis)
{
  # Remove any statical comparisons that are not significant (p < 0.05, adjusted p value if
  # multiple comparisons).
  # @param readout list of statistical analysis using ggpubr
  # @return readout list of only significant (p <0.05) statistical comparisons. 
  for (i in 1:nrow(statistics.analysis[,5])){
    p <- statistics.analysis[i,5]
    if (p >= 0.05){
      remove <- append(remove,i)
    } 
  }

  # all comparisons with p >= 0.05 are removed from stats table
  if (length(remove < 0)){
    statistics.analysis <- statistics.analysis[-as.numeric(remove),]
  }
  
  return(statistics.analysis)
 
}


getComparisons <- function(signif)
{
  # ggpubr requires a list of all comparisons to be made. 
  # @param signif list of all signigicant (p < 0.05) comparisions. 
  # @return comp 2D list of all significant comparisions. 
  comp <- list()
  if(nrow(signif) > 0){
    col1 <- signif$group1
    col2 <- signif$group2
    for (i in 1:nrow(signif)){
      comp[[i]] <- c(col1[i], col2[i])
    }
  }
  return(comp)
}


getSize <- function(values){
  # Determine values for graph width and column width
  # @param values list all values to be plotted.
  # @ return width numeric list for width and column width fo graph required.
  
  numchar <- 0
  count <- 0
  
  # get unique names (i.e name of each experimental condition).
  experiment.conditions <- (unique(values[1]))
  
  # determine width of graph based on number of letters in each name. 
  for(i in 1:length(experiment.conditions[[1]])){
    numchar <- as.numeric(nchar(as.character(experiment.conditions[i, 1])))
    if (numchar > 24){
      num <- numchar - 24
      count <- count + num
    }
  }
  
  count <-  count * 0.25
  count2 <- length(experiment.conditions[[1]]) * 2
  width <- count + count2
  
  # determine column width. 
  if (width > 4){
    colwidth <- (width - 4)/94
    colwidth <- 0.6 - colwidth
  } else{
    colwidth <- 0.6
  }
 
  width <- append(width, colwidth)
  return(width)
  
}


statsGraph <- function(results, comparisons, comps2, ylab2, max2, size)
{
  
  # make bar graph using ggpubr and save as a .tiff using the nature publishing group (NPG)
  # palette from ggsci. 
  # @param results dataframe for values to be plotted and experiment name.
  # @param comparisons 2D list for all significant (p < 0.05) comparisons. 
  # @param comps2 string data.frame column name for results to plot (in results). 
  # @param ylab2 string for y axis label. 
  # @param max2 numeric for largest value to be plotted.
  # @param sizes numeric list for graph and column width. 
  
  tiff(paste(comps2, "tiff", sep = "."), height = 5, width = size[1], units = 'in', res = 600, compression = 'lzw')
  p <- ggbarplot(results, x = "experiment.name", y = comps2, xlab = FALSE, width = size[2], legend = "none", add = c("mean_se", "jitter"), color = "experiment.name", ylab = ylab2)+
    stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "fdr", label = "p.signif")
  set_palette(p, "NPS")
  ggpar(p, ylim = c(0,30))
  print(p) # inside loop or function have to use print. Otherwise graph will not be saved!!
  dev.off()
}

statsGraph.control <- function(results, comparisons, comps2, ylab2, max2, size, control)
{
  
  # make bar graph using ggpubr and save as a .tiff using the nature publishing group (NPG)
  # palette from ggsci. 
  # @param results dataframe for values to be plotted and experiment name.
  # @param comparisons 2D list for all significant (p < 0.05) comparisons. 
  # @param comps2 string data.frame column name for results to plot (in results). 
  # @param ylab2 string for y axis label. 
  # @param max2 numeric for largest value to be plotted.
  # @param size numeric list for graph and column width. 
  # @param control string defines which is experimental conditions is the control. 
  
  comps3 <- paste(comps2, "compared to control only", sep = "")
  tiff(paste(comps3, "tiff", sep = "."), height = 5, width = size[1], units = 'in', res = 600, compression = 'lzw')
  p <- ggbarplot(results, x = "experiment.name", y = comps2, xlab = FALSE, width = size[2], legend = "none", add = c("mean_se", "jitter"), color = "experiment.name", ylab = ylab2)+
    stat_compare_means(comparisons = comparisons, method = "t.test", p.adjust.method = "fdr", label = "p.signif")
  set_palette(p, "NPS")
  ggpar(p, ylim = c(0,30))
  print(p) # inside loop or function have to use print. Otherwise graph will not be saved!!
  dev.off()
}


distribution.graph <- function(results,comps2, ylab2, max2, size)
{
  comps3 <- paste(comps2, "distribution", sep = "")
  tiff(paste(comps3, "tiff", sep = "."), height = 3.08, width = size[1], units = 'in', res = 600, compression = 'lzw')
  p <- ggbarplot(results, x = "allExperiment.name", y = comps2, ylab = ylab2, legend = "none", xlab = FALSE, width = size[2],
                 color = "allExperiment.name",
                 add = c("jitter", "mean_se"))
  set_palette(p, "NPS")
  ggpar(p, ylim = c(0,max2))
  print(p) # inside loop or function have to use print. Otherwise graph will not be saved!!
  dev.off()
}

# Get a list of subdirectories and .csv files
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get the name of all subdirectories (each subdirectory is one experimental conditionsition,
# which contains a .csv file for each experimental repeat). 
for (subdir in list.dirs(recursive=FALSE)) {
  subdirectories <- as.character(append(subdirectories, subdir))
}
# get a list of subdictories without ./ in the string (gets in the way)
subdirectories <- sapply(strsplit(subdirectories, split='/', fixed=TRUE), function(x) (x[2]))

# for each subdirectory get a list of all csv files (2D list)
csv.list <- list()
for (i in 1:length(subdirectories)){
  csv.list[[i]] = list.files(subdirectories[i], pattern="*.csv")
}


# determine whether the user has specific a control - if (control) in name.
# remove "(control)" from name so it does not appear in the graph.
# make a new list of subdirectories with the corrected name

for (i in 1:length(subdirectories)){
  if (grepl("control", subdirectories[i]) == TRUE){
    control <- strsplit(subdirectories[i], split = "[(]")[[1]]
    control <- control[1]
    subdirectories2 <- subdirectories
    subdirectories2[i] <- control
  }
}


for (i in 1:length(csv.list)){
  if (length(control) != 0){
  experiment.name <- append(experiment.name,rep(subdirectories2[i], length(csv.list[[i]]))) # make sure all of the means from one subdirectory will have the same name (as they are the same experimental conditionsition)
  } else{
    experiment.name <- append(experiment.name,rep(subdirectories[i], length(csv.list[[i]])))
  }
  
  # for each subdirectory get mean myelination and neurite density for each .csv file
  for (x in 1:length(csv.list[[i]])){
      data <- read.csv(paste(subdirectories[i], csv.list[[i]][[x]], sep = "/"), header = TRUE)
      data <- data[-1]
      neurite.average <- append(neurite.average, mean(as.numeric(data[1,])))
      neurite.raw <- append(neurite.raw, as.numeric(data[1,]))
      myelination.average <- append(myelination.average, mean(as.numeric(data[2,])))
      myelination.raw <- append(myelination.raw, as.numeric(data[2,]))
      allExperiment.name <- append(allExperiment.name, rep(csv.list[[i]][[x]], ncol(data)))
  }
}

# find the largest number and add 10% - this will be used in the graph so the y axis is correct

myelin.Max <- getMax(myelination.average)
neurite.Max <- getMax(neurite.average)
  
# make results the correct format for ggpubr i.e one list for means and one list for experimental conditionsition
myelinaverage.ggpubr <- data.frame(experiment.name)
myelinaverage.ggpubr["myelination.average"] <- as.numeric(myelination.average)
  
neuriteaverage.ggpubr <- data.frame(experiment.name)
neuriteaverage.ggpubr["neurite.average"] <- as.numeric(neurite.average)

myelinRaw.ggpubr <- data.frame(allExperiment.name)
myelinRaw.ggpubr["myelination.average"] <- as.numeric(myelination.raw)

neuriteRaw.ggpubr <- data.frame(allExperiment.name)
neuriteRaw.ggpubr["neurite.average"] <- as.numeric(neurite.raw)

# Perform statistical analysis

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(csv.list) > 1){ # safegard, needs to be at least two experiments 
  
  # Perform statistical analysis
  myelinstats <- compare_means(data = myelinaverage.ggpubr, myelination.average ~ experiment.name , method = "t.test", p.adjust.method = "fdr", paired = FALSE)
  neuritestats <- compare_means(data = neuriteaverage.ggpubr, neurite.average ~ experiment.name , method = "t.test", p.adjust.method = "fdr", paired = FALSE)
  # Save resutls from statistical analysis as .csv file
  write.csv(myelinstats, file = "%myelination statisticsSummary.csv",row.names=FALSE, na="")
  write.csv(neuritestats, file = "%neuritedensity statisticsSummary.csv",row.names=FALSE, na="")
  
  # is user has defined a control then all experimental conditionsitions will be compared to 
  # control only
  if (length(control) != 0){
    control.myelinstats <- compare_means(data = myelinaverage.ggpubr, myelination.average ~ experiment.name , method = "t.test", p.adjust.method = "fdr", ref.group = control, paired = FALSE)
    control.neuritestats <- compare_means(data = neuriteaverage.ggpubr, neurite.average ~ experiment.name , method = "t.test", p.adjust.method = "fdr", ref.group = control, paired = FALSE)
    # Save resutls from statistical analysis as .csv file
    write.csv(control.myelinstats, file = "%myelination statisticsSummary compared to control only.csv",row.names=FALSE, na="")
    write.csv(control.neuritestats, file = "%neuritedensity statisticsSummary compared to control only.csv",row.names=FALSE, na="")
  }
      
  # so that all comparisons are not put on the graph (which can make it very busy) 
  # comparisons with p values > 0.05 are removed.
  
  myelinstats <- removeNS(myelinstats)
  neuritestats <- removeNS(neuritestats)
  
  if (length(control) != 0){
     control.myelinstats <- removeNS(control.myelinstats)
     control.neuritestats <- removeNS(control.neuritestats)
  }
  
  myelinsignif.comparisons <- getComparisons(myelinstats)
  neuritesignig.comparisons <- getComparisons(neuritestats)
  
  if (length(control) != 0){
      myelinsignif.comparisonsVscontrol <- getComparisons(control.myelinstats)
      neuritesignif.comparisonsVscontrol <- getComparisons(control.neuritestats)
  }
  # make bar graphs using T test to compare all experimental conditions
  myelingraph.size <- getSize(myelinaverage.ggpubr)
  statsGraph(myelinaverage.ggpubr, myelinsignif.comparisons, "myelination.average", "% myelination", myelin.Max, myelingraph.size)
  neuritegraph.size <- getSize(neuriteaverage.ggpubr)
  statsGraph(neuriteaverage.ggpubr, neuritesignig.comparisons, "neurite.average", "% neurite density", neurite.Max, neuritegraph.size)
  
  # make bar graphs using T test to compare all experimental conditions to control. 
  if (length(control) != 0){
    statsGraph.control(myelinaverage.ggpubr, myelinsignif.comparisonsVscontrol, "myelination.average", "% myelination", myelin.Max, myelingraph.size, control)
    statsGraph.control(neuriteaverage.ggpubr, neuritesignif.comparisonsVscontrol, "neurite.average", "% neurite density", neurite.Max, neuritegraph.size, control)
  }
}    
  
myelin.Max <- getMax(myelination.raw)
neurite.Max <- getMax(neurite.raw)

# make bar graph illustrating distribution of values for each image.
neurite.distributionSize <- getSize(neuriteRaw.ggpubr)
distribution.graph(neuriteRaw.ggpubr, "neurite.average", "% neurite density", neurite.Max, neurite.distributionSize)

myelin.distributionSize <- getSize(myelinRaw.ggpubr)
distribution.graph(myelinRaw.ggpubr, "myelination.average", "% myelination", myelin.Max, myelin.distributionSize)

