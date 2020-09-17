### load necessary packages ###
### package used to create labels ###
library(calibrate)


### create function to plot local genetic covariances on manhattan style plot ###
### major changes: set significance threshold, allow highlighting significant regions, designate highlight color, spread chromosomes and correct log scales ###
supergnovaplot <- function(x, chr="CHR", bp="BP", rho="RHO", pval="PVAL",
                           pointcol="gray90", sigcol="red",
                           highlight=NULL, logp=TRUE, annotatePval = NULL, ...) {

  ### Set to NULL, will get warning without ###
  CHR=BP=RHO=index=NULL

  ### Make sure you have chr, bp and rho columns ###
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(rho %in% names(x))) stop(paste("Column", rho, "not found!"))
  ## warn if you don't have a region column
  #if (!(region %in% names(x))) warning(paste("No REGION column found. Needed for highlighting."))
  ## make sure chr, bp, and rho columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[rho]])) stop(paste(rho, "column should be numeric."))

  ### Add a region column to the input data frame ###
  ### column format is chr:start-end ###
  ### first two columns to paste together ###
  cols <- c( 'chr' , 'start')
  ### create a new column `regionstart` with chr and start pasted together, separated by ":" ###
  x$regionstart <- apply( x[ , cols ] , 1 , paste , collapse = ":" )
  ### next columns to paste together ###
  cols2 <- c( 'regionstart' , 'end')
  ### create "REGION" column by pasting together chr:start and end, separated by "-" ###
  x$region <- apply( x[ , cols2 ] , 1 , paste , collapse = "-" )

  ### create "pval" column, used to determine significance later ###
  x$pval <- x$p

  ### Create dataframe with
  d = data.frame(CHR=x[[chr]], BP=x[[bp]], RHO=x[[rho]], PVAL=x$pval, pos = NA, index = NA ,REGION=x$region, stringsAsFactors = FALSE)


  # Set positions, ticks, and labels for plotting
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$RHO)
  } else {
    d$logp <- d$RHO
  }


  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$REGION,d$CHR,length))  # replace the for loop of line 92-96 to improve efficiency

  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromosome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    #  ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    #  labs = ticks          ## unused, from code line: 169
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
        lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
        d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
        d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
        # lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        # d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase

      }
    }
    ticks <-tapply(d$pos,d$index,quantile,probs=0.5)   # replace line 135
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }

  ### Initialize plot ###
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  ### See http://stackoverflow.com/q/23922130/654296
  ### Set default arguments ###
  ### ylims are set to min/max roh values +/- 0.0001 to make sure the points and labels fit ###
  def_args <- list(xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20,
                   xlim=c(xmin,xmax), ylim=c(min(d$RHO)-0.0001,max(d$RHO)+0.0001),
                   xlab=xlabel, ylab=expression('local genetic covariance'))
  ### Get a list of ... arguments ###
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))


  # Add an axis.
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }


  ### Add all covariance estimates to the plot ###
  with(d, points(pos, logp, pch=20, col=pointcol, ...))



  ### Highlight regions with significant local covariance (p < set in annotatePval) ###
  topHits = subset(d, PVAL <= annotatePval)
  d.highlight=d[which(d$REGION %in% topHits$REGION), ]
  with(d.highlight, points(pos, RHO, col=sigcol, pch=20), ...)
  with(d.highlight, textxy(pos, RHO, offset = 0.57, labs=REGION, cex = 0.7))
  par(xpd = FALSE)
}


### location variables:
  ### chr = chromosome
  ### bp = start (I just used the start of the region because it was easy and the scale is small enough I figured it didn't make a difference)
### estimate variables
  ### rho = covariance estimates
  ### pval = covariance p value
### aesthetics
  ### pointcol = color of all estimates (default is gray)
  ### sigcol = color of regions with significant covariance (default is red)
### other arguments
  ### logp = keep set to false unless covariance estimates are log transformed
  ### annotatePval = significance threshold for p-value after bonferonni correction (no default)

traits <- c("B", "C", "D", "E")

for (k in traits) {
  #### read in data ####
  data_t <- read.table(paste(k,sep = "_","local_cov.txt"), header = TRUE)
  pdf(paste("A",k,"local_cov_plot.pdf",sep = "_"), width = 10, height = 5)
  supergnovaplot(data_t, chr="chr", bp="start", rho="rho", pval="p", cex.axis = 0.8, sigcol = "red", pointcol = "gray90", main = paste("Local Genetic Covariance between", k,"and Trait 1"), logp = FALSE, annotatePval = 0.00000076)
  dev.off()
}
