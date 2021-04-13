# ==============================================================================
# 02_load_libraries_and_customfxns.R
# ==============================================================================



# load R packages
# ==============================================================================

# general i/o & data cleaning tools
library(data.table)
library(magrittr)
library(tidyverse)
library(readxl)
library(writexl)
library(stringi)
library(convertr)

# misc visualization
library(eulerr)
library(ComplexHeatmap)
library(cowplot)
library(ggConvexHull)
library(ggthemes)
library(ggsignif)
library(ggrepel)
library(png)
library(grid)
library(viridis)
library(circlize)

# bioconductor / RNA-seq analysis
library(biomaRt)
library(RUVSeq)
library(tximport)
library(edgeR)
library(GSVA)
library(WebGestaltR)
library(vegan)
library(ggbeeswarm)

# # also used for TvC comparisons (script 16)
# # commented out here because interfere with tidyverse (e.g., select)
# library(GEOquery)
# library(oligo)
# library(limma)



# some custom functions used throughout analysis
# ==============================================================================

summaryBlankLine <- function() {
  ""
}

summaryMedianIQR <-
  function(x, digits = 1, na.rm = F,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = T,
           symbol = NULL) {
    if (!(typeof(x) %in% c("double", "integer"))) {
      warning("input vector not numeric")
      x <- as.numeric(x)
    }


    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100

    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }

    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }

    output <- quantile(
      as.numeric(x),
      c(0.25, 0.5, 0.75),
      na.rm = na.rm
    )

    output <- formatC(output, format = "f", digits = digits)

    output <- ifelse(is.null(symbol),
      paste0(output[2], " (", output[1], ", ", output[3], ")"),
      paste0(output[2], symbol, " (", output[1], symbol, ", ", output[3], symbol, ")")
    )

    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }

    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }


    return(output)
  }

summaryMedianRange <-
  function(x, digits = 1, na.rm = F,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = T,
           symbol = NULL) {
    if (!(typeof(x) %in% c("double", "integer"))) {
      warning("input vector not numeric")
      x <- as.numeric(x)
    }

    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100

    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }

    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }

    output <- quantile(
      as.numeric(x),
      c(0, 0.5, 1),
      na.rm = na.rm
    )

    output <- formatC(output, format = "f", digits = digits)

    output <- ifelse(is.null(symbol),
      paste0(output[2], " (", output[1], ", ", output[3], ")"),
      paste0(output[2], symbol, " (", output[1], symbol, ", ", output[3], symbol, ")")
    )

    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }

    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }


    return(output)
  }

summaryCountPercent <-
  function(x, values, count_NAs = F,
           digits = 1,
           NA_percent_in_brackets = F,
           NA_count_in_brackets = T) {
    
    n_NA <- sum(is.na(x))
    ratio_miss <- n_NA / length(x) * 100

    if (n_NA != 0) {
      warning(paste0(
        n_NA,
        " values are missing/NA in input vector"
      ))
    }

    if (n_NA == length(x)) {
      warning("all values missing!")
      return("-")
    }

    n <- sum(x %in% values, na.rm = T)
    tot <- ifelse(count_NAs, length(x), length(x) - n_NA)

    output <- paste0(
      n, " (", formatC(n / tot * 100, format = "f", digits = digits),
      "%)"
    )

    if (NA_count_in_brackets & n_NA != 0) {
      output <-
        paste0(output, " [", n_NA, "]")
    }

    if (NA_percent_in_brackets & n_NA != 0) {
      output <-
        ratio_miss %>%
        formatC(., format = "f", digits = digits) %>%
        paste0(output, " [", ., "%]")
    }

    return(output)
  }


fxnMakeQuantiles <-
  function(prob = c(0.25, 0.5, 0.75)) {
    function(x) {
      tmp <- quantile(x, prob)
      tmp <- as.data.frame(t(tmp))
      names(tmp) <- c("ymin", "y", "ymax")
      tmp
    }
  }

fxnMedianLine <-
  function(x) {
    tmp <- rep(median(x), 3)
    tmp <- as.data.frame(t(tmp))
    names(tmp) <- c("ymin", "y", "ymax")
    tmp
  }




ggboxplotWrapper <- function(dataframe, xval, yval,
                             xlabel = gsub("`", "", xval), ylabel = yval,
                             seed = 1,
                             font_size = 16,
                             beesize = 2, beecex = 0.2) {
  set.seed(seed)

  ggplot(data = dataframe,
                  aes_string(x = xval, y = yval, color = xval, fill = xval)) +
    geom_quasirandom(alpha = 0.4, width = beecex) +
    stat_summary(fun.data = fxnMakeQuantiles(),
                 geom = "errorbar", color = "black",
                 width = 0.25, size = 0.8, alpha = 0.8) +
    stat_summary(fun.data = fxnMedianLine,
                 geom = "errorbar", color = "black",
                 width = 0.15, size = 1.5, alpha = 0.8) +
    xlab(xlabel) +
    ylab(ylabel) +
    theme_classic() +
    theme(text = element_text(size = font_size), legend.position = "none")
}



