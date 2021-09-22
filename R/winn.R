# WiNN version 0.4

library(hwwntest)
library(gam)
library(stringr)
library(lawstat) # <= Form Leveve test. Temporary until we get the "car" library working

#####################################################
# Utility function for test.wn
# Find greatest power of 2 that is closest to "n"
# ("closest.pw2", could be > or < "n") and the greatest
# power of 2 that is closest to "n" and < "n"
# ("greatest.pw2")
#####################################################
greatest.closest.power.2 <- function(n) {
  i <- 0
  
  repeat {
    if (n >= 2 ^ i & n < 2 ^ (i + 1)) {
      break
    }
    else{
      i <- i + 1
    }
  }
  
  i.closest <- i
  
  # Is this really the power od 2 closest to "n"?
  # (n=255 would return a "power of 2" i value = 7, while it is closer to 2^8)
  if (abs(n - 2 ^ (i + 1)) < (n - 2 ^ i)) {
    i.closest <- (i + 1)
  }
  
  return(list(closest.pw2 = i.closest, greatest.pw2 = i))
}


######################################
# White Noise test
# The test combines hwwntest and Box
# tests
######################################
test.wn <- function(x, box.test.lag = 20) {
  # The p-value for WN to return
  p <- NA
  
  # Remove NAs
  x.no.na <- x[!is.na(x)]
  
  # The WN test can be applied for series lenghts >= 16
  if (length(x.no.na) >= 16) {
    # Determine greatest power of 2 which is < length of data
    z <- greatest.closest.power.2(length(x.no.na))[["greatest.pw2"]]
    
    # Compute delta
    delta <- length(x.no.na) - 2 ^ z
    
    # If delta < 20 , delete "delta" values randomly from the series
    p <- c()
    if (delta < 20) {
      if (delta > 0) {
        x.no.na <- x.no.na[-sample(1:length(x.no.na), delta)]
      }
      
      # Perform the test and return p.value for w.noise (Ho: WN)
      p <- genwwn.test(x.no.na)[["p.value"]]
    } else{
      #print("here 1")
      
      # Create 2 series of length 2^z:
      #   - from position 1 to 2^z
      #   - from position (n - 2^z) to n
      x.no.na.1 <- x.no.na[1:2 ^ z]
      n <- length(x.no.na)
      x.no.na.2 <- x.no.na[(n - 2 ^ z + 1):n]
      
      p1 <- genwwn.test(x.no.na.1)[["p.value"]]
      p2 <- genwwn.test(x.no.na.2)[["p.value"]]
      
      p <- p1
      if (p2 <= p1) {
        p <- p2
      }
    }
  } else{
    print("Error: WN test cannot applied, length of series < 16")
  }
  
  # Add the Box test
  p.box <- Box.test(x, box.test.lag)$p.value
  
  # Return the smallest of the 2 p-values
  p.ret <- p
  if (!is.na(p.box) &
      !is.na(p) &
      p.box <= p) {
    p.ret <- p.box
  } #print(" ===> Selected p.box")}
  
  if (is.na(p.ret)) {
    stop("WN test cannot applied, length of series < 16")
  }
  
  return(p.ret)
}

###################
# smoothing spline
###################
apply.gam.hastie <- function(y, fun = "spline", k) {
  x <- 1:length(y)
  
  # Use spline or loess in gam
  fun.type <- "s"
  if (fun == "loess") {
    fun.type -> "lo"
  }
  #form <- paste("y ~ ", fun.type,"(x)", sep="") s(x, k = -1, bs = "cs"
  form <- paste("y ~ ", fun.type, "(x, ", k , ")", sep = "")
  m <- gam(as.formula(form), na.action = na.exclude)
  
  return(predict(m, newdata = data.frame(x = x)))
}


########
# Plot
########
plot.uncorr.corr <- function(orig.m,
                             unresid.metab,
                             signal,
                             corrected.signal,
                             seq.order,
                             plate.coord,
                             n.plates,
                             file.name,
                             is.wn,
                             p,
                             k)
{
  correction        <- signal - corrected.signal
  range  <- c(0, as.numeric(plate.coord[length(plate.coord)]))
  
  pdf(file.name, width = 15, height = 10)
  layout(matrix(
    c(1, 2, 3),
    nrow = 3,
    ncol = 1,
    byrow = TRUE
  ))
  
  plot(
    seq.order,
    unresid.metab,
    type = "l",
    xlim = range,
    col = "black",
    xlab = "Reading Sequence",
    ylab = "Measurement",
    main = paste(orig.m, ": Untransformed")
  )
  
  # Add vertical lines demarcating the plates
  draw.plates <- function(n, plate.coord, y.pos) {
    start = 0
    for (j in 1:n) {
      abline(v = plate.coord[j], lty = 2)
      text(
        x = start + (plate.coord[j] - start) / 2 ,
        y = y.pos,
        j,
        col = "blue",
        cex = 1.2
      )
      # paste("Plate:", j), col = "blue", cex = 1.2)
      start <- plate.coord[j]
    }
  }
  draw.plates(
    n = n.plates,
    plate.coord = as.numeric(plate.coord),
    y.pos = max(unresid.metab) * 0.85
  )
  
  plot(
    seq.order,
    unlist(signal),
    type = "l",
    xlim = range,
    col = "darkgrey",
    xlab = "Reading Sequence",
    ylab = "Residuals",
    main = paste(orig.m, ": Normalized by plate sd + residualized by plate")
  )
  draw.plates(
    n = n.plates,
    plate.coord = as.numeric(plate.coord),
    y.pos = max(signal) * 0.85
  )
  lines(seq.order, unlist(correction), col = "red")
  legend(
    "bottomleft",
    lty = c(1, 1),
    col = c("darkgrey", "red"),
    legend = c("Residualized", "Correction")
  )
  
  title <-
    paste("Corrected (WN Test p-val:",
          round(p, 4),
          " - ",
          is.wn,
          " - Best smooth. param.:",
          k,
          ")")
  plot(
    seq.order,
    corrected.signal,
    type = "l",
    xlim = range,
    xlab = "Reading Sequence",
    ylab = "Corrected Residuals",
    main = paste(orig.m, ":", title)
  )
  draw.plates(
    n = n.plates,
    plate.coord = as.numeric(plate.coord),
    y.pos = max(corrected.signal) * 0.85
  )
  
  dev.off()
}

plot.uncorr.corr.2 <- function(orig.m,
                               unresid.metab,
                               signal,
                               corrected.signal,
                               seq.order,
                               plate.coord,
                               n.plates,
                               file.name,
                               is.wn,
                               p,
                               k,
                               code)
{
  plot.uncorrected <- FALSE
  if (is.na(corrected.signal) &
      is.na(is.wn) & is.na(p) & is.na(k)) {
    plot.uncorrected <- TRUE
  }
  
  if (!plot.uncorrected) {
    correction <- signal - corrected.signal
  }
  
  range  <- c(0, as.numeric(plate.coord[length(plate.coord)]))
  
  pdf(file.name, width = 15, height = 10)
  layout(matrix(
    c(1, 2, 3),
    nrow = 3,
    ncol = 1,
    byrow = TRUE
  ))
  
  plot(
    seq.order,
    unresid.metab,
    type = "l",
    xlim = range,
    col = "black",
    xlab = "Reading Sequence",
    ylab = "Measurement",
    main = paste(orig.m, ": Untransformed")
  )
  
  # Add vertical lines demarcating the plates
  draw.plates <- function(n, plate.coord, y.pos) {
    start = 0
    for (j in 1:n) {
      abline(v = plate.coord[j], lty = 2)
      text(
        x = start + (plate.coord[j] - start) / 2 ,
        y = y.pos,
        j,
        col = "blue",
        cex = 1.2
      )
      # paste("Plate:", j), col = "blue", cex = 1.2)
      start <- plate.coord[j]
    }
  }
  draw.plates(
    n = n.plates,
    plate.coord = as.numeric(plate.coord),
    y.pos = max(unresid.metab) * 0.85
  )
  
  ###############################
  ylab <- c()
  title <- c()
  if (code == "00") {
    ylab = "Measurements (no sd. norm., no resid)"
    title = "No sd. normalization, no residualization"
  }
  if (code == "10") {
    ylab = "Measurements - sd. normalized"
    title = "Normalized by plate SD, no residualization"
  }
  if (code == "01") {
    ylab = "Measurements - no sd. normalized, residualized)"
    title = "Not normalized by plate SD, residualized"
  }
  if (code == "11") {
    ylab = "Measurements (sd. normalized + residualized)"
    title = "Normalized by plate SD and residualized by plate"
  }
  
  plot(
    seq.order,
    unlist(signal),
    type = "l",
    xlim = range,
    col = "darkgrey",
    xlab = "Reading Sequence",
    ylab = ylab,
    main = paste(orig.m, ":", title)
  )
  draw.plates(
    n = n.plates,
    plate.coord = as.numeric(plate.coord),
    y.pos = max(signal) * 0.85
  )
  
  if (!plot.uncorrected) {
    lines(seq.order, unlist(correction), col = "red")
    legend(
      "bottomleft",
      lty = c(1, 1),
      col = c("darkgrey", "red"),
      legend = c("Before correction", "Correction")
    )
    
    
    #################################################
    title <-
      paste("Corrected (WN Test p-val:",
            round(p, 4),
            " - ",
            is.wn,
            " - Best smooth. param.:",
            k,
            ")")
    plot(
      seq.order,
      corrected.signal,
      type = "l",
      xlim = range,
      xlab = "Reading Sequence",
      ylab = "Corrected Residuals",
      main = paste(orig.m, ":", title)
    )
    draw.plates(
      n = n.plates,
      plate.coord = as.numeric(plate.coord),
      y.pos = max(corrected.signal) * 0.85
    )
  }
  dev.off()
}


###########################
# Validate input dataset
###########################
validate.input <- function(df, group.name) {
  # Return value
  ret <- NA
  
  # Verify df is matrix or df
  is.matrix.or.df <- TRUE
  
  if (!is.matrix(df) & !is.data.frame(df)) {
    stop("\n- \"met.dat\" must be a matrix or a data frame")
    #quit(save="no")
  }
  
  # All Missing values ?
  mis <- unlist(lapply(df, function(x) {
    return(sum(is.na(x)))
  }))
  if (length(mis[mis == 0]) == 0) {
    stop("\n- missing values detected in ALL columns of \"met.dat\"")
  }
  
  flag.rownames                  <- FALSE
  flag.rownames.group.name       <- FALSE
  flag.rownames.group.number     <-
    FALSE # is.numeric("group numbers")?
  flag.rownames.group.number.seq <-
    FALSE # is.numeric("group numbers") AND in sequential order?
  flag.rownames.order.name       <- FALSE # is this == "order"
  flag.rownames.order.number     <-
    FALSE # is.numeric("order numbers")?
  flag.rownames.order.number.seq <-
    FALSE # is.numeric("order numbers") AND in sequential order?
  
  # Verify rownames structure
  z <- sapply(rownames(df), function(x) {
    return(strsplit(x, "_"))
  })
  
  # (1) verify that the length of each vector in the list z is
  # either 4 ("plate_<>_order_<>") or 6 ("plate_<>_order_<>_id_<>")
  z.l <- unlist(lapply(z, function(x) {
    return(length(x))
  }))
  if (!(sum(z.l) == 4 * length(z.l) | sum(z.l) == 6 * length(z.l))) {
    flag.rownames <- TRUE
  }
  
  #### If the row names are wrong STOP here and return an error
  # stopifnot(!flag.rownames)
  if (flag.rownames) {
    stop(
      "\n - wrong row names format must be: \"<group name>_<group number>_order_<order number>\" or \"<group name>_<group number>_order_<order number>_id_<participant id>\""
    )
  }
  
  # (2) Validate fields
  group.names     <- unlist(lapply(z, function(x) {
    return(x[1])
  }))
  group.numbers   <- unlist(lapply(z, function(x) {
    return(x[2])
  }))
  order.names     <- unlist(lapply(z, function(x) {
    return(x[3])
  }))
  order.numbers   <- unlist(lapply(z, function(x) {
    return(x[4])
  }))
  
  # Check group names are all the same
  if (sum(group.names == group.name) != length(group.names)) {
    flag.rownames.group.name <- T
  }
  
  # Check the "group numbers" are actually alphanumeric
  is.num <-
    sapply(group.numbers, function(x) {
      return(!grepl("^[A-Za-z]+$", x , perl = T))
    })
  if (sum(is.num) != length(group.numbers)) {
    flag.rownames.group.number <- T
  }
  # print(flag.rownames.group.number)
  
  # Check that the group numbers are in sequential order
  if (!flag.rownames.group.number) {
    if (is.unsorted(as.numeric(group.numbers))) {
      flag.rownames.group.number.seq <- T
    }
  } else{
    flag.rownames.group.number.seq <- T
  }
  
  # Check "order" names
  if (sum(order.names == "order") != length(order.names)) {
    flag.rownames.order.name <- T
  }
  
  # Check "order numbers" are alphanumeric
  is.num <-
    sapply(order.numbers, function(x) {
      return(!grepl("^[A-Za-z]+$", x , perl = T))
    })
  if (sum(is.num) != length(order.numbers)) {
    flag.rownames.order.number <- T
  }
  # print(flag.rownames.order.number)
  if (!flag.rownames.order.number) {
    if (is.unsorted(as.numeric(order.numbers))) {
      flag.rownames.order.number.seq <- T
    }
  } else{
    flag.rownames.order.number.seq <- T
  }
  
  
  # Check that the order numbers are in sequential order
  # if(is.unsorted(as.numeric(order.numbers))){flag.rownames.order.number.seq <- T}
  
  messages <- c()
  if (flag.rownames.group.name) {
    messages <-
      c(
        messages,
        paste(
          "- wrong group_name field in one or more row names (expected:\"",
          group.name,
          "\")",
          sep = ""
        )
      )
  }
  if (flag.rownames.group.number) {
    messages <-
      c(messages,
        "- group_number field not alphanumeric in one or more row names")
  }
  if (flag.rownames.group.number.seq) {
    messages <- c(messages, "- group_number(s) not sorted in row names")
  }
  
  if (flag.rownames.order.name) {
    messages <-
      c(messages,
        "- wrong order_name field in one or more row names (expected: \"order\")")
  }
  if (flag.rownames.order.number) {
    messages <-
      c(messages,
        "- order number field not alphanumeric in one or more row names")
  }
  if (flag.rownames.order.number.seq) {
    messages <-
      c(messages, "- order_number(s) not sorted in row names")
  }
  
  # Should we issue an error
  issue.error <-
    any(
      flag.rownames.group.name,
      flag.rownames.group.number,
      flag.rownames.group.number.seq,
      flag.rownames.order.name,
      flag.rownames.order.number,
      flag.rownames.order.number.seq
    )
  
  if (issue.error) {
    # Create stop message
    stop.message <- paste("\n", paste(messages, collapse = "\n"))
    stop(stop.message)
  }
  
  # If we made it so far, then we can return the dataset with complete
  ret <- df
  if (length(names(mis[mis == 0])) != dim(df)[2]) {
    complete <- names(mis[mis == 0])
    ret <- df[complete]
    print(
      paste(
        "Warning: found",
        length(complete),
        "complete (i.e. with no NA) variables out of",
        dim(df)[2]
      )
    )
  }
  return(ret)
}

#################################################
# Create a summary statistics for the PP dataset
# for each PP type
#################################################
summary.pp <- function(pp.type, pooled.plasma) {
  # Get the metabolite names
  met.names <- names(pooled.plasma)
  
  # Create the summary data frame to populate
  mean.sd.median <-
    data.frame(
      PP.type = rep(pp.type, length(met.names)),
      n = rep(NA, length(met.names)),
      not.miss = rep(NA, length(met.names)),
      metab = met.names,
      mean = rep(NA, length(met.names)),
      sd = rep(NA, length(met.names)),
      median = rep(NA, length(met.names)),
      max = rep(NA, length(met.names)),
      min = rep(NA, length(met.names)),
      cv = rep(NA, length(met.names)),
      stringsAsFactors = FALSE
    )
  
  ns <- c()
  means <- c()
  sds <- c()
  medians <- c()
  
  maxs <- c()
  mins <- c()
  not.nas <- c()
  
  # Select the PP data rows corresponding to the pp.type
  pp.types <-
    sapply(rownames(pooled.plasma), function(x) {
      strsplit(x, "_id_")[[1]][2]
    })
  dt <- pooled.plasma[pp.types == pp.type, ]
  
  # Get the "n"
  n.qc <- dim(dt)[1]
  mean.sd.median$n <- rep(n.qc, length(met.names))
  
  # Now loops through the metabolites and collect # of not missing, mean, sd, max, min, median
  for (j in met.names) {
    not.nas <- c(not.nas, n.qc - sum(is.na(dt[[j]])))
    means   <- c(means, mean(dt[[j]], na.rm = T))
    sds     <- c(sds,     sd(dt[[j]], na.rm = T))
    maxs    <- c(maxs,   max(dt[[j]], na.rm = T))
    mins    <- c(mins,   min(dt[[j]], na.rm = T))
    medians <- c(medians, median(dt[[j]], na.rm = T))
  }
  
  # Now add to the data frame
  mean.sd.median$mean     <- means
  mean.sd.median$sd       <- sds
  mean.sd.median$median   <- medians
  mean.sd.median$max      <- maxs
  mean.sd.median$min      <- mins
  mean.sd.median$not.miss <- not.nas
  mean.sd.median$cv       <-
    sds / means # Compute the CVs of the metabolites for this QC type
  
  return(mean.sd.median)
}




########################################################
# Transform Pooled Plasma data
#
# For each metabolite:
# (1) Get the mean baseline metabolite value: mean(met)
# (2) For each PP group:
#     Apply a linear transformation to the values in
#     the PP group so that the mean is the same as
#     "m.bsln" and the sd is such that CV (CV_PP) is preserved.
#
#      met_pp' = a*met_pp + b         (2.1)
#
#      We want:
#         - mean(met_pp') == mean(met_bsln)   (2.2)
#         - sd(met_pp') such that cv' == cv, ie
#               sd(met_pp')/mean(met_pp') == sd(met_pp)/mean(met_pp) (2.3)
#
#     It can be shown, that conditions (2.2) and (2.3) are satisfied
#     if the coeffs of the linear transformation (2.1) are:
#
#          - a = cv*mean(met_bsln)/sd(met_pp)
#          - b = mean(met_bsln) - (cv*mean(met_bsln)/sd(met_pp))*mean(met_pp)
#
#######################################################
transform.pp <- function(untrnsf.dt, pooled.plasma) {
  # A list for collecting the summary statistics relative to each PP type
  summary.list <- list()
  
  # Collect the summary statistics for each PP type
  # Get the PP types in the PP dataset
  pp.types <-
    unique(sapply(rownames(pp.dt), function(x) {
      strsplit(x, "_id_")[[1]][2]
    }))
  for (pp.type in pp.types) {
    summary.list[[pp.type]] <-
      summary.pp(pp.type = pp.type, pooled.plasma = pooled.plasma)
  }
  
  # Create a summary statistics data frame
  pp.stats.summ <- do.call("rbind", summary.list)
  str(pp.stats.summ)
  
  # Create a new PP data frame to store the scaled and shifted data set
  new.pp.dt <-
    data.frame(matrix(
      nrow = dim(pooled.plasma)[1],
      ncol = dim(pooled.plasma)[2]
    ))
  names(new.pp.dt)      <- names(pooled.plasma)
  rownames(new.pp.dt)   <- rownames(pooled.plasma)
  
  i <- 1
  for (met in names(untrnsf.dt)) {
    #(1) Get mean from the untransformed data
    mean.met.bsln <- mean(untrnsf.dt[[met]], na.rm = T)
    
    if (!i %% 100) {
      print(paste("Processed", i, "metabolites ..."))
    }
    i <- i + 1
    
    # Define vectors to store the transformed met
    new.met <- c()
    new.order <- c()
    
    for (pp.gr in pp.types) {
      # print(pp.gr)
      # Extract cv, mean and sd for this metabolite in the pp.gr PP type
      cv.pp <-
        pp.stats.summ[pp.stats.summ$PP.type == pp.gr &
                        pp.stats.summ$metab == met, "cv"]
      mn.pp <-
        pp.stats.summ[pp.stats.summ$PP.type == pp.gr &
                        pp.stats.summ$metab == met, "mean"]
      sd.pp <-
        pp.stats.summ[pp.stats.summ$PP.type == pp.gr &
                        pp.stats.summ$metab == met, "sd"]
      
      # Select the rows in the "met" column corresponding to this PP type
      pp.type.ids <-
        sapply(rownames(pooled.plasma), function(x) {
          strsplit(x, "_id_")[[1]][2]
        })
      tmp.dt <- pp.dt[pp.type.ids == pp.gr, met, drop = FALSE]
      
      # Coefficients of the linear transformation that would shift the data to the mean value
      # of the metabolite in the baseline experimental data AND scale so that the CV is preserved
      a <- cv.pp * mean.met.bsln / sd.pp
      b <- mean.met.bsln - a * mn.pp
      
      # Apply transformation
      new.met   <- c(new.met, a * tmp.dt[[met]] + b)
      
      # Crucial step - BE VERY CAREFUL HERE: keep track of the position in the sequence of wells!!
      seq.ord <-
        as.numeric(sapply(rownames(tmp.dt), function(x) {
          return(strsplit(strsplit(x, "_id_")[[1]][1], "_order_")[[1]][2])
        }))
      new.order <- c(new.order, seq.ord)
      
    } # End loop through the PP types
    
    names(new.met)   <- new.order
    new.met.ordered  <- new.met[order(new.order)]
    new.pp.dt[[met]] <- new.met.ordered
    
  } # End loop through the metabolites
  
  return(new.pp.dt)
}


#################################################
# Applying WiNN correction to Pooled Plasma data
#################################################
#' @title Apply WiNN correction to pooled plasma data and returns
#' CVs before and after correction.
#'
#' @description Applies WiNN correction to PP data and returns CVs
#' before and after correction.
#' @param winn.corrected The WiNN corrected metabolite dataset.
#' @param untrnsf.dt The uncorrected metabolite dataset (with *NO* pooled plasma in it).
#' @param pooled.plasma The Pooled Plasma dataset
#' @keywords correct.pp
#' @return list with CVs before and after correction
#' @export
#' @examples
#' correct.pp()
correct.pp <- function(winn.corrected,
                       untrnsf.dt,
                       pooled.plasma,
                       stdy.id) {

  #str(winn.corrected)
  
  # The transformed PP data
  pp.dt <- transform.pp(untrnsf.dt, pooled.plasma)
  
  cv.before <- c()
  cv.after  <- c()
  
  # FG: debug
  #str(pp.dt)
  # FG: debug
  
  # Get the wells sequence from rownames of the corrected data, which could be either
  # "plate_<>_order<order in sequence>_id<>" or just "plate_<>_order_<order in sequence>"
  
  corr.well.order <-
    as.numeric(unlist(lapply(strsplit(rownames(winn.corrected), "_id_"), function(x) {
      z <- x[[1]]
      if (length(z == 2)) {
        split <- z[1]
      } else{
        split <- z
      }
      return(strsplit(split, "_order_")[[1]][2])
    })))
  
  
  
  # Get the well sequence in the PP data
  pp.well.order <-
    as.numeric(sapply(rownames(pp.dt), function(x) {
      return(strsplit(strsplit(x, "_id_")[[1]][1], "_order_")[[1]][2])
    }))
  
  # Corrected
  all.corrected.met     <- names(winn.corrected)
  
  # loop through the mets
  for (met in all.corrected.met) {
   pp.met <- get.orig.met.name(met)
   
   pp.met.mn <-
      mean(pp.dt[[pp.met]], na.rm = T)     # Same as mean of experimental sample untranformed and uncorrected
   
    # Shift and scale the corrected experimental sample to its original mean and sd
    shift.corr.met <-
      (sd(untrnsf.dt[[pp.met]], na.rm = T) / sd(winn.corrected[[met]], na.rm = T)) *
      (winn.corrected[[met]] - mean(winn.corrected[[met]], na.rm = T)) + 
      mean(untrnsf.dt[[pp.met]], na.rm = T)
    
    #print("**********************************************")
    #print(met)
    #print(pp.met)
    #print(sd(untrnsf.dt[[pp.met]], na.rm = T))
    #print(sd(winn.corrected[[met]], na.rm = T))
    #print(mean(winn.corrected[[met]], na.rm=T))
    #print(mean(untrnsf.dt[[pp.met]], na.rm = T))
    #print("**********************************************")
    
    # Correction ratio
    corr.ratio <- as.numeric(shift.corr.met / untrnsf.dt[[pp.met]])
    
    # Interpolate
    pp.interp.ratio <-
      approx(corr.well.order, corr.ratio, xout = pp.well.order)
    new.pp <- pp.dt[[pp.met]] * pp.interp.ratio[["y"]]
   
    # Compute CV before and after
    cv.before.corr <-
      sd(pp.dt[[pp.met]], na.rm = T) / mean(pp.dt[[pp.met]], na.rm = T)
    cv.after.corr  <- sd(new.pp, na.rm = T) / mean(new.pp, na.rm = T)
    
    # Add to the before and after vectors
    cv.before <- c(cv.before, cv.before.corr)
    cv.after  <- c(cv.after,  cv.after.corr)
  }
  
  ret.df           <- data.frame(cv.before = cv.before, cv.after = cv.after)
  rownames(ret.df) <- all.corrected.met
  
  here <- paste(getwd(), "/", stdy.id, sep = "")
  pp.dir <- paste(here,"/pooled_plasma_cvs",sep="")
  dir.create(pp.dir)
  #ret.df["cv.improved"] <- ifelse(ret.df$cv.after < ret.df$cv.before, TRUE, FALSE)
  ret.df["cv.improved"] <- ifelse(ret.df$cv.after < ret.df$cv.before, TRUE, ifelse(ret.df$cv.after > ret.df$cv.before, FALSE, NA))
  write.table(ret.df, file=paste(pp.dir,"/pp_cv.csv",sep=""), quote=F, row.names = T, col.names = T, sep = "\t")
  
  return(ret.df)
}


##########################
# Pooled Plasma: CV plot
##########################
#' @title Creates before and after correction CV plot
#'
#' @description Creates before and after correction CV plot.
#' @param cv.pp The pooled plasma CV object returned by correct.pp
#' @param cv.thr The CV threshold value for the plot
#' @keywords cv.plot
#' @return
#' @export
#' @examples
#' cv.plot()
cv.plot <- function(cv.pp, cv.thr=NA, stdy.id){
  
  # Get the corrected ones
  cv.pp.corr <- cv.pp[str_starts(rownames(cv.pp), pattern="^detr."),]
  
  # Get the uncorrected but either normalized ("norm.<met>") or 
  # residualized ("resid.<met>")or normalized + residualized ("resid.norm.<met>")
  cv.pp.uncorr <- cv.pp[str_starts(rownames(cv.pp), pattern="^resid.") | 
                  str_starts(rownames(cv.pp), pattern="^norm.") ,]
 
  here <- paste(getwd(), "/", stdy.id, sep = "")
  pp.dir <- paste(here,"/pooled_plasma_cvs",sep="")
  dir.create(pp.dir)

  if(dim(cv.pp.corr)[1] > 0){
    pdf(paste(pp.dir,"/cvs_before_after_detrend.pdf",sep=""),  width = 12, height = 5)
    create.cv.plot(cv.before=cv.pp.corr$cv.before, cv.after=cv.pp.corr$cv.after, cv.thr=cv.thr, 
                 title="Corrected")
    dev.off()
  }
  if(dim(cv.pp.uncorr)[1] > 0){
    pdf(paste(pp.dir,"/cvs_before_after_no_detrend.pdf",sep=""),  width = 12, height = 5)
    create.cv.plot(cv.before=cv.pp.uncorr$cv.before, cv.after=cv.pp.uncorr$cv.after, cv.thr=cv.thr, 
                 title="Uncorrected (norm or resid or norm+resid)")
    dev.off()
  }
  
  # Create Delta-CVs boxplots
  pdf(paste(pp.dir,"/delta_cvs.pdf",sep=""),  width = 8, height = 8)
  delta.no.detrend <- cv.pp.uncorr$cv.after - cv.pp.uncorr$cv.before
  delta.detrend    <- cv.pp.corr$cv.after - cv.pp.corr$cv.before
  
  cv.all <- data.frame(cv.delta=c(delta.no.detrend, delta.detrend), detrended=c(rep(0, length(delta.no.detrend)),
                                                                    rep(1, length(delta.detrend))))

  if(sum(cv.all$detrended) == dim(cv.all)[1]){
    boxplot(cv.all$cv.delta, outline=F, 
            names=c("Delta CVs: with detrend"), xlab="", ylab="Delta CVs(after - before)")
    q.25.d <- round(quantile(delta.detrend)["25%"],2)
    q.50.d <- round(quantile(delta.detrend)["50%"],2)
    q.75.d <- round(quantile(delta.detrend)["75%"],2)
    txt.d  <- paste(q.50.d,"(", q.25.d,"-", q.75.d,")",sep="")
    
    legend("topright", legend=txt.d, cex=0.6, bty="n")
    
  }else{
    boxplot(cv.delta ~ detrended, data=cv.all, outline=F,
            names=c("Delta CVs: no detrend", "Delta CVs: with detrend"), xlab="", ylab="Delta CVs(after - before)")
  
    q.25.n <- round(quantile(delta.no.detrend)["25%"],2)
    q.50.n <- round(quantile(delta.no.detrend)["50%"],2)
    q.75.n <- round(quantile(delta.no.detrend)["75%"],2)
    txt.n  <- paste(q.50.n,"(", q.25.n,"-", q.75.n,")",sep="")
  
    q.25.d <- round(quantile(delta.detrend)["25%"],2)
    q.50.d <- round(quantile(delta.detrend)["50%"],2)
    q.75.d <- round(quantile(delta.detrend)["75%"],2)
    txt.d  <- paste(q.50.d,"(", q.25.d,"-", q.75.d,")",sep="")
  
    legend("topright", legend=c(txt.n,txt.d), cex=0.6, bty="n")
  }
  dev.off()
  
}



######################################################
# Pooled Plasma: CV plot
#   - Creates before and after correction CV plot
#   - cv.before CVs before correction.
#   - cv.after CVs after correction.
#   - cv.thr The CV threshold value for the plot
#   - title The plot title
######################################################
create.cv.plot <- function(cv.before, cv.after, cv.thr, title) {
  
  col.after <- "goldenrod1"
  col.before <- "skyblue3"
  
  layout(matrix(
    c(1, 2),
    nrow = 1,
    ncol = 2,
    byrow = TRUE
  ))
  
  
  plot(
    sort(cv.before),
    1:length(cv.before),
    type = "l",
    col = col.before,
    ylab = "Number of metabolites",
    xlab = "CV",
    cex.axis = 0.9,
    lwd = 3,
    main = title
  )
  
  lines(sort(cv.after),
        1:length(cv.after),
        col = col.after,
        lwd = 3)
  
  cv.thr.leg.b <- NA
  cv.thr.leg.a <- NA
  symb <- NA
# Draw CV threshold dashed lines
  if(!is.na(cv.thr)){
    after.x.coord.y.eq.thr <-
      length(sort(cv.after)[sort(cv.after) < cv.thr])
    segments(
      x0 = cv.thr,
      y0 = 0,
      y1 = after.x.coord.y.eq.thr,
      x1 = cv.thr,
      lty = 2
    )
    segments(
      y0 = after.x.coord.y.eq.thr,
      x0 = 0,
      y1 = after.x.coord.y.eq.thr,
      x1 = cv.thr,
      lty = 2
    )
    before.x.coord.y.eq.thr <-
      length(sort(cv.before)[sort(cv.before) < cv.thr])
    segments(
      y0 = before.x.coord.y.eq.thr,
      x0 = 0,
      y1 = before.x.coord.y.eq.thr,
      x1 = cv.thr,
      lty = 2
    )
  
    # Add 2 intersection points
    points(
      cv.thr,
      before.x.coord.y.eq.thr,
      pch = 4,
      lwd = 3,
      cex = 1.5,
      col = col.before
    )
    points(
      cv.thr,
      after.x.coord.y.eq.thr,
      pch = 4,
      lwd = 3,
      cex = 1.5,
      col = col.after
    )
    cv.thr.leg.b <- paste("# metab before correction with cv <=", 
                          cv.thr,":",before.x.coord.y.eq.thr )
    cv.thr.leg.a <- paste("# metab after correction with cv <=",
                          cv.thr, ":", after.x.coord.y.eq.thr)
    
    symb <- 4 # this is "X"
  }
### End draw CV threshold
  legend(
    "bottomright",
    lty = c(1, 1, NA, NA),
    pch = c(NA, NA, symb, symb),
    bty = "n",
    legend = c(
      "Before Correction",
      "After Correction",
      cv.thr.leg.b,
      cv.thr.leg.a
    ),
    col = c(col.before, col.after),
    cex = 0.9,
    lwd = c(3, 3)
  )
  
  
  # Add delta_cv boxplots
  cvs <- data.frame(cvs=c(cv.before, cv.after), before.after=
                      as.factor(c(rep(0,length(cv.before)),rep(1,length(cv.after)))))
  
  boxplot(cvs ~ before.after, data=cvs, outline=F, names=c("Before Correction", "After Correction"), xlab="", col=c(col.before,col.after))
  
  if(TRUE){
    q.25.a <- round(quantile(cv.after)["25%"],2)
    q.50.a <- round(quantile(cv.after)["50%"],2)
    q.75.a <- round(quantile(cv.after)["75%"],2)
    txt.a <- paste(q.50.a,"(", q.25.a,"-", q.75.a,")",sep="")
  
    q.25.b <- round(quantile(cv.before)["25%"],2)
    q.50.b <- round(quantile(cv.before)["50%"],2)
    q.75.b <- round(quantile(cv.before)["75%"],2)
    txt.b <- paste(q.50.b,"(", q.25.b,"-", q.75.b,")",sep="")
  
    legend("topright", legend=c(txt.b,txt.a), fill=c(col.before, col.after), cex=0.6, bty="n")
  }
}

##############################################
#
# 
#
##############################################
cv.plot.2 <- function(cv.pp, cv.thr=NA, stdy.id){
  
  # Get the corrected ones
  cv.pp.corr <- cv.pp[str_starts(rownames(cv.pp), pattern="^detr."),]
  
  # Get the uncorrected but either normalized ("norm.<met>") or 
  # residualized ("resid.<met>")or normalized + residualized ("resid.norm.<met>")
  cv.pp.uncorr <- cv.pp[str_starts(rownames(cv.pp), pattern="^resid.") | 
                          str_starts(rownames(cv.pp), pattern="^norm.") ,]
  
  here <- paste(getwd(), "/", stdy.id, sep = "")
  pp.dir <- paste(here,"/pooled_plasma_cvs",sep="")
  dir.create(pp.dir)
  
  
  if(dim(cv.pp.corr)[1] > 0 | dim(cv.pp.uncorr)[1] > 0){
    pdf(paste(pp.dir,"/cvs_before_after.pdf",sep=""),  width = 10, height = 10)
    n.row <- 3
    n.col <- 3
  
    if(dim(cv.pp.corr)[1] == 0 | dim(cv.pp.uncorr)[1] == 0){
      n.row <- 2
      n.col <- 3
    }
    
    layout(matrix(
      1:(n.row*n.col),
     nrow = n.row,
      ncol = n.col,
      byrow = TRUE ))
    
    if(dim(cv.pp.corr)[1] > 0){
      create.cv.plot.2(cv.before=cv.pp.corr$cv.before, cv.after=cv.pp.corr$cv.after, cv.thr=cv.thr, 
                   title="With Detrend")
      # dev.off()
    }
    if(dim(cv.pp.uncorr)[1] > 0){
      #pdf(paste(pp.dir,"/cvs_before_after_no_detrend.pdf",sep=""),  width = 12, height = 5)
      create.cv.plot.2(cv.before=cv.pp.uncorr$cv.before, cv.after=cv.pp.uncorr$cv.after, cv.thr=cv.thr, 
                   title="No detrend (norm or resid or norm+resid)")
      # dev.off()
    }
    
    create.cv.plot.2(cv.before=cv.pp$cv.before, cv.after=cv.pp$cv.after, cv.thr=cv.thr, 
                   title="All (with detrend + no detrend)")
    dev.off()
  } 
}


######################################################
# Pooled Plasma: CV plot
#   - Creates before and after correction CV plot
#   - cv.before CVs before correction.
#   - cv.after CVs after correction.
#   - cv.thr The CV threshold value for the plot
#   - title The plot title
######################################################
create.cv.plot.2 <- function(cv.before, cv.after, cv.thr, title) {
  
  col.after <- "goldenrod1"
  col.before <- "skyblue3"
  
  n <- length(cv.before)
  plot(
    sort(cv.before),
    1:length(cv.before),
    type = "l",
    col = col.before,
    ylab = "Number of metabolites",
    xlab = "CV",
    cex.axis = 0.9,
    lwd = 3,
    main = title
  )
  
  lines(sort(cv.after),
        1:length(cv.after),
        col = col.after,
        lwd = 3)
  
  cv.thr.leg.b <- NA
  cv.thr.leg.a <- NA
  symb <- NA
  # Draw CV threshold dashed lines
  if(!is.na(cv.thr)){
    after.x.coord.y.eq.thr <-
      length(sort(cv.after)[sort(cv.after) < cv.thr])
    segments(
      x0 = cv.thr,
      y0 = 0,
      y1 = after.x.coord.y.eq.thr,
      x1 = cv.thr,
      lty = 2
    )
    segments(
      y0 = after.x.coord.y.eq.thr,
      x0 = 0,
      y1 = after.x.coord.y.eq.thr,
      x1 = cv.thr,
      lty = 2
    )
    before.x.coord.y.eq.thr <-
      length(sort(cv.before)[sort(cv.before) < cv.thr])
    segments(
      y0 = before.x.coord.y.eq.thr,
      x0 = 0,
      y1 = before.x.coord.y.eq.thr,
      x1 = cv.thr,
      lty = 2
    )
    
    # Add 2 intersection points
    points(
      cv.thr,
      before.x.coord.y.eq.thr,
      pch = 4,
      lwd = 3,
      cex = 1.5,
      col = col.before
    )
    points(
      cv.thr,
      after.x.coord.y.eq.thr,
      pch = 4,
      lwd = 3,
      cex = 1.5,
      col = col.after
    )
    cv.thr.leg.b <- paste("# metab before correction with cv <=", 
                          cv.thr,":",before.x.coord.y.eq.thr )
    cv.thr.leg.a <- paste("# metab after correction with cv <=",
                          cv.thr, ":", after.x.coord.y.eq.thr)
    
    symb <- 4 # this is "X"
  }
  ### End draw CV threshold
  legend(
    "bottomright",
    lty = c(1, 1, NA, NA),
    pch = c(NA, NA, symb, symb),
    bty = "n",
    legend = c(
      paste("Before Correction (n=", n,")", sep=""),
      "After Correction",
      cv.thr.leg.b,
      cv.thr.leg.a
    ),
    col = c(col.before, col.after),
    cex = 0.9,
    lwd = c(3, 3)
  )
  
  
  # Add boxplots CV before and after 
  cvs <- data.frame(cvs=c(cv.before, cv.after), before.after=
                      as.factor(c(rep(0,length(cv.before)),rep(1,length(cv.after)))))
  
  boxplot(cvs ~ before.after, data=cvs, outline=F, names=c("Before Correction", "After Correction"), xlab="", col=c(col.before,col.after))
  
  # Add legend with median CV and IQR before and after correction
  q.25.a <- round(quantile(cv.after)["25%"],2)
  q.50.a <- round(quantile(cv.after)["50%"],2)
  q.75.a <- round(quantile(cv.after)["75%"],2)
  txt.a <- paste(q.50.a,"(", q.25.a,"-", q.75.a,")",sep="")
    
  q.25.b <- round(quantile(cv.before)["25%"],2)
  q.50.b <- round(quantile(cv.before)["50%"],2)
  q.75.b <- round(quantile(cv.before)["75%"],2)
  txt.b <- paste(q.50.b,"(", q.25.b,"-", q.75.b,")",sep="")
    
  legend("topright", legend=c(txt.b,txt.a), fill=c(col.before, col.after), cex=0.9, bty="n")
  
  # Create boxplot of Delta CVs
  cvs[["delta"]] <- cv.after - cv.before
  boxplot(cvs$delta, outline=F, xlab="Delta CVs", ylab="Delta CVs(after - before)")
  q.25.d <- round(quantile(cvs[["delta"]])["25%"],2)
  q.50.d <- round(quantile(cvs[["delta"]])["50%"],2)
  q.75.d <- round(quantile(cvs[["delta"]])["75%"],2)
  txt.d  <- paste(q.50.d,"(", q.25.d,"-", q.75.d,")",sep="")
    
  legend("topright", legend=txt.d, cex=0.9, bty="n")
}


#################################################################
# Returns named boolean vector
# - TRUE : failed to reject homogeneity of variance
# - FALSE: reject homogeneity of variance
#################################################################
homogen.var <- function(met.dat)
{
  # Create a vector to store results of homogeneity variance test
  is.hom.var <- rep(FALSE, dim(met.dat)[2])
  names(is.hom.var) <- names(met.dat)
  
  #str(met.dat)
  #print(rownames(met.dat))

  for (i in names(met.dat)) {
    #print(i)
    met.df <-  met.dat[i]
    met.df[["group"]] <- sapply(rownames(met.df),
                                function(x) {
                                  strsplit(x, "_")[[1]][2]
                                })
    # Shapiro test of normality
    x <- shapiro.test(met.df[, i])
    
    # If normality  use Levine test for homogeneity of variance
    # Else us Fligner-Killeen test
    p.val <- c()
    if (x$p.value < 0.05) {
      # no normality : Fligner-Killeen
      flig.tst <-
        fligner.test(as.formula(paste(i, " ~ group", sep = "")), data = met.df)
      p.val <- flig.tst$p.value
    } else{
      # Change this once "car" library is available
      lev.tst <- levene.test(met.df[, i],
                             met.df[, "group"],
                             location = "median",
                             correction.method = "zero.correction")
      
      p.val <- lev.tst$p.value
    }
    
    if (p.val > 0.05) {
      is.hom.var[i] <- TRUE
    }
  }
  return(is.hom.var)
}


################################################
# Residualize
################################################
residualize <- function(met.df, met, group.var, plates) {
  for (i in 1:(length(plates) - 1)) {
    on.plate <- paste("in.", i, sep = "")
    met.df[[on.plate]] <- 0
    met.df[grep(paste(group.var, "_", plates[i], "_", sep = ""),
                rownames(met.df)),][[on.plate]] <- 1
  }
  
  # Formula
  form <-
    paste(met, " ~ ", paste(paste(
      "as.factor(in.", 1:(length(plates) - 1), ")", sep = ""
    ), collapse = " + "))
  
  # Run regression and residualize if no errors
  ret <- tryCatch({
    reg <- lm(as.formula(form), data = met.df, na.action = na.exclude)
    
    # Calling "resid" on the regression object created with the na.action=na.exclude
    # option takes care of the missing values. "resid(reg)" return NA if the
    # observation is also NA (remember that reg$resid will only return the residuals for the
    # non-missing values)
    list(res = resid(reg), status = 0)
  },
  error = function(e) {
    print(paste("Error processing", met))
    print(e)
    
    # Residualization failed, return a vector of "res.failed"
    return(list(res = rep("res.failed", dim(met.df)[1]), status = -1))
  }) # End of try-catch
  return(ret)
}

######################################################
# Determine number of independent dimensions by PCA
######################################################
do.pca <- function(res.norm.met.dt) {
  # PCA
  n.90.pc <- tryCatch({
    pc <- prcomp(res.norm.met.dt, scale = TRUE)
    cum.prop    <- summary(pc)$importance["Cumulative Proportion",]
    print(paste("PCA - number of dimensions:", min(which(cum.prop >= 0.9))))
    min(which(cum.prop >= 0.9))
  }, warning = function(w) {
    print(paste("WARNING @ PCA:", w))
    return(NA)
  }, error = function(e) {
    print(paste("ERROR @ PCA", e))
    return(NA)
  })
  
  ret.n.90.pc <- c()
  if (is.na(n.90.pc)) {
    ret.n.90.pc <- dim(res.norm.met.dt)[2]
  }
  else{
    ret.n.90.pc <- n.90.pc
  }
  
  return(ret.n.90.pc)
}

#######################################
# Returns original name of metabolite
######################################
get.orig.met.name <- function(met) {
  met.original <- c()
  if(str_starts(string = met, pattern = "resid\\.norm\\.")) {
    met.original <-str_remove(string = met, pattern = "^resid\\.norm\\.")
  }else if(str_starts(string = met, pattern = "resid\\.")) {
    met.original <- str_remove(string = met, pattern = "^resid\\.")
  }else if(str_starts(string = met, pattern = "norm\\.")) {
    met.original <- str_remove(string = met, pattern = "^norm\\.")
  }else if(str_starts(string = met, pattern = "detr\\.resid\\.norm\\.")) {
    met.original <- str_remove(string = met, pattern = "^detr\\.resid\\.norm\\.")
  }else if(str_starts(string = met, pattern = "detr\\.resid\\.")) {
      met.original <- str_remove(string = met, pattern = "^detr\\.resid\\.")
  }else if(str_starts(string = met, pattern = "detr\\.norm\\.")) {
      met.original <- str_remove(string = met, pattern = "^detr\\.norm\\.")
  }else if(str_starts(string = met, pattern = "detr\\.")) {
      met.original <- str_remove(string = met, pattern = "^detr\\.")
  }else{
      met.original <- met
  }
  
  return(met.original)
}

########
# WiNN
########
#' @title WiNN in development
#'
#' @description A function for metabolite correction
#' @param input.dat input dataset as data frame or matrix
#' @param input.dat Study ID.
#' @param group.var Grouping variable. Defaults to "plate".
#' @param n.start.smooth Starting value for smoothing parameter. Defaults to 1.
#' @param n.stop.smooth  Last value for smoothing parameter. Defaults to 10.
#' @param debug Debug mode. Defaults to TRUE.
#' @param runall Apply correction regrdless of white noise test.Defaults to FALSE.
#' @param save.corrected.data Save corrected dataset. Defaults to TRUE
#' @param save.pdfs Save pdfs of corrected/uncorrected data. Defaults to TRUE
#' @return A list (corrected.strict, corrected.lenient,corrected.resid, correction.summary)
#' @keywords WiNN
#' @export
#' @examples
#' winn()
winn <-
  function(input.dat,
           stdy.id,
           group.var = "plate",
           n.start.smooth = 1,
           n.stop.smooth = 10,
           debug = T,
           runall = F,
           save.corrected.data = T,
           save.uncorrected.data = T,
           save.pdfs = T) {
    
    #### FOR DEVELOPMENT ONLY ##############
    if (FALSE) {
      dir.pp <-
        "/an/vital/vital200/QC/correction_gam_2020/JUNE_2020_replace_outliers/7_APPLY_TO_OTHER_DSETS/OTHER_DSETS/CVD_CACO"
      input.dat      <-
        get(load(paste(
          dir.pp, "/cvdcaco_met_uncorrected.RData", sep = ""
        )))
      pp.dt               <-
        get(load(paste(dir.pp, "/pp_no_outl.RData", sep = "")))
      stdy.id <- "CVD_CACO_12152020"
      group.var = "plate"
      n.start.smooth = 1
      n.stop.smooth = 10
      debug = T
      runall = F
      save.corrected.data = T
      save.uncorrected.data = T
      save.pdfs = T
    }
    ###########################
    
    print("Running WiNN with following parameters:")
    print(paste("met.dat =", deparse(substitute(input.dat))))
    print(paste("stdy.id =", stdy.id))
    print(paste("n.start.smooth =", n.start.smooth))
    print(paste("n.stop.smooth =", n.stop.smooth))
    print(paste("debug =", debug))
    print(paste("runall =", runall))
    print(paste("save.corrected.data =", save.corrected.data))
    print(paste("save.pdfs =", save.pdfs))
    print("########################################")
    
    ################################
    # Validate input
    ################################
    print("=> Validating input data...")
    met.dat <-
      tryCatch(
        validate.input(df = input.dat, group.name = group.var),
        error = function(e) {
          #message("An error occurred:\n", e)
          message("", e)
          quit(save = "default")
        },
        warning = function(w) {
          message("A warning occured:\n", w)
        }
      )
    
    #################################
    # Create study directory
    #################################
    if (save.corrected.data | save.pdfs) {
      here <- paste(getwd(), "/", stdy.id, sep = "")
      dir.create(here)
      
      # If save.pdf == TRUE, create also the plot directory in the study directory
      if (save.pdfs){
        pdfs.d <- paste(here, "/pdfs", sep = "")
      
        untr.d               <- paste(pdfs.d, "/untransformed", sep = "")
        sd.norm.only.d       <- paste(pdfs.d, "/normalized", sep = "")
        resid.only.d         <- paste(pdfs.d, "/residualized", sep = "")
        sd.norm.resid.d      <- paste(pdfs.d, "/residualized_normalized", sep = "")
        
        detr.untr.d          <- paste(pdfs.d, "/detrended_untransformed", sep = "")
        detr.sd.norm.only.d  <- paste(pdfs.d, "/detrended_normalized", sep = "")
        detr.resid.only.d    <- paste(pdfs.d, "/detrended_residualized", sep = "")
        detr.sd.norm.resid.d <- paste(pdfs.d, "/detrended_residualized_normalized", sep = "")
        
        dir.create(pdfs.d)
        dir.create(untr.d)
        dir.create(sd.norm.only.d)
        dir.create(resid.only.d)
        dir.create(sd.norm.resid.d)
        dir.create(detr.untr.d)
        dir.create(detr.sd.norm.only.d)
        dir.create(detr.resid.only.d)
        dir.create(detr.sd.norm.resid.d)
      }
    }
  
    ########################################################
    # Create a summary data frame to keep track of
    # the transformations/corrections performed on the data
    ########################################################
    n.mets <- dim(met.dat)[2]
    summary.transf <- data.frame(
      metabolite = names(met.dat),
      normalized = rep(FALSE, n.mets),
      residualized = rep(FALSE, n.mets),
      pass.wn.1 = rep(TRUE, n.mets),
      detrended = rep(FALSE, n.mets),
      pass.wn.2 = rep(NA, n.mets)
    )
    
    ########################################################
    # Step 1: Test for homogeneity of variance among plates
    ########################################################
    print("=> Testing homogeneity of variance and normalizing ...")
    is.hom.var <- homogen.var(met.dat)
    print(table(is.hom.var))
    
    ##############################################################
    # Step 2: Normalize by plate SD if homogeneity of variance
    # fails
    ##############################################################
    # Get metabs that do not pass hom. var. test
    mets.no.hom.var <- names(is.hom.var[!is.hom.var])
  
    # Loop through the metabolites and normalize by sd of groups (e.g. "plates")
    plates <-
      unique(unlist(lapply(strsplit(
        rownames(met.dat), "_"
      ), function(x) {
        return(x[2])
      })))
    
    for (i in mets.no.hom.var) {
      met.plates.sd <- c()
      for (j in plates) {
        grep.plates <-
          grep(paste(group.var, "_", j, "_", sep = ""), row.names(met.dat))
        sd.plate    <-  sd(met.dat[grep.plates, i], na.rm = T)
        if (sd.plate == 0) {
          sd.plate <- NA
        }
        met.plates.sd <-
          c(met.plates.sd, rep(sd.plate, length(grep.plates)))
      }
      
      # Normalize
      met.dat[[paste("norm.", i, sep = "")]] <- met.dat[[i]] / met.plates.sd
    
      # Update summary
      summary.transf[summary.transf$metabolite == i, "normalized"] <-TRUE
    }
    
    #############################################
    # Step 3:  Residualize by plate number
    #############################################
    # First loop through the original metabolite names and determine
    # which one should go through the ANOVA test and be potentially
    # residualized. If the metabolite has been normalized it is 
    # the norm.<original metab name>" that should be tested,
    # otherwise "<original metab name>":
    #
    # sd.normalized
    #       F             "<original metab name>"
    #       T             "norm.<original metab name>"
    print("=> Anova and residualization ...")
    metabs.to.anova.test <- c()
    for (i in as.character(summary.transf$metabolite)) {
      x <- summary.transf[as.character(summary.transf$metabolite) == i, ]
      if (x$normalized == F) {
        metab.to.test <- i
      }
      if (x$normalized == T) {
        metab.to.test <- paste("norm.", i, sep = "")
      }
      metabs.to.anova.test <- c(metabs.to.anova.test, metab.to.test)
    }
    
    # First let's determine the presence of a plate effect
    n.anova <- 0
    n.residualized <- 0
    for (met in metabs.to.anova.test) {
      # Extract the original nameif this metabolite was normalized
      met.original <- met
      if (str_starts(string = met, pattern = "norm.")) {
        met.original <- str_remove(string = met, pattern = "^norm\\.")
      }
      
      # Create the group variable for the ANOVA
      met.df            <- met.dat[met]
      met.df[["group"]] <-
        as.factor(paste(group.var, "_", unlist(lapply(strsplit(rownames(met.df), "_"),
          function(x) {return(x[2])})), sep = ""))
      
      # ANOVA test
      a <- aov(as.formula(paste(met, "~ group", sep = "")), data = met.df)
      a.pval <- summary(a)[[1]][, "Pr(>F)"][1]
      
      # Residualize if ANOVA test is significant
      if (a.pval < 0.05) {
        n.anova <- n.anova + 1
        residualized <- residualize(met.df, met, group.var, plates)
        
        # Was the residualization successfull, i.e. no errors?
        # If so, replace the metabolite value with  the residuals and
        # add the "resid" prefix to the name
        if (residualized$status == 0) {
          n.residualized <- n.residualized + 1
          met.dat[[paste("resid.", met, sep = "")]] <-residualized$res
          
          # Update summary
          summary.transf[summary.transf$metabolite == met.original, "residualized"] <-
            TRUE
        }
      }
    }
    
    print(paste(n.residualized, " metabolites have been residualized", sep=""))
    
    ##########################################################################
    # Step 4: WN test (first)
    ##########################################################################
    #
    # First loop through the original metabolite names and determine
    # which one should go through the WN test and be potentially
    # corrected. This is based on the transformation the metab
    # went through so far.
    #
    # sd.normilized residualized   metab
    #       F           F         "<original metab name>"
    #       F           T         "resid.<original metab name>"
    #       T           F         "norm.<original metab name>"
    #       T           T         "resid.norm.<original metab name>"
    #########################################################################
    metabs.to.wn.test <- c()
    for (i in as.character(summary.transf$metabolite)) {
      x <- summary.transf[as.character(summary.transf$metabolite) == i, ]
      if (x$normalized == F &
          x$residualized == F) {
        metab.to.test <- i
      }
      if (x$normalized == F &
          x$residualized == T) {
        metab.to.test <- paste("resid.", i, sep = "")
      }
      if (x$normalized == T &
          x$residualized == F) {
        metab.to.test <- paste("norm.", i, sep = "")
      }
      if (x$normalized == T &
          x$residualized == T) {
        metab.to.test <- paste("resid.norm.", i, sep = "")
      }
      metabs.to.wn.test <- c(metabs.to.wn.test, metab.to.test)
    }
    
    # Test for WN
    print("=> Testing for WN...")
    p.wn.tst <- sapply(metabs.to.wn.test, function(x) {
      ret <- NA; p <- NA
      p <- tryCatch({
        test.wn(met.dat[, x])
      }, warning = function(w) {
        print(paste("WARNING @ test.wn > ", x, ":", w))
      }, error = function(e) {
        print(paste("ERROR @ test.wn > ", x, ":", e))
        return(NA)
      })
      
      if (!is.na(p)) {
        ret <- p
      }
      
      return(ret)
    })
    
    names(p.wn.tst)  <- metabs.to.wn.test
    
    # Determine number of independent dimensions by PCA
    is.wn.before.corr <- c()
    
    
    # The number of independent dimensions
    print("=> Running first PCA ...")
    n.90.pc.before.corr <-
      do.pca(res.norm.met.dt = met.dat[metabs.to.wn.test])
    if (debug) {
      print(paste(
        "DBG 1 - number independent components:",
        n.90.pc.before.corr
      ))
    }
    
    # Set the p-value threshold to Bonferroni:
    p.thr.before.corr <- 0.05 / n.90.pc.before.corr
    if (debug) {
      print(paste("DBG 2 - p-value threshold:",  p.thr.before.corr))
    }
    
    # Determine which metab is WN (Ho: WN)
    is.wn.before.corr <-
      ifelse(p.wn.tst <= p.thr.before.corr, FALSE, TRUE)
    names(is.wn.before.corr) <- metabs.to.wn.test
    
    # Select the metabolites that cannot be tested for WN because
    # errors in the test for white noise function (eg length < 16)
    mets.cannot.test <-
      names(is.wn.before.corr[is.na(is.wn.before.corr)])
    if (length(mets.cannot.test) > 0) {
      print("Metabolites that could not be tested for WN:", log.file = cannot.correct)
      for (i in mets.cannot.test) {
        print(i, log.file = cannot.correct)
      }
    }
    
    # Select the metabolites that do *not* pass the WN test and must be corrected
    mets.to.correct <-
      names(is.wn.before.corr[!is.na(is.wn.before.corr) &
                                is.wn.before.corr == FALSE])
    
    print(paste(
      "Number of metabolites not passing WN to be corrected:",
      length(mets.to.correct)
    ))
    #for(i in mets.to.correct[1:5]){print(i)}
    
    # Perform a second PCA restricted to the metabs to be corrected
    print("=> Running second PCA ...")
    n.90.pc.after.corr <-
      do.pca(res.norm.met.dt = met.dat[mets.to.correct])
    if (debug) {
      print(
        paste(
          "DBG 3 - number independent components among mets to correct:",
          n.90.pc.after.corr
        )
      )
    }
    
    # n.90.pc.var <- 185
    p.thr.after.corr <- 0.05 / n.90.pc.after.corr
    if (debug) {
      print(paste(
        "DBG 4 - p-value threshold for mets to correct:",
        p.thr.after.corr
      ))
    }
    
    # Update the transf summary with the results of the WN test
    original.names <- c()
    for (i in mets.to.correct) {
      met.original <- c()
      if (str_starts(string = i, pattern = "resid\\.norm\\.")) {
        met.original <- str_remove(string = i, pattern = "^resid\\.norm\\.")
      } else if (str_starts(string = i, pattern = "resid\\.")) {
        met.original <- str_remove(string = i, pattern = "^resid\\.")
      } else if (str_starts(string = i, pattern = "norm\\.")) {
        met.original <- str_remove(string = i, pattern = "^norm\\.")
      } else{
        met.original <- i
      }
      
      original.names <- c(original.names, met.original)
      summary.transf[summary.transf$metabolite == met.original, "pass.wn.1"] <-
        FALSE
      
      # print()
    }
    
   # names(original.names) <- mets.to.correct
    
    ####################################
    # Step 5: apply correction
    ####################################
    plates <- as.numeric(plates)
    
    # If this is true, just run them all regardless of WN test
    if (runall) {
      mets.to.correct <- metabs.to.wn.test
      p.thr.after.corr <- 0.05 / (length(metabs.to.wn.test))
      n.90.pc.after.corr <- length(metabs.to.wn.test)
      
      # Also, we are going to correct them all, so let's update
      summary.transf["detrended"] <- rep(TRUE, n.mets)
    }
    
    # Post-correction WN test
    is.wn.post.corr <- c()
    count <- 0
    best.pvals <- c()
    
    # Get plate coordinates. It will be used later for plotting purposes
    plate.coord <- c()
    for (i in 1:length(plates)) {
      tmp.df <-
        met.dat[grep(paste(group.var, "_", i, "_order_", sep = ""),
                     rownames(met.dat),
                     value = T), ]
      plate.coord <-
        c(plate.coord, tail(sapply(rownames(tmp.df), function(x) {
          plate.order <-
            strsplit(x, "_id_")[[1]][1]
          return(strsplit(plate.order, "_order_")[[1]][2])
        }), n = 1))
    }
    
    # print(plate.coord)
    
    # Get also the sequence order for later
    seq.order <- as.numeric(sapply(rownames(met.dat), function(x) {
      plate.order <-
        strsplit(x, "_id_")[[1]][1]
      return(strsplit(plate.order, "_order_")[[1]][2])
    }))
    
    
    #q()
    
    
    # If we have metabolites to correct
    if (length(mets.to.correct) > 0) {
      # Get plate coordinates. It will be used later for plotting purposes
      print("=> Detrending ")
      
      if(FALSE){
      plate.coord <- c()
      for (i in 1:length(plates)) {
        tmp.df <-
          met.dat[grep(paste(group.var, "_", i, "_order_", sep = ""),
                       rownames(met.dat),
                       value = T), ]
        plate.coord <-
          c(plate.coord, tail(sapply(rownames(tmp.df), function(x) {
            plate.order <-
              strsplit(x, "_id_")[[1]][1]
            return(strsplit(plate.order, "_order_")[[1]][2])
          }), n = 1))
      }
      
      
      # Get also the sequence order for later
      seq.order <- as.numeric(sapply(rownames(met.dat), function(x) {
        plate.order <-
          strsplit(x, "_id_")[[1]][1]
        return(strsplit(plate.order, "_order_")[[1]][2])
      }))
    }
      #q()
      
      #### DEBUG
      n.dbg <- 0
      ##### DEBUG
      
      # Loop through the metabolites
      for (met in mets.to.correct) {
        # Get the original name
        #met.original <- original.names[met]
        met.original <- get.orig.met.name(met)
        
        # Extract the metabolite
        met.dt <- met.dat[met]
        
        print("##########################")
        print(paste("Detrending...", met))
        
        # Update count
        count <- count + 1
        if (count %% 20 == 0) {
          print(count)
        }
        
        pvals           <- c()
        z.best          <- c()
        p.best          <- 0
        k.best          <- c()
        pred.trend.best <- c()
        
        # Loop through the "smoothing parameter" values and select the
        # one which gives the greatest p-value for WN test (i.e. the one that
        # best detrends the signal to WN)
        for (k in n.start.smooth:n.stop.smooth) {
          # The vector that stores the trends for this smoothing parameter
          pred.trend <- c()
          
          # Loop through the plates
          for (i in plates) {
            # Subset to this plate
            plate.dt.met <-
              met.dt[grep(paste(group.var, "_", i, "_", sep = ""),
                          rownames(met.dt)), ]
            
            # Get the metabolite measurements in this plate
            y <- plate.dt.met
            
            # Get the trend for this metabolite for this plate and append it to pred.trend
            # ret.trend <- apply.gam.hastie(y=y, k=k)
            
            ret.trend <- tryCatch({
              apply.gam.hastie(y = y, k = k)
            }, warning = function(w, z = length(plate.dt.met)) {
              print(
                paste(
                  "WARNING @ apply.gam.hastie > ",
                  met,
                  " - smooth=",
                  k,
                  " - plate=",
                  i,
                  ":",
                  w
                )
              )
              return(rep(0, z))
            }, error = function(e, z = length(plate.dt.met)) {
              print(
                paste(
                  "ERROR @ apply.gam.hastie > ",
                  met,
                  " - smooth=",
                  k,
                  " - plate=",
                  i,
                  ":",
                  e
                )
              )
              return(rep(0, z))
            })
            
            # append to pred.trend
            pred.trend <- c(pred.trend, ret.trend)
          } # closes loop through plates
          
          # Correct the metabolite by subtracting  the trend relative to this
          # smoothing parameter
          z <- met.dt[[met]] - pred.trend
          
          # Get the white noise test pvalue relative to this smoothing parameter
          p.wn <- test.wn(z)
          
          # Add the WN p-value to the vector of WN test pvalues
          pvals <- c(pvals, p.wn)
          
          # Is this the "best" p-val so far? If so, set z.best, p.best, pred.trend.best
          if (p.wn >= p.best) {
            z.best <- z
            
            p.best <- p.wn
            
            k.best <- k
            
            pred.trend.best <- pred.trend
          }
        }  # closes for k in 1:n.smooth
        
        # Define a new "corrected metabolite" variable set to the "z.best" value obtained
        # after looping through the smoothing parameters
        corr.met <- paste("detr.", met, sep = "")
        met.dt[[corr.met]] <- z.best
        
        # Add the corrected metabolite to the original dset
        met.dat[[corr.met]] <- met.dt[[corr.met]]
        
        # Now compare the "p.best" value for this metabolite with the
        # threshold p.value determined in step 4
        result <- "Keep Ho: WN"
        if (p.best <= p.thr.after.corr) {
          result <- "Reject Ho: not WN"
        }
        print(paste(corr.met, " :  p.wn=", p.best, " - ", result, sep = ""))
        #q()
        # Add the test result to the vector of WN test results
        is.wn.post.corr <-
          c(is.wn.post.corr, p.best <= p.thr.after.corr)
        
        # Add the p.best value to the vector of "best" pvalues
        best.pvals <- c(best.pvals, p.best)
        
        # Finally update status
        # print(paste("HERE 1:", met.original, "HERE 1"))
        summary.transf[summary.transf$metabolite == met.original, "detrended"] <-
          TRUE
        # q()
        summary.transf[summary.transf$metabolite == met.original, "pass.wn.2"] <-
          p.best > p.thr.after.corr
        #q()
        
        #############
        # Plot
        #############
        save.pdfs <- TRUE
        if (save.pdfs) {
          #file.name <- paste(pdf.dir, "/corr.", met, ".pdf", sep = "")
          
          # Get the transformed metabolite that enters the correction ...
          code <-
            paste(as.numeric(summary.transf[summary.transf$metabolite == met.original, "normalized"]),
                  as.numeric(summary.transf[summary.transf$metabolite == met.original, "residualized"]),
                  sep = "")
          
          file.name <- c()
          if (code == "00") {
            file.name <- paste(detr.untr.d, "/", corr.met, ".pdf", sep = "")
          }
          if (code == "10") {
            file.name <-
              paste(detr.sd.norm.only.d, "/", corr.met, ".pdf", sep = "")
          }
          if (code == "01") {
            file.name <- paste(detr.resid.only.d, "/", corr.met, ".pdf", sep = "")
          }
          if (code == "11") {
            file.name <-
              paste(detr.sd.norm.resid.d, "/", corr.met, ".pdf", sep = "")
          }
        
          plot.uncorr.corr.2(
            orig.m = met.original,
            unresid.metab = met.dat[[met.original]],
            signal = met.dat[[met]], # <- this is what enters the correction step ....
            corrected.signal = met.dat[[corr.met]],
            seq.order = seq.order,
            plate.coord = plate.coord,
            n.plates = length(plates),
            file.name = file.name,
            is.wn = result,
            p = p.best,
            k = k.best,
            code = code
          )
        }
        
      } # Closes "for(met in ...)"
      
    } # Closes "if (length(mets.to.correct) > 0)"
    
    
    ##################################################
    # Create also pdfs of the metabolites that do not
    # go through correction but are either:
    # - uncorrected
    # - sd-plate normalized
    # - residualized
    # - sd-normalized and residualized
    #dir.create(pdfs.d)
    #dir.create(untr.d)
    #dir.create(sd.norm.only.d)
    #dir.create(resid.only.d)
    #dir.create(sd.norm.resid.d)
    ##################################################
    if (save.pdfs) {
      
      uncorr.plot <- function(mets,
                 summary.transf,
                 met.dat,
                 seq.order,
                 plate.coord,
                 plates) {
          #print(summary.transf)
          #print(mets)
          
          for (i in mets) {
            met.original.name <- get.orig.met.name(i)
            #print(i)
            #print(met.original.name)
            
            code <-
              paste(as.numeric(summary.transf[summary.transf$metabolite == met.original.name, "normalized"]),
                    as.numeric(summary.transf[summary.transf$metabolite == met.original.name, "residualized"]),
                    sep = "")
            
            file.name <- c()
            if (code == "00") {
              file.name <- paste(untr.d, "/", met.original.name, ".pdf", sep = "")
            }
            if (code == "10") {
              file.name <-
                paste(sd.norm.only.d, "/norm.", met.original.name, ".pdf", sep = "")
            }
            if (code == "01") {
              file.name <- paste(resid.only.d, "/resid.", met.original.name, ".pdf", sep = "")
            }
            if (code == "11") {
              file.name <-
                paste(sd.norm.resid.d, "/resid.norm.", met.original.name, ".pdf", sep = "")
            }
            
            plot.uncorr.corr.2(
              orig.m = get.orig.met.name(i),
              unresid.metab = met.dat[[get.orig.met.name(i)]],
              signal = met.dat[[i]], # <- this is what enters the correction step ....
              corrected.signal = NA,
              seq.order = seq.order,
              plate.coord = plate.coord,
              n.plates = length(plates),
              file.name = file.name,
              is.wn = NA,
              p = NA,
              k = NA,
              code = code
            )
          }
        }
      
      # Get the uncorrected metabs
      sel.untr         <-
        summary.transf$detrended == FALSE &
        summary.transf$normalized == FALSE &
        summary.transf$residualized == FALSE
      sel.norm         <-
        summary.transf$detrended == FALSE &
        summary.transf$normalized == TRUE  &
        summary.transf$residualized == FALSE
      sel.resid        <-
        summary.transf$detrended == FALSE &
        summary.transf$normalized == FALSE &
        summary.transf$residualized == TRUE
      sel.norm.resid   <-
        summary.transf$detrended == FALSE &
        summary.transf$normalized == TRUE  &
        summary.transf$residualized == TRUE
      
      uncorr.untransf.mets   <- summary.transf[sel.untr, "metabolite"]
      uncorr.norm.mets       <-
        paste("norm.", summary.transf[sel.norm, "metabolite"], sep = "")
      uncorr.resid.mets      <-
        paste("resid.", summary.transf[sel.resid, "metabolite"], sep = "")
      uncorr.norm.resid.mets <-
        paste("resid.norm.", summary.transf[sel.norm.resid, "metabolite"], sep =
                "")
      
      if (sum(sel.untr))
        uncorr.plot(
          mets = uncorr.untransf.mets,
          summary.transf = summary.transf,
          met.dat = met.dat,
          seq.order = seq.order,
          plate.coord = plate.coord,
          plates = plates
        )
      if (sum(sel.norm))
        uncorr.plot(
          mets = uncorr.norm.mets,
          summary.transf = summary.transf,
          met.dat = met.dat,
          seq.order = seq.order,
          plate.coord = plate.coord,
          plates = plates
        )
      if (sum(sel.resid))
        uncorr.plot(
          mets = uncorr.resid.mets,
          summary.transf = summary.transf,
          met.dat = met.dat,
          seq.order = seq.order,
          plate.coord = plate.coord,
          plates = plates
        )
      if (sum(sel.norm.resid))
        uncorr.plot(
          mets = uncorr.norm.resid.mets,
          summary.transf = summary.transf,
          met.dat = met.dat,
          seq.order = seq.order,
          plate.coord = plate.coord,
          plates = plates
        )
    }
    
    ############################################
    # Save summary
    ############################################
    correction.summary <- summary.transf
    save(
      list = c("summary.transf"),
      file = paste(here, "/correction_summary.RData", sep = "")
    )
    
    ############################################
    # Save transformed and corrected dataset
    ############################################
    # Loop through the metabolites and select the corresponding untransformed/transformed/corrected
    # columns according to the following criteria
    #
    #    corrected  sd.normalized   residualized    col name
    #       F           F               F           <original metab name>
    #       F           F               T           resid.<original metab name>
    #       F           T               F           norm.<original metab name>
    #       F           T               T           resid.norm.<original metab name>
    #       T           F               F           corr.<original metab name>
    #       T           F               T           corr.resid.<original metab name>
    #       T           T               F           corr.norm.<original metab name>
    #       T           T               T           corr.resid.norm.<original metab name>
    #
    
    col.names.sel <- c()
    for (i in summary.transf$metabolite) {
      x <- summary.transf[as.character(summary.transf$metabolite) == i, ]
      
      if (x$detrended == F &
          x$normalized == F & x$residualized == F) {
        col.name <- i
      }
      if (x$detrended == F &
          x$normalized == F &
          x$residualized == T) {
        col.name <- paste("resid.", i, sep = "")
      }
      if (x$detrended == F &
          x$normalized == T &
          x$residualized == F) {
        col.name <- paste("norm.", i, sep = "")
      }
      if (x$detrended == F &
          x$normalized == T &
          x$residualized == T) {
        col.name <- paste("resid.norm.", i, sep = "")
      }
      
      if (x$detrended == T &
          x$normalized == F &
          x$residualized == F) {
        col.name <- paste("detr.", i, sep = "")
      }
      if (x$detrended == T &
          x$normalized == F &
          x$residualized == T) {
        col.name <- paste("detr.resid.", i, sep = "")
      }
      if (x$detrended == T &
          x$normalized == T &
          x$residualized == F) {
        col.name <- paste("detr.norm.", i, sep = "")
      }
      if (x$detrended == T &
          x$normalized == T &
          x$residualized == T) {
        col.name <- paste("detr.resid.norm.", i, sep = "")
      }
      
      col.names.sel <- c(col.names.sel, col.name)
    }
    
    # The dataset to save and return
    transf.corrected <- met.dat[col.names.sel]
    save(
      list = c("transf.corrected"),
      file = paste(here, "/transf_corrected.RData", sep = "")
    )
    
    return(list(summary.transf = summary.transf, transf.corrected = transf.corrected))
  }
