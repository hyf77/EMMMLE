#prepare function
geom_mean <- function(x, poscounts = TRUE, na.rm = TRUE) {

  x_log <- log(x)

  if (poscounts) x_log <- x_log[x > 0]

  exp(mean(x_log, na.rm = na.rm))

}

#offset function
offset_tss <- function(counts) {
  rowSums(counts)
}

offset_none <- function(counts) {
  return(rep(1,nrow(counts)))
}

offset_gmpr <- function(counts) {

  if (nrow(counts) == 1) stop("GMPR is not defined when there is only one sample.")

  ## median of pairwise ratios between counts of samples i and j, limited to positive counts

  pairwise_ratio <- function(i, j) {

    c_i <- counts[i, ]; c_j <- counts[j, ]

    ratio <- c_i / c_j

    median(ratio[c_i > 0 & c_j > 0])

  }

  ## Matrix of pairwise ratios

  n <- nrow(counts)

  mat_pr <- matrix(NaN, nrow = n, ncol = n)

  for (i in 1:n) {

    for (j in 1:n) {

      if (i != j) mat_pr[i, j] <- pairwise_ratio(i, j)

    }

  }

  ## Geometric mean of pairwise ratio

  size_factor <- apply(mat_pr, 1, geom_mean)

  if (any(size_factor == 0 | !is.finite(size_factor))) stop("Some sample(s) do not share any species with other samples, GMPR normalization failed.")

  size_factor

}

offset_css <- function(counts, reference = median) {

  ## special treatment for edge case of one-column matrix (1 OTU, many samples)

  if (ncol(counts) == 1) return( counts[, 1] / median(counts) )

  ## remove 0s and check that all samples have at least two positive counts

  counts[counts == 0] <- NA

  if (any(rowSums(!is.na(counts)) < 2)) {

    warning("Some samples only have 1 positive values. Can't compute quantiles and fall back to TSS normalization")

    return(rowSums(counts, na.rm = TRUE))

  }

  ## compute sample-specific quantiles and cumulative sums up to quantiles

  cumsum_up_to <- function(counts, quantiles) {

    colSums((counts * outer(counts, quantiles, `<=`)),na.rm = TRUE)

  }

  mat_sample_quant <- t(apply(counts, 1, quantile, probs = seq(0, 1, length.out = ncol(counts)), na.rm = TRUE))

  mat_sample_cumsum <- t(sapply(1:nrow(counts), function(i) { cumsum_up_to(counts[i, ], mat_sample_quant[i, ]) }))

  ## reference quantiles, computed as median (nature article) or mean (metagenomeSeq::cumNormStat[Fast]) of sample_specific quantiles

  ## and MAD around the reference quantiles

  ref_quant <- apply(mat_sample_quant, 2, reference)

  ref_quant_mad <- apply(abs(sweep(mat_sample_quant, 2, ref_quant)),2, median)

  ## find smallest quantile for which high instability is detected

  ## instability for quantile l is defined as ref_quant_mad[l+1] - ref_quant_mad[l] >= 0.1 * ref_quant_mad[l]

  instable <- (diff(ref_quant_mad) >= 0.1 * head(ref_quant_mad, -1))

  if (any(instable)) {

    ## Hack to mimick package implementation: never choose quantile below 50%

    lhat <- max(min(which(instable)), ceiling(ncol(counts)/2))

  } else {

    warning("No instability detected in quantile distribution across samples, falling back to scaled TSS normalization.")

    lhat <- ncol(counts)

  }

  ## scaling factors are cumulative sums up to quantile lhat, divided by their median

  size_factors <- mat_sample_cumsum[ , lhat] / median(mat_sample_cumsum[ , lhat])

  unname(size_factors)

}

############################
compute_offset <- function(counts, offset = c("TSS", "GMPR", "CSS", "none"), ...) {

  ## special behavior for data.frame

  if (inherits(offset, "data.frame")) {

    stop(

      "You supplied a data.frame to compute_offset(). Did you mean to supply a numeric matrix?

  Try converting your data.frame to a matrix with as.matrix()."

    )

  }

  ## special behavior for numeric offset

  if (is.numeric(offset)) {

    return(offset_numeric(counts, offset, ...))

  }

  ## Choose offset function

  offset <- match.arg(offset)

  offset_function <- switch(offset,

                            "TSS"    = offset_tss,

                            "GMPR"   = offset_gmpr,

                            "CSS"    = offset_css,

                            "none"   = offset_none

  )

  ## Ensure that counts is a matrix

  counts <- data.matrix(counts)

  ## Compute offset (with optional parameters)

  offset_function(counts, ...)

}
