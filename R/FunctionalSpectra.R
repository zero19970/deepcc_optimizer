#' Preprocess Gene List
#'
#' This function preprocess gene list for futhur process
#'
#' @param geneList a named vector containing the value of gene expression
#' @return a named vecter containing the value of gene expression
#' @examples
#' preprocessGeneList(geneList)
preprocessGeneList <- function(geneList) {
  geneList <- geneList[which((!is.na(geneList)) & (names(geneList)!="") & (!is.na(names(geneList))))]
  geneList <- tapply(t(geneList), names(geneList), max)
  geneList[order(geneList, decreasing = TRUE)]
}

#' Calculate Enrichment Score
#'
#' This function calculates enrichment score of a gene list on a specific gene set.
#'
#' @param geneList a named vector containing the values of gene expression
#' @param geneSet a vector containing genes to represent a gene set
#' @return a numeric indicating enrichment score
#' @useDynLib DeepCCv2
#' @import Rcpp
#' @export
#' @examples
#' newCalcEnrichmentScore(geneList, geneSet)
newCalcEnrichmentScore <- function(geneList, geneSet)
{
  newCalcEnrichmentScoreCPP((names(geneList) %in% geneSet), geneList, 1)
}


#' Calculate single sample Enrichment Score
#'
#' This function calculates enrichment score of a gene list on a specific gene set.
#'
#' @param geneList a named vector containing the values of gene expression
#' @param geneSet a vector containing genes to represent a gene set
#' @return a numeric indicating enrichment score
#' @useDynLib DeepCCv2
#' @import Rcpp
#' @export
#' @examples
#' single_newCalcEnrichmentScore(geneList, geneSet)
single_newCalcEnrichmentScore <- function(geneList, geneSet)
  {
  single_newCalcEnrichmentScoreCPP((names(geneList) %in% geneSet), geneList, 1, times = 1000)
}
#' Generate Functional Spectra
#'
#' This function generates functional spectra for given gene expressin profiles.
#'
#' @param eps a data.frame containing gene expression profiles (each row presents one sample)
#' @param geneSet a List containing gene sets (default: MSigDB v6)
#' @param cores a integer indicating cpu cores used in parallel computing (default = all cores -2 )
#' @return a data.frame containing functional spectra
#' @param show_progress a bool value, to show the running progress, default: FALSE
#' @seealso  \code{\link{getFunctionalSpectrum}} for a single expression profile.
#' @importFrom foreach foreach %dopar%
#' @export
#' @examples
#' get_functional_spectra(eps)
get_functional_spectra <- function(eps, geneSets = 'MSigDBv7', scale = T, cores = parallel::detectCores() - 2) {
  if (geneSets == 'MSigDBv5') {
    data(MSigDBv5)
    geneSets = MSigDBv5
  } else if (geneSets == 'MSigDBv6') {
    data(MSigDBv6)
    geneSets = MSigDBv6
  } else if(geneSets == 'MSigDBv7') {
    data(MSigDBv7)
    geneSets = MSigDBv7
  }

  if(scale) eps <- scale(eps, scale = FALSE)

  doParallel::registerDoParallel(cores)

  res <- foreach(idx = 1:nrow(eps), .combine = rbind) %dopar% {
    geneList <- preprocessGeneList(eps[idx, ])
    if(show_progress) {
      sink("progress_log.txt.", append = T)
      cat("processing sample: ", idx, ":", rownames(eps)[idx], "finish: ", (idx-1)/nrow(eps), "\n")
    }
    sapply(geneSets, function(x) newCalcEnrichmentScore(geneList, x))
  }

  rownames(res) <- rownames(eps)
  res
}

#' Generate Functional Spectrum
#'
#' This function generates functional spectrum for a single gene expression profile.
#'
#' @param expressionProfile a named numeric vector containing gene expression profile
#' @param geneSets a List containing gene sets (default: MSiDB v6)
#' @return a numeric vector containing functional spectrum
#' @note You can generate the reference expression profile from your previous data or public data, which is the same(similiar) cancer type and platform.
#' In DeepCC we also prepared average expression profiles of each cancer types in TCGA project as references. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.
#' If your single sample is microarray data, we strongly sugguest turn the parameter \code{inverseRescale} on, since TCGA is RNA-Seq, which has very small expression value for low expressed genes, compared with microarray.
#' @seealso \code{\link{get_functional_spectra}} for a batch of gene expression profiles.
#' @export
#' @examples
#' getFunctionalSpectrum(ep)
getFunctionalSpectrum <- function(expressionProfile, geneSets = 'MSigDBv7', show_progress = F, No=NULL) {

  geneList <- preprocessGeneList(expressionProfile)

  if (geneSets == 'MSigDBv5') {
    data(MSigDBv5)
    geneSets = MSiDBv5
  } else if (geneSets == 'MSigDBv6') {
    data(MSigDBv6)
    geneSets = MSigDBv6
  } else if (geneSets == 'MSigDBv7') {
    data(MSigDBv7)
    geneSets = MSigDBv7
  }

  res <- sapply(1:length(geneSets), function(idx) {
    if(show_progress) {
      cat("sample:", paste0("No", No), "\t")
      cat("progressidx:", idx/length(geneSets)*100, "%\n")
    }

    single_newCalcEnrichmentScore(geneList, geneSets[[idx]])
  })
  names(res) <- names(geneSets)
  res
}


#' Generate Functional Spectrum for multi samples
#'
#' This function generates functional spectrum for a single gene expression profile for multiple samples.
#'
#' @param eps a data.frame containing gene expression profiles (each row presents one sample)
#' @param geneSets a List containing gene sets (default: MSiDB v6)
#' @param cores a integer indicating cpu cores used in parallel computing (default = all cores -2 )
#' @param show_progress a bool value, to show the running progress, default: FALSE
#' @return a numeric vector containing functional spectrum
#' @note You can generate the reference expression profile from your previous data or public data, which is the same(similiar) cancer type and platform.
#' In DeepCC we also prepared average expression profiles of each cancer types in TCGA project as references. To use them, just use the TCGA identifier (COADREAD, BRCA, OV, etc.) to indicate the cancer type.
#' If your single sample is microarray data, we strongly sugguest turn the parameter \code{inverseRescale} on, since TCGA is RNA-Seq, which has very small expression value for low expressed genes, compared with microarray.
#' @seealso \code{\link{get_functional_spectra}} for a batch of gene expression profiles.
#' @export
#' @examples
#' multi_sample_ssgsea(ep, geneSet = "MsigDBv7", show_progress = T)
multi_sample_ssgsea <- function(eps, geneSet = "MSigDBv7", cores = parallel::detectCores() - 2, show_progress = FALSE){
  doParallel::registerDoParallel(cores)

  res <- foreach(idx = 1:nrow(eps), .combine = rbind) %dopar% {
    geneList <- preprocessGeneList(eps[idx, ])
    if(show_progress) {
      sink("progress_log.txt.", append = T)
      cat("processing sample: ", idx, ":", rownames(eps)[idx], "finish: ", (idx-1)/nrow(eps), "\n")
      getFunctionalSpectrum(geneList, geneSet, show_progress, idx)
    } else {
      getFunctionalSpectrum(geneList, geneSet)
    }
  }


#  res = res / diff(res)
  rownames(res) <- rownames(eps)
  res
}
