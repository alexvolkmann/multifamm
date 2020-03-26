#' Phonetic data
#'
#' The data are part of a large study on consonant assimilation, which is
#' the phenomenon that the articulation of two consonants becomes
#' phonetically more alike when they appear subsequently in fluent speech.
#' The data set contains the audio signals of nine different speakers which
#' repeated the same sixteen German target words each five times. In addition to
#' these acoustic signals, the data set also contains the electropalatographic
#' data. The target words are bisyllabic noun-noun compound words which
#' contained the two abutting consonants of interest, s and sh, in either order.
#' Consonant assimilation is accompanied by a complex interplay of
#' language-specific, perceptual and articulatory factors. The aim in the study
#' was to investigate the assimilation of the two consonants as a function of
#' their order (either first s, then sh or vice-versa), syllable stress
#' (stressed or unstressed) and vowel context, i.e. which vowels are immediately
#' adjacent to the target consonants of interest. The vowels are either of the
#' form ia or ai. For more details, see references below.
#'
#' @format A data.frame with 50644 observations and 12 variables:
#' \describe{
#'   \item{\code{dim}}{factor for identifying the acoustic (aco) and
#'      electropalatographic (epg) dimensions.}
#'   \item{\code{subject_long}}{unique identification number for each speaker.}
#'   \item{\code{word_long}}{unique identification number for each target word.}
#'   \item{\code{combi_long}}{number of the repetition of the combination of the
#'      corresponding speaker and target word.}
#'  \item{\code{y_vec}}{the response values for each observation point.}
#'  \item{\code{n_long}}{unique identification number for each curve.}
#'  \item{\code{t}}{the observations point locations.}
#'  \item{\code{covariate.1}}{order of the consonants, reference category first
#'      /s/ then /sh/.}
#'  \item{\code{covariate.2}}{stress of the final syllable of the first
#'      compound, reference category 'stressed'.}
#'  \item{\code{covariate.3}}{stress of the initial syllable of the second
#'      compound, reference category 'stressed'.}
#'  \item{\code{covariate.4}}{vowel context, reference category ia.}
#'  \item{\code{word_names_long}}{names of the target words}
#' }
#'
#' @source Pouplier, Marianne and Hoole, Philip (2016): Articulatory and
#'   Acoustic Characteristics of German Fricative Clusters, Phonetica, 73(1),
#'   52--78.
#' @source Cederbaum, Pouplier, Hoole, Greven (2016): Functional Linear Mixed
#'   Models for Irregularly or Sparsely Sampled Data. Statistical Modelling,
#'   16(1), 67-88.
#' @source Jona Cederbaum (2019). sparseFLMM: Functional Linear Mixed Models for
#'   Irregularly or Sparsely Sampled Data. R package version 0.3.0.
#'    \url{https://CRAN.R-project.org/package=sparseFLMM}
#'
"phonetic"

#' Phonetic data (subset)
#'
#' A small subset of the phonetics data set \code{\link[multifamm]{phonetic}}
#' with observations from two speakers and two items only. This will not produce
#' meaningful results but can be used as a toy data set when testing the code.
#' The variables are as in the full data set, see
#' \code{\link[multifamm]{phonetic}}
#'
#' @format A data.frame with 1336 observations and 12 variables.
#'
#' @source Pouplier, Marianne and Hoole, Philip (2016): Articulatory and
#'   Acoustic Characteristics of German Fricative Clusters, Phonetica, 73(1),
#'   52--78.
#' @source Cederbaum, Pouplier, Hoole, Greven (2016): Functional Linear Mixed
#'   Models for Irregularly or Sparsely Sampled Data. Statistical Modelling,
#'   16(1), 67-88.
#' @source Jona Cederbaum (2019). sparseFLMM: Functional Linear Mixed Models for
#'   Irregularly or Sparsely Sampled Data. R package version 0.3.0.
#'    \url{https://CRAN.R-project.org/package=sparseFLMM}
"phonetic_subset"

