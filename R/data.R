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
#'   \item{\code{dim}}{Factor for identifying the acoustic (aco) and
#'      electropalatographic (epg) dimensions.}
#'   \item{\code{subject_long}}{Unique identification number for each speaker.}
#'   \item{\code{word_long}}{Unique identification number for each target word.}
#'   \item{\code{combi_long}}{Number of the repetition of the combination of the
#'      corresponding speaker and target word.}
#'  \item{\code{y_vec}}{The response values for each observation point.}
#'  \item{\code{n_long}}{Unique identification number for each curve.}
#'  \item{\code{t}}{The observations point locations.}
#'  \item{\code{covariate.1}}{Order of the consonants, reference category first
#'      /s/ then /sh/.}
#'  \item{\code{covariate.2}}{Stress of the final syllable of the first
#'      compound, reference category 'stressed'.}
#'  \item{\code{covariate.3}}{Stress of the initial syllable of the second
#'      compound, reference category 'stressed'.}
#'  \item{\code{covariate.4}}{Vowel context, reference category ia.}
#'  \item{\code{word_names_long}}{Names of the target words}
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


#' Snooker data
#'
#' The data are part of a study on the impact of a muscular training program on
#' snooker technique. 25 recreational snooker players were split into treatment
#' (receiving instructions for a training program) and control group (no
#' training program). The data set contains the movement trajectories of the
#' snooker players in two sessions (before and after the training period), where
#' each snooker player repeated a snooker shot of maximal force six times. The
#' interest lies in the movement of hand, elbow, and shoulder on a
#' two-dimensional grid (called X and Y). The trajectories are normalized on a
#' [0,1] time grid and the beginning of the hand trajectories are centered to
#' the origin.
#'
#' @format A data.frame with 56910 observations and 11 variables:
#' \describe{
#'  \item{\code{y_vec}}{The response values for each observation point.}
#'  \item{\code{t}}{The observations point locations.}
#'  \item{\code{n_long}}{Unique identification number for each curve.}
#'   \item{\code{subject_long}}{Unique identification number for each snooker
#'      player.}
#'   \item{\code{word_long}}{Integer specifying the session. 1: Before the
#'      training, 2: After the training.}
#'   \item{\code{dim}}{Factor for identifying the univariate dimensions.}
#'   \item{\code{combi_long}}{Number of the repetition of the snooker shot.}
#'  \item{\code{covariate.1}}{Skill level of the snooker player. 0: Unskilled,
#'      1: Skilled.}
#'  \item{\code{covariate.2}}{Group of the snooker player. 0: Control group,
#'      1: Treatment group.}
#'  \item{\code{covariate.3}}{Session indicator. 0: Before the treatment, 1:
#'      After the treatment.}
#'  \item{\code{covariate.4}}{Interaction of group and session, i.e. the
#'      treatment effect indicator.}
#' }
#'
#' @source Enghofer, T. (2014). Überblick über die Sportart Snooker, Entwicklung
#'    eines Muskeltrainings und Untersuchung dessen Einflusses auf die
#'    Stoßtechnik. Unpublished Zulassungsarbeit, Technische Universität
#'    München.
#'
"snooker"

