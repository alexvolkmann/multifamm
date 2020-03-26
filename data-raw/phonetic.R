################################################################################
################################################################################
##                                                                            ##
##            Preparation of the Speech Production Data for Analysis          ##
##                                                                            ##
################################################################################
################################################################################


# Aim:
# Combine the raw data supplied for each dimension to a common data.frame.

################################################################################
# First domain: acoustic data
################################################################################

# Original data
load("prep_crossed_indiv_mirrored_s_sh_acoustic.Rdata")

# Extract data.table
phonetic_a <- out$curve_info
# phonetic_a$subject_long:   Subject ID
# phonetic_a$n_long:         Curve ID
# phonetic_a$combi_long:     Number of repetition for word
# phonetic_a$number_long:    Number of observations in curve
# phonetic_a$group_long:     0: s  -> sh
#                       1: sh -> s


################################################################################
# Second domain: EPG data
################################################################################

# Original data
load("prep_crossed_indiv_mirrored_s_sh_EPG.Rdata")
# Same structure as above

# Extract data.table
phonetic_e <- out$curve_info
# Same structure as above


################################################################################
# Combine domains
################################################################################

# Combine the data.tables and create identifier
phonetic <- rbind(phonetic_a, phonetic_e)
phonetic$dim <- factor(c(rep("aco", nrow(phonetic_a)),
                         rep("epg", nrow(phonetic_e))))

# Rename data.table for modeling
# covariate.1 : group_long
# covariate.2 : stress1_long
# covariate.3 : stress2_long
# covariate.4 : vowel_long
names(phonetic)[grep("group_long", names(phonetic))] <- "covariate.1"
names(phonetic)[grep("stress1_long", names(phonetic))] <- "covariate.2"
names(phonetic)[grep("stress2_long", names(phonetic))] <- "covariate.3"
names(phonetic)[grep("vowel_long", names(phonetic))] <- "covariate.4"
names(phonetic)[grep("ident_words_real_long",
                     names(phonetic))] <- "word_names_long"

# Reorder the variables
phonetic <- phonetic[, c("dim", "subject_long", "word_long", "combi_long",
                         "y_vec", "n_long", "t", "covariate.1", "covariate.2",
                         "covariate.3", "covariate.4", "word_names_long")]
usethis::use_data(phonetic)



################################################################################
# Create a subset
################################################################################

# Similar to acoustic_subset
index <- phonetic$subject_long %in% c(1, 2) & phonetic$word_long %in% c(1, 2) &
  phonetic$combi_long %in% 1:5
phonetic_subset <- phonetic[index, ]
usethis::use_data(phonetic_subset)
