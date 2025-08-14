#################################################################################
# FILENAME   : 01_generate_seeds.R
# AUTHOR     : Gulcin Ozer <ozer_gulcin@lilly.com>
# DATE       : 2025-08-11
# PROJECT    : LGM siRNA RNAseq Reports
# DESCRIPTION: Generate all possible 7mer seeds that can be utilized in reports
#              Generate G/U wobble version of seeds 
#################################################################################

library(Biostrings)


############################################################################
# Generate all possible 7mer seeds
# Generate ambiguity code to handle wobbles
# Reverse complement sequence to be searched


# Generate all combinations of 7 bases
bases = c("A", "U", "C", "G")
combinations = expand.grid(rep(list(bases), 7))

# Concatenate each row into a single string
seeds = apply(combinations, 1, paste0, collapse = "")

# Reverse complement of seeds will be searched in transcriptome
seq.searched = lapply(seeds, function(x) reverseComplement(DNAString(gsub("U", "T", x))))
seq.searched =  unlist(lapply(seq.searched, as.character))

# dataframe to be saved
dt.seed = data.frame(seed.7mer = seeds, seq.searched = seq.searched)


# Allow G/U wobble
# Transcript is already fixed
# To allow G on siRNA to be able to match to a U on mRNA, create ambiguity 
# Change G to G/A so that it could also match to U 
# IUPAC ambiguity code R for G/A
# To allow U on siRNA to be able to match to a G on mRNA, create ambiguity 
# Change U to U/C so that it could also match to G 
# IUPAC ambiguity code Y for C/U

# Example, G on siRNA match to U on mRNA:
# mRNA         5' CCAGGCCTCAGCTCCGGGC 3'
#                    |||||||
# seed            3' CCGGAGU 5'

# mRNA         5' CCAGGCCTCAGCTCCGGGC 3'
#                    ||||:||
# seed w/wobble   3' CCGGGGU 5'


# Example, U on siRNA match to G on mRNA: 
# mRNA         5' CCAGGCCTCAGCTCCGGGC 3'
#                    |||||||
# seed            3' CCGGAGU 5'

# mRNA         5' CCAGGCCTCAGCTCCGGGC 3'
#                    :||||||
# seed w/wobble   3' UCGGAGU 5'


seeds.wobble = gsub("G", "R", seeds)
seeds.wooble = gsub("U", "Y", seeds.wobble)
# reverse complement changes R to Y (C or U), changes Y to R (G or A)
seq.searched.wobble = lapply(seeds.wobble, function(x) reverseComplement(DNAString(RNAString(x))))
seq.searched.wobble =  unlist(lapply(seq.searched.wobble, as.character))

# dataframe to be saved
dt.seed.wobble = data.frame(seed.7mer.wobble = seeds.wobble, seq.searched.wobble = seq.searched)


# save to prelim data
save(dt.seed, file = "../data_prelim/Seeds_7mer.RData")
save(dt.seed.wobble, file = "../data_prelim/Seeds_7mer_w_wobble.RData")
