################################################################################

# read file
file = "data/signal.tsv"
v = readdlm(file)

# setup hidden Markov model
hmm = setup(v)

################################################################################
