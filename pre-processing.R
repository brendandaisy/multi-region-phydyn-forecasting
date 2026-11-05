## This pre-processing script gets everything set up to make the
## BEAST2 XML (via BEAUti and a small amount of hand tweaking).

library(timtamslamR)
library(stringr)

## -------------------------------------------------------------------
## Distribute the sequenced samples across the dates that they
## correspond to so that there are no collisions as we are modelling
## these as unscheduled observations.

dated_seqs <- read_fasta("sample-sequences.fasta")

plot_dates(dated_seqs)

timed_seqs <- rename_dates_to_times_a(dated_seqs)

plot_times(timed_seqs)

write_fasta(timed_seqs, "sample-sequences-timed.fasta")

## -------------------------------------------------------------------
## Because we need to put the disaster times and sizes in the XML from
## BEAUti we compute these now and print them to standard output so
## they can be used later.

present <- get_present(dated_seqs, timed_seqs)
ts_dated <- read.csv("aggregated-occurrences.csv")
ts_timed <- rename_time_series(present, ts_dated)


cat("\n\nHere are the DISASTER TIMES\n\n")
cat(str_flatten(string = ts_timed$bwd_times, collapse = " "))
cat("\n\nHere are the DISASTER SIZES\n\n")
cat(str_flatten(string = ts_timed$size, collapse = " "))
cat("\n\nSaving \"present\" into last-sequence.Rds\n\n")
saveRDS(object = present, file = "last-sequence.Rds")
