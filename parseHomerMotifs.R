#! /usr/bin/env Rscript

parseHomerMotif <- function(motif) {
  res <- read_lines(motif)
  
  head <- res[1]
  head <- read.table(text=head, header = FALSE, sep = "\t")
  colnames(head) <- c("motif","details","logOddsThreshold","logPvalue","gap","info")
  head <- head %>% 
  separate(details, into = c("consensus","details"), sep = ",") %>%
  separate(info, into = c("target","background"), sep = ",") %>%
  separate(target, into = c("counts_target", "frequency_target"), sep = "\\(") %>%
  separate(background, into = c("counts_background", "frequency_background"), sep = "\\(") %>%
  mutate(motif = gsub(">", "", motif),
         counts_target = as.numeric(gsub("T:", "", counts_target)), 
         counts_background = as.numeric(gsub("B:", "", counts_background)),
        frequency_target = as.numeric(gsub("\\%\\)" ,"", frequency_target))/100,
        frequency_background = as.numeric(gsub("\\%\\)" ,"", frequency_background))/100,
        gap = NULL,
        #logo = sprintf('<img src="%s/%s" width="104"></img>', motif_dir, "motif2.logo.svg")
        ) %>%
  mutate(log2Ratio = log2(frequency_target/frequency_background))
  
  pwm <- res[-1]
  pwm <- t(as.matrix(read.table(text=pwm, sep ="\t", header = FALSE)))
  rownames(pwm) <- c("A", "C", "G", "T")
  
  return(list(details = head, pwm = pwm))

}
