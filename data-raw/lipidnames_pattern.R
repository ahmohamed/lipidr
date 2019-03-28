lipidnames_pattern <- list()
lipidnames_pattern$itsd_list <- c(
  "15:0-18:1(d7) PG", "15:0-18:1(d7) PI", "Sa1P 17:0",
  "So1P 17:1", "15:0-18:1(d7) PE", "18:1-d9 SM",
  "Cer d18:1/C12:0", "PS 33:1 d7", "15:0-18:1(d7) PC",
  "18:1(d7) Lyso PC", "18:1(d7) Lyso PE", "Cer1P d18:1/12:0",
  "GlucCer d18:1/12:0", "LacCer d18:1/12:0", "Sa 17:0", "So 17:1"
)
lipidnames_pattern$class <- "([[:alnum:]]{2,15})"
lipidnames_pattern$chain <- "(\\d{1,2}:\\d{1,2})"
lipidnames_pattern$isotope <- "(\\((d\\d)|(IS)\\))"
lipidnames_pattern$notes <- "(\\(.+\\))"
lipidnames_pattern$mol <- paste0(
  lipidnames_pattern$class, "[ -]",
  lipidnames_pattern$chain, "([/-_]", lipidnames_pattern$chain,
  ")?\\s*", "([/-_]", lipidnames_pattern$chain,
  ")?\\s*", lipidnames_pattern$notes, "?"
)
lipidnames_pattern$itsd <- paste0(
  "^", lipidnames_pattern$mol, ".*", lipidnames_pattern$isotope
)
# Examples: PE 32:2, PE 16:0/16:2, GlucCer 18:0/18:0
lipidnames_pattern$matching <- paste0(
  "^", lipidnames_pattern$mol,
  "(\\s*/\\s*", lipidnames_pattern$mol, ")?$"
)
usethis::use_data(lipidnames_pattern, overwrite = TRUE)
