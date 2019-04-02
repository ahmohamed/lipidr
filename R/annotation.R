#' Parse molecule names to extract lipid class and chain information.
#'
#' Parse lipid names to return a data.frame containing lipid class,
#' chain length and unsaturation. Lipids should follow the pattern
#' 'class xx:x/yy:y', with class referring to the abbreviated lipid class,
#' xx:x as the composition of the first chain and yy:y as the second chain.
#' Alternatively, lipids can be supplied following the pattern 'class zz:z',
#' where zz:z indicates the combined chain length and unsaturation information.
#'
#' @param molecules A character vector containing lipid molecule names.
#' @return A data.frame with lipid annotations as columns. Input lipid names
#'   are given in a column named "Molecule".
#'
#' @importFrom dplyr %>% filter left_join full_join
#' @export
#'
#' @examples
#' lipid_list <- c(
#'   "Lyso PE 18:1(d7)",
#'   "PE(32:0)",
#'   "Cer(d18:0/C22:0)",
#'   "TG(16:0/18:1/18:1)"
#' )
#' annotate_lipids(lipid_list)
annotate_lipids <- function(molecules) {
  .data_internal("lipidDefaults")
  def <- .myDataEnv$lipidDefaults$clean_mols
  not_in_db <- molecules[!molecules %in% def$Molecule]

  if (length(not_in_db) == 0) {
    def %>%
      filter(Molecule %in% molecules) %>%
      return()
  }

  clean_ <- .clean_molecule_name(not_in_db)
  if (any(clean_$not_matched)) {
    warning(
      "Some lipid names couldn't be parsed because they don't follow ",
      "the pattern 'CLS xx:x/yy:y' \n    ",
      clean_$Molecule[clean_$not_matched]
    )
  }

  clean_ %>%
    filter(!not_matched) %>%
    .parse_lipid_info() %>%
    .left_join_silent(.myDataEnv$lipidDefaults$class_info) %>%
    .full_join_silent(def %>% filter(Molecule %in% molecules)) %>%
    return()
}

#### Internal functions ############################################

.clean_molecule_name <- function(lipids_list) {
  .data_internal("lipidnames_pattern")
  p <- .myDataEnv$lipidnames_pattern
  olipids <- trimws(lipids_list)

  # PC(O-32:0) --> PCO-32:0
  p2 <- paste0(p$class, "[ -]*\\(", "([OP]-\\d{1,2}:\\d{1,2}[^)]*)", "\\)")
  olipids <- gsub(p2, "\\1\\2", olipids)

  # Cer(d18:0/C18:0) --> Cer d18:0/18:0
  p2 <- paste0(
    p$class, "[ -]*\\(",
    "([[:alpha:]]{0,1}\\d{1,2}:\\d{1,2}[^)]*)", "\\)"
  )
  olipids <- gsub(p2, "\\1 \\2", olipids)

  # Cer d18:0/C18:0 --> Cer 18:0/18:0
  chain_p2 <- "([[:alpha:]]{0,1})(\\d{1,2}:\\d{1,2})"
  olipids <- gsub(chain_p2, "\\2", olipids)

  # TG 14:1 18:1 18:1 --> TG 14:1/18:1/18:1
  p2 <- paste0(p$chain, "[-_ ]", p$chain)
  olipids <- sub(p2, "\\1/\\2", olipids)
  olipids <- sub(p2, "\\1/\\2", olipids)

  # TG 14:1 18:1 18:1 --> TG 14:1/18:1/18:1
  p2 <- paste0(p$class, "[ -]", p$chain)
  olipids <- gsub(p2, "\\1 \\2", olipids)

  # Lyso PC --> LPC
  olipids <- sub("Lyso\\s?P", "LP", olipids)

  # 18:1(d7) LPC --> LPC 18:1(d7)
  olipids <- sub(
    "(^[[:digit:]+][^[:blank:]]*)[[:blank:]?](.*$)",
    "\\2 \\1", olipids
  )

  # LPC 18:1-d7 --> LPC 18:1(d7)
  olipids <- sub("[^(](d\\d)(\\D|$)", "(\\1)\\2", olipids)
  # Trim
  olipids <- sub(" NEG$", "", olipids)
  olipids <- sub(" ID\\d+$", "", olipids)

  is_itsd <- grepl(p$itsd, olipids) |
    lipids_list %in% p$itsd_list |
    olipids %in% p$itsd_list

  return(data.frame(
    Molecule = lipids_list, clean_name = olipids,
    ambig = grepl(paste0("^", p$mol, "(\\s*/\\s*", p$mol, ")$"), olipids),
    not_matched = (!grepl(p$matching, olipids) & !grepl(p$itsd, olipids)),
    itsd = is_itsd
  ))
}

#' @importFrom dplyr %>% mutate rowwise ungroup
#' @importFrom tidyr separate
.parse_lipid_info <- function(clean_df) {
  .data_internal("lipidnames_pattern")
  p <- .myDataEnv$lipidnames_pattern

  clean_df %>%
    mutate(
      first_mol = sub(
        paste0("^(", p$mol, ")(\\s*/\\s*", p$mol, ")?$"),
        "\\1", clean_name
      ),
      first_mol = sub(
        paste0(
          p$class, "[ -]",
          p$chain, "([/-]", p$chain, ")?",
          "([/-]", p$chain, ")?.*$"
        ),
        "\\1#$#\\2#$#\\4#$#\\6", first_mol
      )
    ) %>%
    separate(
      first_mol, c("class_stub", "chain1", "chain2", "chain3"),
      sep = "#\\$#"
    ) %>%
    separate(
      chain1, c("l_1", "s_1"),
      sep = "\\:", remove = FALSE, convert = TRUE
    ) %>%
    separate(
      chain2, c("l_2", "s_2"),
      sep = "\\:", remove = FALSE, convert = TRUE, fill = "right"
    ) %>%
    separate(
      chain3, c("l_3", "s_3"),
      sep = "\\:", remove = FALSE, convert = TRUE, fill = "right"
    ) %>%
    rowwise() %>%
    mutate(
      total_cl = sum(l_1, l_2, l_3, na.rm = TRUE),
      total_cs = sum(s_1, s_2, s_3, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Defined as lipid annotations
utils::globalVariables(c(
  "first_mol",
  "chain1", "chain2", "chain3",
  "l_1", "s_1", "l_2", "s_2", "l_3", "s_3",
  "total_cl", "total_cs",
  "Molecule", "clean_name", "ambig", "not_matched", "itsd",
  "Class"
  ))
