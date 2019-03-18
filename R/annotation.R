#' Parse Molecule names to extract class and chain information.
#'
#' @param data Skyline data.frame created by \code{\link{read_skyline}}
#' @return
#' 
#' @importFrom dplyr %>% filter left_join full_join
#' @export
#'
#' @examples
annotate_lipids <- function(data){
  def = .myDataEnv$lipidDefaults$clean_mols
  mols = unique(data$Molecule)
  not_in_db = mols[!mols %in% def$Molecule]
  
  if(length(not_in_db) == 0) {
    def %>% filter(Molecule %in% mols) %>% 
      return ()
  }
  
  clean_ = .clean_molecule_name(not_in_db)
  if(any(clean_$not_matched)) {
    warning("Some lipid names couldn't be parsed because they don't follow the pattern 'CLS xx:x/yy:y' \n    ",
            clean_$Molecule[clean_$not_matched]
    )
  }
  
  clean_ %>% filter(!not_matched) %>% .parse_lipid_info() %>%
    .left_join.silent(.myDataEnv$lipidDefaults$class_info) %>%
    .full_join.silent(def %>% filter(Molecule %in% mols)) %>%
    return()
}

#############################################################################################################

.en <- function(...) paste0("(", ..., ")")
csep = "[/-_]"

itsd_list = c(
  '15:0-18:1(d7) PG','15:0-18:1(d7) PI','Sa1P 17:0','So1P 17:1','15:0-18:1(d7) PE',
  '18:1-d9 SM','Cer d18:1/C12:0','PS 33:1 d7','15:0-18:1(d7) PC','18:1(d7) Lyso PC',
  '18:1(d7) Lyso PE','Cer1P d18:1/12:0','GlucCer d18:1/12:0','LacCer d18:1/12:0','Sa 17:0',
  'So 17:1')
class_p = "[[:alnum:],-]{2,35}"
dbond_config = "\\((?:\\d{1,2}[ZE][,]*)+\\)"
chain_notes = "\\((?:\\d{1,2}[^)]*)+\\)"
chain_p = paste0("[dth]?\\d{1,2}:\\d{1,2}",.en("?:", chain_notes), "*")
chain_multi_p = paste0(.en("?:", chain_p, csep, "*"), "+")
isotope_p = "\\(IS\\)|\\(d\\d\\)"
notes_p = "\\(.+\\)\\W*"
notes_multi = paste0(.en("?:", notes_p), "+")
mol_p = paste0("^", class_p, "[ -]*\\(", chain_multi_p, "\\)", .en("?:", notes_multi), "*")
mol_multi = paste0(mol_p, "\\s*/\\s*", mol_p)

# PE 32:2, PE 16:0/16:2, GlucCer 18:0/18:0
p = paste0("^", mol_p, "(\\s*/\\s*",mol_p,")?$") #"^([[:alnum:]]{2,7}) (\\d{2}:\\d{1,2})(/\\d{2}:\\d{1,2})?$"
itsd = paste0("^", mol_p, ".*", isotope_p)#"^([[:alnum:]]{2,7}) (\\d{2}:\\d{1,2})(/\\d{2}:\\d{1,2})?(\\(d\\d\\))$"

.preprocess_molecule_names <- function(lipids_list) {
  olipids = trimws(lipids_list)
  
  # PC(O-32:0) --> PCO-32:0
  # PC(O-16:0_18:1)/PC(P-16:2_18:1) --> PCO-16:0_18:1/PCP-16:2_18:1
  # PC(O-16:0/18:2(9Z))/PC(P-16:2(9Z,12E)/18:1(12Z)) --> PCO-16:0_18:1/PCP-16:2_18:1
  p2 = paste0(.en(class_p), "[ -]*\\(", "([OP])-",.en(chain_multi_p), "\\)")
  olipids = gsub(p2, "\\1\\2(\\3)", olipids)
  olipids
}


.parse_molecule <- function(lipids_list) {
  mol_p_en = paste0(.en(class_p), "[ -]*\\(", .en(chain_multi_p), "\\)",.en(notes_multi), "*")
  matches = str_match(lipids_list, mol_p_en)
  colnames(matches) = c("match", "class", "chains", "notes")
  data.frame(input=lipids_list, matches)
}


.clean_molecule_name <- function(lipids_list) {
  olipids = trimws(lipids_list)
  
  # PC(O-32:0) --> PCO-32:0
  # PC(O-16:0_18:1)/PC(P-16:2_18:1) --> PCO-16:0_18:1/PCP-16:2_18:1
  # PC(O-16:0/18:2(9Z))/PC(P-16:2(9Z,12E)/18:1(12Z)) --> PCO-16:0_18:1/PCP-16:2_18:1
  p2 = paste0(.en(class_p), "[ -]*\\(", .en("[OP]-", chain_multi_p), "\\)")
  olipids = gsub(p2, "\\1\\2", olipids)
  
  # Cer(d18:0/C18:0) --> Cer d18:0/18:0
  p2 = paste0(.en(class_p), "[ -]*\\(", chain_multi_p, "\\)")
  olipids = gsub(p2, "\\1 \\2", olipids)
  
  # Cer d18:0/C18:0 --> Cer 18:0/18:0
  chain_p2 = "([[:alpha:]]{0,1})(\\d{1,2}:\\d{1,2})"
  olipids = gsub(chain_p2, "\\2", olipids)
  
  #TG 14:1 18:1 18:1 --> TG 14:1/18:1/18:1
  p2 = paste0(chain_p, "[-_ ]", chain_p)
  olipids = sub(p2, "\\1/\\2", olipids)
  olipids = sub(p2, "\\1/\\2", olipids)
  
  #TG 14:1 18:1 18:1 --> TG 14:1/18:1/18:1
  p2 = paste0(class_p, "[ -]", chain_p)
  olipids = gsub(p2, "\\1 \\2", olipids)
  
  # Lyso PC --> LPC
  olipids = sub("Lyso\\s?P", "LP", olipids)
  
  # 18:1(d7) LPC --> LPC 18:1(d7)
  p2 = paste0("^", chain_p, "(\\s*.*)\\W", class_p, "(.*$)")
  olipids = sub(p2, "\\3 \\1\\2\\4",  olipids)
  
  # LPC 18:1-d7 --> LPC 18:1(d7)
  olipids = sub("[^(](d\\d)(\\D|$)", "(\\1)\\2",  olipids)
  # Trim
  olipids = sub(" NEG$", "", olipids)
  olipids = sub(" ID\\d+$", "", olipids)
  
  #sort(olipids[!grepl(p, olipids) & !grepl(itsd, olipids)])
  
  
  return (data.frame(
    Molecule=lipids_list, clean_name=olipids, 
    ambig=grepl(paste0("^", mol_p, "(\\s*/\\s*",mol_p,")$"), olipids),
    not_matched=(!grepl(p, olipids) & !grepl(itsd, olipids)),
    itsd=grepl(itsd, olipids) | lipids_list %in% itsd_list | olipids %in% itsd_list
  ))
}

#' @importFrom dplyr %>% mutate rowwise ungroup
#' @importFrom tidyr separate
.parse_lipid_info <- function(clean_df) {
  clean_df %>% mutate(
    first_mol = sub(paste0("^(", mol_p, ")(\\s*/\\s*",mol_p,")?$"), "\\1", clean_name),
    first_mol = sub(paste0(class_p, "[ -]", 
                           chain_p, "([/-]", chain_p, ")?", 
                           "([/-]", chain_p, ")?.*$"), 
                    "\\1#$#\\2#$#\\4#$#\\6", first_mol)
  ) %>% 
    separate(first_mol, c("class_stub", "chain1", "chain2", "chain3"), sep="#\\$#") %>%
    separate(chain1, c("l_1", "s_1"), sep="\\:", remove = F, convert=T) %>% 
    separate(chain2, c("l_2", "s_2"), sep="\\:", remove = F, convert=T, fill="right") %>% 
    separate(chain3, c("l_3", "s_3"), sep="\\:", remove = F, convert=T, fill="right") %>% 
    rowwise() %>%
    mutate(
      total_cl = sum(l_1, l_2, l_3, na.rm=T),
      total_cs = sum(s_1, s_2, s_3, na.rm=T)
    ) %>% ungroup
}
