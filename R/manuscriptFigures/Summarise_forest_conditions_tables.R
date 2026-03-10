library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(readr)

# This script reads in corrected files digitized from historic forest health conditions reports from the USDA forest service for R9
# Script overview:
#     -set up states, pest names, etc and helper functions
#     -
`%||%` <- function(a,b) if (!is.null(a) && length(a)>0 && !all(is.na(a))) a else b

##########################################################################
# 1) set up helper functions to standardize state names for the region----

new_england_states <- c("Connecticut","Maine","Massachusetts","New Hampshire","Rhode Island","Vermont")

mid_atlantic_states <- c("New York","New Jersey","Pennsylvania","Delaware","Maryland")

lake_states <- c("Michigan","Minnesota","Wisconsin")

# midwest states not typically called "lake states"
midw_states <- c("Illinois","Indiana","Ohio", "Iowa", "Missouri")



# separate for northeastern states flags:
northeast_states <- c(
  "Connecticut","Delaware","Maine","Maryland","Massachusetts",
  "New Hampshire","New Jersey","New York","Pennsylvania",
  "Rhode Island","Vermont","West Virginia","Ohio"
)

# full "region-wide" set of states:
region_states <- unique(c(new_england_states, mid_atlantic_states, northeast_states, midw_states))


# define aspen range, as one record has this:
aspen_states <- region_states

abbr_to_name <- c(
  "CT"="Connecticut","DE"="Delaware","ME"="Maine","MD"="Maryland","MA"="Massachusetts",
  "NH"="New Hampshire","NJ"="New Jersey","NY"="New York","PA"="Pennsylvania","RI"="Rhode Island",
  "VT"="Vermont","WV"="West Virginia","OH"="Ohio",
  "MI"="Michigan","MN"="Minnesota","WI"="Wisconsin","IL"="Illinois","IN"="Indiana"
)

# helper functions to normalize state names
normalize_state_token <- function(x) {
  x2 <- str_to_upper(str_squish(x))
  ifelse(x2 %in% names(abbr_to_name), unname(abbr_to_name[x2]),
         str_to_title(str_to_lower(str_squish(x))))
}

# function to remove blank spaces and standardize the separators in the location column
std_sep <- function(x) {
  x %>%
    tidyr::replace_na("") %>%
    str_replace_all("[;|/|.|:]", ",") %>%
    str_replace_all("\\s*,\\s*", ",") %>%
    str_replace_all("\\bAnd\\b|\\band\\b", ",") %>%
    str_squish()
}

# split_tokens helper
split_tokens <- function(x) {
  toks <- str_split(std_sep(x), ",", simplify = FALSE)[[1]]
  toks <- str_squish(toks)
  toks[toks != ""]
}

# the `expand_location()` function is the workhorse function
# to separate a single row entry into multi-row entries for each pest reflecting one state per row
expand_location <- function(loc) {
  toks <- split_tokens(loc)
  low  <- str_to_lower(toks)

  has_regionwide <- any(str_detect(low, "\\b(region[- ]?wide|area[- ]?wide|throughout (the )?region)\\b"))
  has_newengland <- any(str_detect(low, "\\bnew england\\b"))
  has_midatl     <- any(str_detect(low, "\\b(mid[- ]?atlantic|middle atlantic)\\b"))
  has_lakest     <- any(str_detect(low, "\\b(lake states|great lakes)\\b"))
  has_northeast     <- any(str_detect(low, "\\b(northeastern region|northeast|north east|(all )?northeastern (united )?states|(throughout )?northeastern states)\\b"))
  has_aspen_rng     <- any(str_detect(low, "\\b(throughout aspen range)\\b"))

  special_pat <- "\\b(region[- ]?wide|area[- ]?wide|new england|mid[- ]?atlantic|middle atlantic|lake states|great lakes|throughout (the )?region|northeastern region|northeast|north east|(all )?northeastern (united )?states|(throughout )?northeastern states|throughout aspen range)\\b"

  toks_clean <- toks[!str_detect(low, special_pat)]

  explicit_states <- map_chr(toks_clean, normalize_state_token)
  explicit_states <- explicit_states[explicit_states != "" & !is.na(explicit_states)]

  # if it is a special case, deal with it here
  unique(c(
    explicit_states,
    if (has_newengland) new_england_states else character(0),
    if (has_midatl)     mid_atlantic_states else character(0),
    if (has_lakest)     lake_states else character(0),
    if (has_northeast)  northeast_states else character(0),
    if (has_aspen_rng)  aspen_states else character(0),
    if (has_regionwide) region_states else character(0)

  ))
}

# state.name and state.abb 
all_state_names <- state.name
all_state_abbr  <- state.abb
abbr_to_state <- setNames(all_state_names, all_state_abbr)

# regex for matching full state names (longest first avoids partial matches)
state_name_re <- paste0("\\b(", paste(sort(all_state_names, decreasing = TRUE), collapse="|"), ")\\b")
state_abbr_re <- paste0("\\b(", paste(all_state_abbr, collapse="|"), ")\\b")

# converting abbreviations to full state names
token_to_state <- function(tok) {
  tok <- str_squish(tok)
  up  <- str_to_upper(tok)
  if (up %in% names(abbr_to_state)) return(unname(abbr_to_state[up]))
  # otherwise assume it's a full name
  str_to_title(str_to_lower(tok))
}

# clean up the state strings
clean_state_string <- function(x) {
  x %>%
    replace_na("") %>%
    str_replace_all("^\\s*and\\s+", "") %>%   # remove leading "and"
    str_replace_all("\\.", " ") %>%           # remove periods
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

# normalizing the state field
normalize_state_field <- function(state_value) {

  s <- clean_state_string(state_value)
  s_low <- str_to_lower(s)

  # ---- Expand region-wide phrases ----
  if (str_detect(s_low, "northeastern united states|all northeastern states|throughout northeastern states")) {
    return(northeast_states)
  }

  # ---- Remove broken fragments ----
  if (s %in% c("New", "")) {
    return(character(0))
  }

  # ---- Split if multiple states jammed together ----
  # e.g., "Minnesota Wisconsin"
  parts <- str_split(s, "\\s{2,}|,|;", simplify = FALSE)[[1]]
  parts <- str_squish(parts)

  # Keep only real NE states
  #parts <- parts[parts %in% northeast_states]

  unique(parts)
}

# ----------------------------------------------------------------
# 2) Function for Host standardization and expanding to long-format
#    - keep a full host species name if present (e.g., "balsam fir")
#    - split up multi-species entries into separate rows
#    - use to genus level groups ("fir"/"pine") or types ("hardwood"/"conifer") if these are the only values present

# set up a table with different possible patterns & standardized
#Note: these will need to be updated if additional species show up as we add condition reports
species_patterns <- tibble::tribble(
  ~pattern, ~host_std,

  # conifers
  "\\bbalsam fir\\b", "balsam fir",
  "\\bwhite spruce\\b", "white spruce",
  "\\bblack spruce\\b", "black spruce",
  "\\bnorway spruce\\b", "norway spruce",
  "\\bred spruce\\b", "red spruce",
  "\\bblue spruce\\b", "blue spruce",
  "\\bwhite cedar\\b|\\bNorthern white cedar\\b","white cedar",
  "\\bred cedar\\b|\\bNorthern red cedar\\b|\\bredcedar\\b","red cedar",
 
  "\\beastern hemlock\\b|\\bhemlock\\b", "eastern hemlock",

  # assign tamarack as eastern larch
 # "\\btamarack\\b", "tamarack",
  "\\beastern larch\\b|\\larche?\\b|\\btamarack\\b", "eastern larch",
  "\\beuropean lar(ch|che)\\b|\\beuropean larch\\b", "european larch",

  # Pines (keep distinct). Putting white pine as eastern white pine
  "\\beastern white pine\\b|\\bwhite pine\\b", "eastern white pine",
  "\\bjack pine\\b", "jack pine",
  "\\bred pine\\b", "red pine",
  "\\bscotch pine\\b|\\bscot's pine\\b|\\bscots pine\\b", "scots pine",
  "\\baustrian pine\\b", "austrian pine",
  "\\bvirginia pine\\b|\\bvirgina pine\\b", "virginia pine",
  "\\bshortleaf pine\\b", "shortleaf pine",
  "\\bloblolly pine\\b|\\bloblolly\\b", "loblolly pine",

  # All hardwoods:

  "\\bamerican beech\\b", "american beech",

 # Birches:
  "\\byellow birch\\b", "yellow birch",
  "\\bpaper birch\\b|\\bwhite birch\\b", "paper birch",

  "\\bgrey birch\\b|\\bgray birch\\b", "grey birch",

 # Oaks:
 "\\bpin oak\\b", "pin oak",
 "\\bwhite oak\\b", "white oak",
 "\\bblack oak\\b", "black oak",
 "\\bchestnut oak\\b", "chestnut oak",

 "\\bred oaks?\\b|\\bred oak group\\b", "red oak",




 # ashes

  "\\bwhite ash\\b", "white ash",
  "\\bgreen ash\\b", "green ash",
  "\\bblack ash\\b", "black ash",

 # maples:

 "\\bsugar maples?\\b", "sugar maple",
 "\\bred maples?\\b", "red maple",
 "\\bsilver maples?\\b", "silver maple",


 "\\balder\\b", "alder",
 "\\bflowering dogwood\\b|\\bdogwood\\b", "dogwood",

 "\\bblack locust\\b|\\blocust\\b", "locust",
  "\\bass?pens?\\b", "aspen",
  "\\bbasswood\\b", "basswood",
  "\\bsycamore\\b", "sycamore",
  "\\bhickory\\b|\\bhickories\\b", "hickory",
  "\\bblack walnut\\b|\\bwalnuts?\\b|\\bwainut\\b", "black walnut",
  "\\bbutternuts?\\b", "butternut",

  "\\bblack cherry\\b", "black cherry",

  "\\bapples?\\b", "apples",

  "\\bpoplars?\\b", "poplar",
  "\\byellow poplars?\\b|\\btuliptrees?\\b", "yellow poplar",
  "\\bsassafras\\b", "sassafras",
  "\\bhackberry\\b", "hackberry",
  "\\bRussian olive\\b", "russian olive"
) %>%
  mutate(pattern = regex(pattern, ignore_case = TRUE))

# set up some string patters for a genus level/general strings
generic_patterns <- tibble::tribble(
  ~pattern, ~host_std,
  "\\bfir\\b", "fir",
  "\\bbirch\\b|\\bbirches\\b", "birch",
  "\\bash\\b|\\bbashes\\b", "ash",
  "\\bmaples?\\b|\\bboxelder\\b", "maple",
  "\\bulmus\\b|\\belm\\b", "elm",
  "\\bspruce\\b", "spruce",
  "\\boaks?\\b|\\boak\\b", "oak",
  "\\byews?\\b", "yew",
  "\\bpine\\b|\\bpines\\b", "pine",
  "\\cedars?\\b","cedar",
  "\\bcherry\\b|\\bcherries\\b", "cherry",
  "\\bconifers?\\b", "conifers",
  "\\b(other )?hardwoods?\\b|\\b(various )?hardwoods?\\b", "hardwoods",
  "\\ball trees\\b", "all trees",
  "\\bvarious species\\b", "various species"
) %>%
  mutate(pattern = regex(pattern, ignore_case = TRUE))

# function to clean and standardize the tree species host names
# removing some extra words and standardizing separators
clean_host_string <- function(x) {
  x %>%
    replace_na("") %>%
    str_to_lower() %>%
    str_replace_all("&", " and ") %>%
    str_replace_all("particularly", " ") %>%
    str_replace_all("mostly", " ") %>%
    str_replace_all("hardwoods:", "hardwoods") %>%
    str_replace_all("[()\\[\\]]", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

# further standardize the host names
tokenize_hosts <- function(x) {
  x <- clean_host_string(x)
  x <- str_replace_all(x, "\\band\\b", ",")
  x <- str_replace_all(x, "[;/]", ",")
  x <- str_replace_all("\\s*,\\s*", ",", x)
  toks <- unlist(str_split(x, ",", simplify = FALSE))
  toks <- str_squish(toks)
  toks <- toks[toks != ""]
  toks <- toks[!toks %in% c("ee", "pennsylvania")]  
  toks
}

match_all <- function(s, patterns_tbl) {
  hits <- patterns_tbl$host_std[str_detect(s, patterns_tbl$pattern)]
  unique(hits)
}

# function to split standardized host strings into a single row for each host_std



standardize_host_row <- function(host_raw) {
  if (is.na(host_raw) || stringr::str_squish(host_raw) == "") return(NA_character_)
  # get a token for each separate host
  toks <- tokenize_hosts(host_raw)
  
  # for each host token see if it is specific or general
  out <- unique(unlist(purrr::map(toks, function(tok) {
    tok_clean <- clean_host_string(tok)
    
    # 1) Try specifics for this token
    spec <- unique(match_all(tok_clean, species_patterns))
    if (length(spec) > 0) return(spec)
    
    # 2) If there are not any matches to species, try generic_patterns table for this token
    gen <- unique(match_all(tok_clean, generic_patterns))
    if (length(gen) > 0) return(gen)
    
    character(0)
  })))
  
  out <- out[!is.na(out) & out != ""]
  if (length(out) == 0) return(NA_character_)
  out
}

# ----------------------------------------------------------
# 3) functions to clean up forest pest agents----

# cleaning agen names, removing anything in parentheses, and replacing some strange character strings
clean_agent_name <- function(x) {
  x %>%
    replace_na("") %>%
    str_to_lower() %>%
    str_replace_all("\\(.*?\\)", "") %>%  # remove scientific names
    str_replace_all("[^a-z0-9\\s-]", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}


# create pest/agent dictionary to correct different pests:
# Note: this will need to be updated with additional new pests and typos that arise
pest_dict <- tibble::tribble(
  ~pattern, ~Agent_std,

  # Defoliators / insects
  "spruce budworm", "Spruce Budworm",
  "bruce spanworm", "Bruce Spanworm",
  "forest tent caterpillar", "Forest Tent Caterpillar",
  "Eastern Caterpillar Tent", "Eastern Tent Caterpillar",
  "Eastern Tent Caterpillar", "Eastern Tent Caterpillar",
  "gypsy moth", "Spongy Moth",
  "fall webworm", "Fall Webworm",
  "fall cankerworm", "Fall Cankerworm",
  "saddled prominent", "Saddled Prominent",
  "linden looper", "Linden Looper",
  "half wing geometer|half-wing geometer", "Half-wing Geometer",
  "cherry scallop shell moth|Cherry moth scallop shell moth", "Cherry Scallop Shell Moth",
  "basswood thrips", "Basswood Thrips",
  "jack pine budworm", "Jack Pine Budworm",
  "Spearmarked Moth Black", "Spearmarked Black Moth",
  "Balsam Shootboring Sawlly", "Balsam Shootboring Sawfly",
  "Balsam Adelges piceae woolly adelgid", "Balsam Woolly Adelgid",
  "Birch casebearer, Coleophora serratella", "Birch Casebearer",
  "Birch leafminer|Birch leaf miner", "Birch Leaf Miner",
  "Marssonina Lef Spot", "Marssonina Leaf Spot",
  "Nantucket Moth Pine Tip", "Nantucket Pine Tip Moth",
  "Pine engraver, ips pini", "Pine engraver",
  "Spring Paleacrita cankerworm vernata", "Spring Cankerworm",
  "Maple leafcutter Paraclemensia acerifoliella", "Maple leafcutter",
  "Twolined Chestnut Borer", "Twolined Chestnut Borer",

  # Diseases / pathogens (incl. variants)
  "beech bark disease.*var\\s*faginata|var\\s*faginata.*beech bark disease|beech bark disease var faginata",
  "Beech Bark Disease",
  "beech bark disease", "Beech Bark Disease",
  "Beech Disease Bark Disease", "Beech Bark Disease",
  "Pinewood nematode|pinewood nematode", "Pine Wood Nematode",
  "Saratoga spittlebug|Pine spittlebug|Spittlebugs", "Spittlebugs",
  "Stillwell's syndrome|Stillwell S Syndrome", "Stillwells Syndrome",
  "Wainut Caterpillar", "Walnut Caterpillar",
  "Walking stick|Walkingsticks", "Walkingsticks",


  "dogwood anthracnose", "Dogwood Anthracnose",
  "^anthracnose$|\\banthracnose\\b", "Anthracnose",

  "oak wilt", "Oak Wilt",
  "dutch elm disease", "Dutch Elm Disease",
  "white pine blister rust", "White Pine Blister Rust",
  "annosus root rot|Annosus root and butt rot", "Annosus Root and Butt Rot",
  "shoestring root rot", "Shoestring Root rot",
  "dwarf mistletoe|eastern dwarf mistletoe", "Dwarf mistletoe",
  "scleroderris canker", "Scleroderris Canker",
  "diplodia|sphaeropsis", "Diplodia Tip Blight",
  "larch needlecast", "Larch Needlecast",
  "european larch canker", "European Larch Canker",
  "pine needle rust", "Pine Needle Rust",

  # phys.needle droop
  "physiological droop needle", "Physiological Needle Droop",

  # Rust oddity you have
  "pine oak.*gall rust|pine-oak gall rust|pine oak gall rust", "Pine-oak Gall Rust",

  # Decline complexes (collapse all variants)
  "ash.*(decline|dieback|mortality)", "Ash decline complex",
  "birch.*(decline|dieback|mortality)", "Birch decline complex",
  "maple.*(decline|dieback|mortality)", "Maple decline complex",
  "oak.*(decline|dieback|mortality)", "Oak decline complex",
  "spruce.*(decline|mortality)", "Spruce decline complex",
  "larch.*(decline|mortality)", "Larch decline complex",
  "hardwood decline", "Hardwood decline complex",

  # Abiotic
  "drought","Drought",
  "frost|Frost damage", "Frost Injury",
  "hail", "Hail damage",
  "snow|Snow and winter injury|winter injury", "Winter injury",
  "tornado|tornadoes", "Tornado",
  "Weather|Weather-related injury", "Weather",

  "ozone|air pollution|oxidant", "Pollution"
) %>%
  mutate(pattern = regex(pattern, ignore_case = TRUE))

# clean up names from common OCR errors:
clean_agent_name <- function(x) {
  x %>%
    replace_na("") %>%
    str_to_lower() %>%
    str_replace_all("\\(.*?\\)", "") %>%                 # remove parenthetical sci names
    str_replace_all("[^a-z0-9\\s-]", " ") %>%            # drop punctuation
    str_replace_all("\\s+", " ") %>%
    str_squish() %>%
    # ---- common OCR/typo fixes you have ----
  str_replace_all("\\bwortality\\b", "mortality") %>%
    str_replace_all("\\bcarkerwoun\\b|\\bcankerworn\\b|\\bcarkerworm\\b", "cankerworm") %>%
    str_replace_all("\\bburopean\\b", "european") %>%
    str_replace_all("\\boperophters\\b", "operophtera") %>%
    str_replace_all("\\bfumi\\s*ferena\\b|\\bfumiferena\\b", "fumiferana") %>%
    str_replace_all("\\bmalecosoma\\b|\\bmalacosom\\b|\\bmalacosona\\b", "malacosoma") %>%
    str_replace_all("\\bgremeniella\\b|\\bgrenmeniella\\b|\\bgremmeniella\\b", "gremmeniella")
}

# function to apply dictionary and agent standarization
standardize_agent <- function(agent_raw) {

  cleaned <- clean_agent_name(agent_raw)

  # First try dictionary matches
  for (i in seq_len(nrow(pest_dict))) {
    if (str_detect(cleaned, pest_dict$pattern[i])) {
      return(pest_dict$Agent_std[i])
    }
  }

  # Fallback: title-case cleaned text
  str_to_title(cleaned)
}

# ------------------------------------------------------------------------------
# 4) REMARKS: pull out state-specific remarks & get some numeric values
#     CAUTION: this is still clunky and likely results in some dropped or misassigned state remarks

# extract number gets any numbers listed
extract_number <- function(text, pattern) {
  t <- str_to_lower(text %||% "")
  m <- str_match(t, pattern)
  ifelse(is.na(m[,2]), NA_real_, as.numeric(str_remove_all(m[,2], ",")))
}
# use extract number to get any numbers associated with "acres"
extract_acres_any <- function(text) {
  # if million acres: 
  has.million <- str_match(text, "million")
  if(!is.na(has.million)){
    acres <- extract_number(text, "([-+]?\\d*\\.?\\d+)\\s*million( acres|acre|ac\\b)")*10^6
  }else{
    acres <- extract_number(text, "([0-9]{1,3}(?:,[0-9]{3})+|[0-9]+)\\s*(acres|acre|ac\\b)")
  }
  return(acres)
}

# get any acres listed for defoliated
extract_acres_defoliated <- function(text) {
  # "12,000 acres defoliated" or "defoliation of 12,000 acres"
  extract_number(text, "(?:defoliat\\w*\\s*(?:of|over)?\\s*)?([0-9]{1,3}(?:,[0-9]{3})+|[0-9]+)\\s*acres\\s*(?:defoliat\\w*)?")
}

# get any percentages
extract_percent <- function(text) {
  extract_number(text, "([0-9]+(?:\\.[0-9]+)?)\\s*%")
}

# get trees dead
extract_dead_trees <- function(text) {
  extract_number(text, "([0-9]{1,3}(?:,[0-9]{3})+|[0-9]+)\\s*(dead\\s+trees?|trees?\\s+dead|trees?\\s+killed|killed?\\s+trees)")
}

extract_trap_count <- function(text) {
  # "X moths per trap", "X/trap", "X caught in traps", etc. (very flexible)
  extract_number(text, "([0-9]{1,3}(?:,[0-9]{3})+|[0-9]+)\\s*(?:per\\s+trap|/\\s*trap|traps?|caught)")
}

extract_egg_masses <- function(text) {
  extract_number(text, "([0-9]{1,3}(?:,[0-9]{3})+|[0-9]+)\\s*(egg\\s*masses?|eggmasses?)")
}

# flags to identify if documented trends: 
flag_increasing <- function(txt) {
  t <- str_to_lower(txt %||% "")
  str_detect(t, "\\bincreas(ing|ed|e)\\b|\\bexpan(d|ding|sion)\\b|\\boutbreak\\b|\\bepidemic\\b|\\bspread\\b|\\bhigh(er|est)?\\b|\\bpeak\\b")
}

flag_decreasing <- function(txt) {
  t <- str_to_lower(txt %||% "")
  str_detect(t, "\\bdecreas(ing|ed|e)\\b|\\bsubsided\\b|\\bdown\\b|\\bdeclin(ed|ing|e)\\b|\\blow(er)?\\b|\\bless\\b")
}
flag_stable <- function(txt) {
  t <- str_to_lower(txt %||% "")
  str_detect(t, "\\bstable\\b|\\blittle change\\b|\\bunchanged\\b|\\bremain(ed)\\b")
}
flag_detected <- function(txt) {
  t <- str_to_lower(txt %||% "")
  str_detect(t, "\\bdetect(ed|ion)\\b|\\bfound\\b|\\breported\\b|\\bobserved\\b")
}
flag_mortality <- function(txt) {
  t <- str_to_lower(txt %||% "")
  str_detect(t, "\\bmortality\\b|\\bdead\\b|\\bdieback\\b|\\bkilled\\b|\\bdecline\\b")
}
flag_not_detected <- function(txt) {
  t <- stringr::str_to_lower(txt %||% "")
  stringr::str_detect(
    t,
    "\\bnot detected\\b|\\bno detections?\\b|\\bno detection\\b|\\bnot found\\b|\\bno (?:sign|signs) of\\b|\\babsent\\b"
  )
}

# splitting up remarks by state (use caution with these)
# state token regex for remark segmentation
state_abbr <- names(abbr_to_name)
state_names <- unique(unname(abbr_to_name))
abbr_re <- paste0("\\b(", paste(state_abbr, collapse="|"), ")\\b")
name_re <- paste0("\\b(", paste(sort(state_names, decreasing = TRUE), collapse="|"), ")\\b")

token_to_state <- function(tok) {
  tok <- str_squish(tok)
  up <- str_to_upper(tok)
  if (up %in% names(abbr_to_name)) return(unname(abbr_to_name[up]))
  str_to_title(str_to_lower(tok))
}

split_remarks_by_state <- function(txt) {
  txt <- replace_na(txt, "")
  txt <- str_replace_all(txt, "\r", "\n")
  txt <- str_squish(txt)
  if (txt == "") return(tibble(state_token=character(), segment=character()))

  marker_re <- paste0("(?i)(?<=\\s|^)", "(", abbr_re, "|", name_re, ")", "\\s*[:\\-—]")
  marked <- str_replace_all(txt, marker_re, "||\\1:")

  parts <- str_split(marked, "\\|\\|", simplify = FALSE)[[1]]
  parts <- parts[parts != ""]

  if (length(parts) == 1 && !str_detect(parts[1], paste0("(?i)^(", abbr_re, "|", name_re, "):"))) {
    return(tibble(state_token=character(), segment=character()))
  }

  tibble(raw = parts) %>%
    mutate(
      state_token = str_match(raw, paste0("(?i)^(", abbr_re, "|", name_re, "):"))[,2],
      segment = str_trim(str_replace(raw, paste0("(?i)^(", abbr_re, "|", name_re, "):"), ""))
    ) %>%
    filter(!is.na(state_token), segment != "")
}


split_remarks_by_state_markers <- function(txt) {
  txt <- replace_na(txt, "")
  txt <- str_replace_all(txt, "\r", "\n")
  txt <- str_squish(txt)
  if (txt == "") return(tibble(State=character(), segment=character(), segment_source=character()))

  marker_re <- paste0("(?i)(?<=\\s|^)", "(", state_abbr_re, "|", state_name_re, ")", "\\s*[:\\-—]")
  marked <- str_replace_all(txt, marker_re, "||\\1:")

  parts <- str_split(marked, "\\|\\|", simplify = FALSE)[[1]]
  parts <- parts[parts != ""]

  if (length(parts) == 1 && !str_detect(parts[1], paste0("(?i)^(", state_abbr_re, "|", state_name_re, "):"))) {
    return(tibble(State=character(), segment=character(), segment_source=character()))
  }

  tibble(raw = parts) %>%
    mutate(
      state_token = str_match(raw, paste0("(?i)^(", state_abbr_re, "|", state_name_re, "):"))[,2],
      segment = str_trim(str_replace(raw, paste0("(?i)^(", state_abbr_re, "|", state_name_re, "):"), "")),
      State = map_chr(state_token, token_to_state),
      segment_source = "marker"
    ) %>%
    select(State, segment, segment_source) %>%
    filter(segment != "")
}

split_remarks_by_state_inline <- function(txt) {
  txt <- replace_na(txt, "")
  txt <- str_squish(txt)
  if (txt == "") return(tibble(State=character(), segment=character(), segment_source=character()))

  # Split into clauses/sentences 
  clauses <- unlist(str_split(txt, "(?<=[\\.;!?])\\s+|\\s*;\\s*"))
  clauses <- str_squish(clauses)
  clauses <- clauses[clauses != ""]

  out <- map_dfr(clauses, function(cl) {
    states1 <- str_extract_all(cl, regex(state_name_re, ignore_case = TRUE))[[1]]
    states2 <- str_extract_all(cl, regex(state_abbr_re, ignore_case = FALSE))[[1]]
    states  <- unique(c(states1, states2))
    if (length(states) == 0) return(NULL)

    tibble(
      State = map_chr(states, token_to_state),
      segment = cl,
      segment_source = "inline"
    )
  })

  out
}

get_state_segments <- function(txt) {
  seg_marker <- split_remarks_by_state_markers(txt)
  if (nrow(seg_marker) > 0) return(seg_marker)

  # fallback to inline mentions
 split_remarks_by_state_inline(txt)
}


# -------------------------------------------------------
# 5) read in the raw, corrected .csv files  

# read in yearly csv files:
files <- list.files("data/FHM_cor_yearfiles", pattern = "\\.csv$", full.names = TRUE)

FHM.df <- map_dfr(files, function(f){
  yr <- str_extract(basename(f), pattern = "\\d+") %>% as.integer()

  read_csv(f, show_col_types = FALSE) %>%
    select(Year, Agent_Type, Agent_Name, Host, Location, Remarks) %>% 
    mutate(year.filename = yr) %>% # if for some reason year is not populated, get it from filename
    mutate(Year = ifelse(is.na(Year), yr, Year))
})
FHM.df %>% filter(is.na(Year))
unique(FHM.df$Year)

# save this file to box:
output.dir <- "C:/Users/KellyHeilman/Box/01. kelly.heilman Workspace/mortality/Eastern-Mortality/mortality_models/"
yr1 <- min(unique(FHM.df$Year))
yr2 <- max(unique(FHM.df$Year))
write.csv(FHM.df, paste0(output.dir, "data/Forest_condition_digitized_", yr1,"_", yr2, ".csv"))

# --------------------------------------------------------------------
# reformat FHM.df from raw, corrected transcription to a long format:
# expand based on hosts and states



# 1) expand each host to multiple rows---
FHM_hosts_long <- FHM.df %>%
  mutate(
    Host_raw = Host,
    Host_std = map(Host_raw, standardize_host_row)
  ) %>%
  unnest(Host_std) %>%
  #filter(!is.na(Host_std), Host_std != "") %>%
  mutate(Host = Host_std) %>%
  select(Year, Agent_Type, Agent_Name, Host, Host_raw, Location, Remarks)

# check the hosts
audit.hosts <- FHM_hosts_long %>% select(Host, Host_raw) %>% distinct()
audit.hosts %>% View() # there may still be some issues, but this is better
unique(FHM_hosts_long$Host_raw)
unique(FHM_hosts_long$Host)
# one issue that we should correct in corrected csvs:
FHM_hosts_long %>% filter(Host_raw %in% "Conifers and hardwood. Particularly black, white, red, pin, and chestnut oak") 

# 2) expand locations to mulitiple State rows
FHM_host_state_long <- FHM_hosts_long %>%
  filter(!is.na(Host)) %>%
  mutate(State = map(Location, expand_location)) %>%
  unnest(State) %>%
  select(Year, Agent_Type, Agent_Name, Host, State, Remarks)

unique(FHM_host_state_long$State)


# 3) build STATE-SEGMENT metrics table (rename metrics with _seg to avoid collisions)
remarks_state_metrics <- FHM_host_state_long %>%
  distinct(Year, Agent_Type, Agent_Name, Remarks) %>%
  mutate(seg = map(Remarks, get_state_segments)) %>%
  unnest(seg) %>%
  
  # combine state segments together, if they exist:
  group_by(Year, Agent_Type, Agent_Name, State, Remarks) %>%
  summarise(segment_collapsed = paste(segment, collapse = "; "), 
            segment_source_collapsed = paste(unique(segment_source), collapse = "; ")) %>%
  
  transmute(
    Year, Agent_Type, Agent_Name,  State,
    segment = segment_collapsed,
    segment_source = segment_source_collapsed,

    increasing_seg = flag_increasing(segment_collapsed),
    decreasing_seg = flag_decreasing(segment_collapsed),
    stable_seg     = flag_stable(segment_collapsed),
    detected_seg   = flag_detected(segment_collapsed),
    mortality_seg  = flag_mortality(segment_collapsed),

    acres_seg            = extract_acres_any(segment_collapsed),
    acres_defoliated_seg = extract_acres_defoliated(segment_collapsed),
    percent_seg          = extract_percent(segment_collapsed),
    dead_trees_seg       = extract_dead_trees(segment_collapsed),
    trap_count_seg       = extract_trap_count(segment_collapsed),
    egg_masses_seg       = extract_egg_masses(segment_collapsed)
  )

# state segments are looking betting
remarks_state_metrics %>% View()

# condense to 1 row
remarks_state_metrics_1row <- remarks_state_metrics %>%
  group_by(Year, Agent_Type, Agent_Name, State) %>%
  summarise(
    # keep text for QA (optional; can be long)
    segment = paste(unique(na.omit(segment)), collapse = " | "),
    segment_source = paste(unique(na.omit(segment_source)), collapse = ","),

    # logical: TRUE if any segment indicates it
    increasing_seg = any(increasing_seg, na.rm = TRUE),
    decreasing_seg = any(decreasing_seg, na.rm = TRUE),
    stable_seg     = any(stable_seg,     na.rm = TRUE),
    detected_seg   = any(detected_seg,   na.rm = TRUE),
    mortality_seg  = any(mortality_seg,  na.rm = TRUE),

    # numeric: take max across segments (safer than sum unless you *know* segments are disjoint)
    acres_seg            = suppressWarnings(max(acres_seg, na.rm = TRUE)),
    acres_defoliated_seg = suppressWarnings(max(acres_defoliated_seg, na.rm = TRUE)),
    percent_seg          = suppressWarnings(max(percent_seg, na.rm = TRUE)),
    dead_trees_seg       = suppressWarnings(max(dead_trees_seg, na.rm = TRUE)),
    trap_count_seg       = suppressWarnings(max(trap_count_seg, na.rm = TRUE)),
    egg_masses_seg       = suppressWarnings(max(egg_masses_seg, na.rm = TRUE)),

    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))

# 4) build table with row flags for the remarks:
row_level_metrics <- FHM_host_state_long %>%
  distinct(Year, Agent_Type, Agent_Name, Host, Remarks) %>%
  transmute(
    Year, Agent_Type, Agent_Name, Host, Remarks,
    increasing_row = flag_increasing(Remarks),
    decreasing_row = flag_decreasing(Remarks),
    stable_row     = flag_stable(Remarks),
    detected_row   = flag_detected(Remarks),
    mortality_row  = flag_mortality(Remarks)#,
   # acres_row      = extract_acres_any(Remarks),
  #  dead_trees_row = extract_dead_trees(Remarks)
  )



# 5) join the FHM_host_state_long file to the state remarks
FHM_host_state_enriched <- FHM_host_state_long %>%
  left_join(remarks_state_metrics,
            by = c("Year","Agent_Type","Agent_Name","State")) %>%
  left_join(row_level_metrics,
            by = c("Year","Agent_Type","Agent_Name","Host","Remarks")) %>%
  
  # use coalesce to use the state-level segments if not NA, but use row if NA
  mutate(
    increasing = coalesce(increasing_seg, increasing_row),
    decreasing = coalesce(decreasing_seg, decreasing_row),
    stable     = coalesce(stable_seg,     stable_row),
    detected   = coalesce(detected_seg,   detected_row),
    mortality  = coalesce(mortality_seg,  mortality_row),
# just use the state level segment values for the numerics
    acres            = acres_seg, 
    acres_defoliated = acres_defoliated_seg,   # row-level version optional
    percent          = percent_seg,
    dead_trees       = dead_trees_seg,
    trap_count       = trap_count_seg,
    egg_masses       = egg_masses_seg
  )


FHM_host_state_long

# cleaning up the agent names and state names:
FHM_host_state_enriched <- FHM_host_state_enriched %>%
  mutate(
    Agent_raw = Agent_Name,
    Agent_std = map_chr(Agent_raw, standardize_agent)
  )
FHM_host_state_enriched <- FHM_host_state_enriched %>%
  mutate(
    State_clean = map(State, normalize_state_field)
  ) %>%
  unnest(State_clean) %>%
  select(-State) %>%
  rename(State = State_clean)


# some additional, post-hoc flags:
FHM_host_state_enriched <- FHM_host_state_enriched %>%
  mutate(
    not_detected = flag_not_detected(coalesce(segment, Remarks))
  )
FHM_host_state_enriched <- FHM_host_state_enriched %>%
  mutate(
    active_any = !(decreasing %in% TRUE) & !(not_detected %in% TRUE)
  )

# check the standardized and raw agents:
audit <- FHM_host_state_enriched %>%
  distinct(Agent_raw, Agent_std) %>%
  arrange(Agent_std, Agent_raw)
#audit %>% View()

# check to make sure there are no duplicated types of an agent
check.agent.type <- FHM_host_state_enriched %>%
  distinct(Agent_Type, Agent_std) %>%
  group_by(Agent_std) %>%
  summarise(should.b.1 = n()) %>%
  filter(should.b.1 > 1)

nrow(check.agent.type)


unique(FHM_host_state_enriched$Agent_std)

# show a few inline-segment examples
FHM_host_state_enriched %>%
  filter(segment_source == "inline") %>%
  select(Year, Agent_Name, Host, State, segment) %>%
  slice_head(n = 10)

# how many rows have state-specific segments?
FHM_host_state_enriched %>%
  summarise(pct_with_segment = mean(!is.na(segment)) * 100)

# look at cases where Location expanded to a state but we had no segment
FHM_host_state_enriched %>%
  filter(is.na(segment)) %>%
  select(Year, Agent_Name, Host, State, Remarks) %>%
  slice_head(n = 10)

# save 
saveRDS(FHM_host_state_enriched, paste0(output.dir, "data/FHM_state_host_long.rds"))

########################################################
# summarize trends in NE:-----

trend_summary_by_host_pest_year <- FHM_host_state_enriched %>%
  filter(State %in% c(northeast_states, "West Virginia")) %>%
  group_by(Year, Host, Agent_std, Agent_Type) %>%
  summarise(
    n_states = n_distinct(State),

    any_increasing = any(increasing %in% TRUE, na.rm = TRUE),
    any_mortality  = any(mortality  %in% TRUE, na.rm = TRUE),
    any_decreasing = any(decreasing %in% TRUE, na.rm = TRUE),
    any_detected   = any(detected   %in% TRUE, na.rm = TRUE),

    
    increasing_sum      = sum(increasing %in% TRUE, na.rm = TRUE),
    mortality_sum      = sum(mortality %in% TRUE, na.rm = TRUE),
    detected_sum      = sum(detected %in% TRUE, na.rm = TRUE),
    
    # optional magnitude context
    acres_sum      = sum(acres, na.rm = TRUE),
    defol_acres_sum = sum(acres_defoliated, na.rm = TRUE),
    dead_trees_sum  = sum(dead_trees, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    across(c(acres_sum, defol_acres_sum, dead_trees_sum),
           ~ ifelse(is.infinite(.x), NA_real_, .x))
  ) %>%
  arrange(Year, Host, desc(any_mortality), desc(any_increasing))


# categorical variables
trend_summary_by_host_pest_year <- trend_summary_by_host_pest_year %>%
  mutate(
    trend_category = case_when(
      any_mortality & any_increasing ~ "Increasing + Mortality",
      any_mortality                  ~ "Mortality",
      any_increasing                 ~ "Increasing",
      any_decreasing                 ~ "Decreasing",
      any_detected                   ~ "Detected",
      TRUE                           ~ "Mentioned"
    )
  )



plot_FHM <- trend_summary_by_host_pest_year %>%
  group_by(Year, Agent_std) %>%
  summarise(
    any_increasing = any(any_increasing, na.rm = TRUE),
    any_mortality  = any(any_mortality,  na.rm = TRUE),
    n_states = max(n_states, na.rm = TRUE),   # optional context if present
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      any_mortality & any_increasing ~ "Increasing + Mortality",
      any_mortality                  ~ "Mortality",
      any_increasing                 ~ "Increasing",
      TRUE                           ~ "Mentioned/Other"
    )
  )

# Order pests by how often they show up as increasing/mortality (so the heatmap is readable)

# alternatively, order by pest types:
pest.type.df <- read.csv("data/northeast_forest_agents_types.csv")

pest_order <- plot_FHM %>%
  group_by(Agent_std) %>%
  summarise(
    years_increasing = sum(any_increasing, na.rm = TRUE),
    years_mortality  = sum(any_mortality,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(years_mortality), desc(years_increasing)) %>%
  pull(Agent_std)


plot_FHM <- plot_FHM %>%
  mutate(Agent_std = factor(Agent_std, levels = pest_order)) %>% 
  left_join(., pest.type.df)


ggplot(plot_FHM, aes(x = Year, y = Agent_std, fill = status)) +
  geom_tile(color = NA) +
  labs(
    title = "Temporal Heatmap of Forest Pests (All Hosts Combined)",
    x = "Year",
    y = "Forest pest (standardized)",
    fill = "Signal"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )






pest_order_class <- plot_FHM %>%
  distinct(Agent_std, Class) %>%
  group_by(Agent_std) %>%
  arrange(desc(Class)) %>%
  pull(Agent_std)

plot_FHM <- plot_FHM %>%
  mutate(Agent_std = factor(Agent_std, levels = pest_order_class)) 

ggplot(plot_FHM, aes(x = Year, y = Agent_std, fill = status)) +
  geom_tile(color = NA) +
  labs(
    title = "Temporal Heatmap of Forest Pests (All Hosts Combined)",
    x = "Year",
    y = "Forest pest (standardized)",
    fill = "Signal"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )









pest_year_increasing <- trend_summary_by_host_pest_year %>%
  group_by(Year, Agent_std) %>%
  summarise(
    any_increasing = any(any_increasing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_increasing) %>%
  mutate(value = 1)

ggplot(pest_year_increasing, aes(x = factor(Year), y = value, fill = Agent_std)) +
  geom_col() +
  labs(
    title = "Forest pests mentioned as increasing by year (presence)",
    x = "Year",
    y = "Number of pests flagged increasing",
    fill = "Pest"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pest_year_increasing_states <- trend_summary_by_host_pest_year %>%
  filter(any_increasing) %>%
  group_by(Year, Agent_std) %>%
  summarise(
    states_sum = sum(n_states, na.rm = TRUE),   # sums across hosts
    .groups = "drop"
  )

ggplot(pest_year_increasing_states, aes(x = factor(Year), y = states_sum, fill = Agent_std)) +
  geom_col() +
  labs(
    title = "Forest pests mentioned as increasing by year (weighted by # states)",
    x = "Year",
    y = "Sum of state-mentions (across hosts)",
    fill = "Pest"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


increasing_counts <- trend_summary_by_host_pest_year %>%
  group_by(Year, Host) %>%
  summarise(
    n_increasing = sum(any_increasing, na.rm = TRUE),
    n_mortality  = sum(any_mortality,  na.rm = TRUE),
    .groups = "drop"
  )

long_term_trends <- trend_summary_by_host_pest_year %>%
  group_by(Host, Agent_std) %>%
  summarise(
    years_increasing = sum(any_increasing, na.rm = TRUE),
    years_mortality  = sum(any_mortality,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Host, desc(years_mortality), desc(years_increasing))


# attempts to quantify intensity--many states-species are missing values so not really infromative
intensity_by_host_pest_year <- FHM_host_state_enriched %>%
  mutate(
    # numeric components (log-scaled so big numbers don't dominate too hard)
    acres_term      = ifelse(is.na(acres), 0, log10(1 + acres)),
    defol_term      = ifelse(is.na(acres_defoliated), 0, log10(1 + acres_defoliated)),
    dead_term       = ifelse(is.na(dead_trees), 0, log10(1 + dead_trees)),
    trap_term       = ifelse(is.na(trap_count), 0, log10(1 + trap_count)),
    egg_term        = ifelse(is.na(egg_masses), 0, log10(1 + egg_masses)),
    percent_term    = ifelse(is.na(percent), 0, percent / 100),

    # boolean components (weights you can tune)
    trend_term =
      2 * as.integer(increasing %in% TRUE) +
      2 * as.integer(mortality  %in% TRUE) +
      1 * as.integer(detected   %in% TRUE) -
      1 * as.integer(decreasing %in% TRUE),

    # per-state intensity
    intensity_state =
      trend_term +
      2.0 * defol_term +
      1.5 * acres_term +
      1.5 * dead_term +
      1.0 * percent_term +
      0.8 * trap_term +
      0.8 * egg_term
  ) %>%
  group_by(Year, Host, Agent_std, Agent_Type) %>%
  summarise(
    n_states = n_distinct(State),

    # aggregate intensity across states (mean = "typical"; sum = "total burden")
    intensity_mean = mean(intensity_state, na.rm = TRUE),
    intensity_sum  = sum(intensity_state,  na.rm = TRUE),

    # keep some interpretable aggregates too
    any_increasing = any(increasing %in% TRUE, na.rm = TRUE),
    any_mortality  = any(mortality  %in% TRUE, na.rm = TRUE),

    acres_sum = sum(acres, na.rm = TRUE),
    defol_acres_sum = sum(acres_defoliated, na.rm = TRUE),
    dead_trees_sum  = sum(dead_trees, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    across(c(acres_sum, defol_acres_sum, dead_trees_sum), ~ ifelse(is.infinite(.x), NA_real_, .x))
  ) %>%
  arrange(Year, Host, desc(intensity_sum))

top10_by_host_year <- intensity_by_host_pest_year %>%
  group_by(Year, Host) %>%
  slice_max(intensity_sum, n = 10, with_ties = FALSE) %>%
  ungroup()

top_overall_by_host <- intensity_by_host_pest_year %>%
  group_by(Host, Agent_std, Agent_Type) %>%
  summarise(
    years_present = n_distinct(Year),
    intensity_total = sum(intensity_sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Host, desc(intensity_total))

top_overall_by_host %>% filter(Host %in% "american beech")
top_overall_by_host %>% filter(Host %in% "black cherry")
top_overall_by_host %>% filter(Host %in% "balsam fir")


library(sf)
library(dplyr)
library(ggplot2)
library(tigris)
library(stringr)


options(tigris_use_cache = TRUE)

states_sf <- tigris::states(cb = TRUE, year = 2022) %>%
  st_transform(5070) %>%                 # US Albers equal-area
  mutate(State = NAME) %>%
  select(State, geometry)

library(dplyr)
library(stringr)
library(tibble)











make_host_meta <- function(hosts) {
  tibble(Host = sort(unique(hosts))) %>%
    mutate(Host = str_to_lower(str_squish(Host)),

      # host "type"
      type_group = case_when(
        Host %in% c("conifers") ~ "conifer",
        Host %in% c("hardwoods") ~ "hardwood",
        Host %in% c("all trees", "various species") ~ "all",

        Host %in% c(
          "fir","balsam fir",
          "spruce","black spruce","red spruce","white spruce","black spruce", "blue spruce", "norway spruce",
          "pine","eastern white pine","jack pine","red pine","scotch pine", "loblolly pine",
          "austrian pine","shortleaf pine", "virginia pine",
          "eastern larch","european larch","tamarack",
          "eastern hemlock", "red cedar", "white cedar"
        ) ~ "conifer",
        Host %in% c(
          "maple","sugar maple","silver maple", "red maple",
          "oak", "pin oak", "white oak", "red oak", "black oak",
          "elm",
          "ash","white ash", "black ash", "green ash",
          "aspen","poplar",
          "birch","yellow birch","paper birch","grey birch",
          "cherry","black cherry",
          "black walnut",
          "american beech",
          "basswood","hickory","sycamore","dogwood", "hardwoods", "butternut"
        ) ~ "hardwood",
        TRUE ~ "unknown"
      ),

      # your "genus" grouping (really: genus/common-group label)
      genus_group = case_when(
        str_detect(Host, "pine") | Host == "pine" ~ "pine",
        str_detect(Host, "spruce") | Host == "spruce" ~ "spruce",
        str_detect(Host, "fir") | Host == "fir" ~ "fir",
        Host %in% c("tamarack","eastern larch","european larch") ~ "larch",
        Host == "birch" | str_detect(Host, "birch") ~ "birch",
        Host == "maple"| Host == "boxelder" | str_detect(Host, "maple") ~ "maple",
        Host == "cherry" | str_detect(Host, "cherry") ~ "cherry",
        Host == "oak" | str_detect(Host, "oak") ~ "oak",
        Host == "butternut" | str_detect(Host, "walnut") ~ "walnut",
        Host == "cedar" | str_detect(Host, "cedar") ~ "cedar",
        Host == "elm" ~ "elm",
        Host == "ash" ~ "ash",
        Host == "aspen" ~ "aspen",
        Host == "american beech" ~ "beech",
        Host == "basswood" ~ "basswood",
        Host == "hickory" ~ "hickory",
        Host == "sycamore" ~ "sycamore",
        Host == "flowering dogwood" ~ "dogwood",
        Host == "eastern hemlock" ~ "hemlock",
        TRUE ~ NA_character_
      ),

      host_level = case_when(
        Host %in% c("conifers","hardwoods","all trees","various species") ~ "type_group",
        # inside make_host_meta() case_when for type_group
        Host %in% c("all trees","various species") ~ "all",
        Host %in% c("pine","spruce","fir","maple","birch","cherry","oak","elm","ash", "cedar", "walnut") ~ "genus_group",
        TRUE ~ "species"
      )
    )
}

host_meta <- make_host_meta(FHM_host_state_enriched$Host)


rollup_hosts <- function(FHM, host_meta,
                         rollup_genus = TRUE,
                         rollup_type  = FALSE,
                         rollup_all   = TRUE) {

  FHM2 <- FHM %>%
    mutate(Host = str_to_lower(str_squish(Host))) %>%
    left_join(host_meta, by = "Host")

  base <- FHM2 %>%
    transmute(Year, State, Host_orig = Host, Host_plot = Host,
              Agent_std, Agent_Type,
              increasing, mortality, detected,
              acres_defoliated, acres, dead_trees)


  out <- base

  # genus-group rows -> all species in that genus
  if (rollup_genus) {
    genus_species_map <- host_meta %>%
      filter(!is.na(genus_group), host_level == "species") %>%
      select(genus_group, Host_plot = Host)

    genus_rows <- FHM2 %>%
      filter(host_level == "genus_group", !is.na(genus_group)) %>%
      select(Year, State, Host_orig = Host, genus_group,
             Agent_std, Agent_Type,
             increasing, mortality, detected,
             acres_defoliated, acres, dead_trees) %>%
      left_join(genus_species_map, by = "genus_group") %>%
      select(-genus_group)

    out <- bind_rows(out, genus_rows)
  }

  # type-group rows (hardwoods/conifers) -> all species in that type
  if (rollup_type) {
    type_species_map <- host_meta %>%
      filter(type_group %in% c("hardwood","conifer"), host_level == "species") %>%
      select(type_group, Host_plot = Host)

    type_rows <- FHM2 %>%
      filter(host_level == "type_group", type_group %in% c("hardwood","conifer")) %>%
      select(Year, State, Host_orig = Host, type_group,
             Agent_std, Agent_Type,
             increasing, mortality, detected,
             acres_defoliated, acres, dead_trees) %>%
      left_join(type_species_map, by = "type_group") %>%
      select(-type_group)

    out <- bind_rows(out, type_rows)
  }

  # areawide rows (all trees/various species) -> ALL species in dataset
  if (rollup_all) {
    all_species <- host_meta %>%
      filter(host_level == "species") %>%
      pull(Host) %>%
      unique()

    all_rows <- FHM2 %>%
      filter(type_group == "all") %>%
      select(Year, State, Host_orig = Host,
             Agent_std, Agent_Type,
             increasing, mortality, detected,
             acres_defoliated, acres, dead_trees) %>%
      tidyr::crossing(Host_plot = all_species)

    out <- bind_rows(out, all_rows)
  }

  out %>% distinct()
}

FHM_roll <- rollup_hosts(FHM_host_state_enriched, host_meta,
                        rollup_genus = TRUE,
                        rollup_type  = FALSE)  # set TRUE if you want "hardwoods" -> all hardwood species

get_host_set <- function(target_host,
                         host_meta,
                         include_species = TRUE,
                         include_genus   = TRUE,
                         include_type    = FALSE,
                         include_alltrees = FALSE) {

  target_host <- str_to_lower(str_squish(target_host))
  if (!target_host %in% host_meta$Host) stop("target_host not found in host_meta$Host: ", target_host)

  row <- host_meta %>% filter(Host == target_host) %>% slice(1)

  out <- character(0)

  if (include_species) {
    out <- c(out, row$Host)
  }

  if (include_genus && !is.na(row$genus_group)) {
    out <- c(out,
             row$genus_group,  # include the group label itself, e.g. "fir"
             host_meta %>% filter(genus_group == row$genus_group) %>% pull(Host))
  }

  if (include_type && !is.na(row$type_group) && row$type_group %in% c("conifer","hardwood")) {
    out <- c(out,
             if (row$type_group == "conifer") "conifers" else "hardwoods",
             host_meta %>% filter(type_group == row$type_group) %>% pull(Host))
  }

  if (include_alltrees) {
    out <- c(out, "all trees", "various species")
  }

  unique(out)
}


get_host_set <- function(target_host, host_meta,
                         include_species = TRUE,
                         include_genus   = TRUE,
                         include_type    = FALSE,
                         include_alltrees = FALSE) {

  target_host <- str_to_lower(str_squish(target_host))
  if (!target_host %in% host_meta$Host) stop("target_host not found in host_meta: ", target_host)

  row <- host_meta %>% filter(Host == target_host) %>% slice(1)
  out <- character(0)

  if (include_species) out <- c(out, row$Host)

  if (include_genus && !is.na(row$genus_group)) {
    out <- c(out, row$genus_group, host_meta %>% filter(genus_group == row$genus_group) %>% pull(Host))
  }

  if (include_type && row$type_group %in% c("conifer","hardwood")) {
    out <- c(out,
             if (row$type_group == "conifer") "conifers" else "hardwoods",
             host_meta %>% filter(type_group == row$type_group) %>% pull(Host))
  }

  if (include_alltrees) out <- c(out, "all trees","various species")

  unique(out)
}

plot_pest_stacked_by_host <- function(FHM_roll, target_host, host_meta,
                                      include_species = TRUE,
                                      include_genus   = TRUE,
                                      include_type    = FALSE,
                                      include_alltrees = FALSE,
                                      signal = c("increasing","mortality","defoliation_high","decline_complex","detected"),
                                      defol_thresh = 1000,
                                      weight = c("count","states"),
                                      top_n_pests = 15, 
                                      fill = c("Agent_std", "Agent_type")) {

  signal <- match.arg(signal)
  weight <- match.arg(weight)

  host_set <- get_host_set(target_host, host_meta,
                           include_species, include_genus, include_type, include_alltrees)

  d <- FHM_roll %>%
    filter(Host_plot %in% host_set)

  # define the event flag we’re plotting
  d <- d %>%
    mutate(
      event = case_when(
        signal == "increasing" ~ (increasing %in% TRUE),
        signal == "mortality"  ~ (mortality  %in% TRUE),
        signal == "detected"   ~ (detected   %in% TRUE),
        signal == "defoliation_high" ~ (!is.na(acres_defoliated) & acres_defoliated >= defol_thresh),
        signal == "decline_complex" ~ str_detect(str_to_lower(Agent_std), "decline complex"),
        TRUE ~ FALSE
      )
    ) %>%
    filter(event)

  # aggregate to Year x Host_plot x Pest
  d_sum <- d %>%
    group_by(Year, Host_plot, Agent_std, Agent_Type) %>%
    summarise(
      n_states = n_distinct(State),
      value = if (weight == "states") n_states else 1,
      .groups = "drop"
    )

  if(fill %in% 'Agent_std'){
  # limit pests for readability
  top_pests <- d_sum %>%
    group_by(Agent_std, Agent_Type) %>%
    summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
    slice_max(total, n = top_n_pests) %>%
    pull(Agent_std)

  d_plot <- d_sum %>%
    mutate(Pest = ifelse(Agent_std %in% top_pests, Agent_std, "Other")) %>%
    group_by(Year, Host_plot, Pest) %>%
    summarise(value = sum(value), .groups = "drop")

  ggplot(d_plot, aes(x = Year, y = value, fill = Pest)) +
    geom_col() +
    facet_wrap(~ Host_plot, scales = "free_y") +
    labs(
      title = paste0("Pest signals over time — host selection: ", target_host),
      subtitle = paste0(
        "Signal = ", signal,
        if (signal == "defoliation_high") paste0(" (≥ ", defol_thresh, " defoliated acres)") else "",
        " | y = ", ifelse(weight == "states", "# states", "count/presence"),
        " | genus roll-up applied"
      ),
      x = "Year",
      y = ifelse(weight == "states", "Number of states with signal", "Count of pests with signal"),
      fill = "Pest"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
    xlim(min(FHM_roll$Year), max(FHM_roll$Year))
  
  }else{
  # otherwise use the agent type to plot
    # limit pests for readability
    top_pests <- d_sum %>%
      group_by(Agent_std, Agent_Type) %>%
      summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
      slice_max(total, n = top_n_pests) %>%
      pull(Agent_std)
    
 
    
    ggplot(d_sum, aes(x = Year, y = value, fill = Agent_Type)) +
      geom_col() +
      facet_wrap(~ Host_plot, scales = "free_y") +
      labs(
        title = paste0("Pest signals over time — host selection: ", target_host),
        subtitle = paste0(
          "Signal = ", signal,
          if (signal == "defoliation_high") paste0(" (≥ ", defol_thresh, " defoliated acres)") else "",
          " | y = ", ifelse(weight == "states", "# states", "count/presence"),
          " | genus roll-up applied"
        ),
        x = "Year",
        y = ifelse(weight == "states", "Number of states with signal", "Count of pests with signal"),
        fill = "Agent type"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
      xlim(min(FHM_roll$Year), max(FHM_roll$Year))
      
  }
}

# try to do this for all of our common names:
nspp <- data.frame(SPCD = c(316, 318, 833, 832, 261, 531, 802, 129, 762,  12, 541,  97, 621, 400, 371, 241, 375))
nspp$Species <- paste(FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$GENUS, FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$SPECIES)
# link up to the species table:
nspp$COMMON <- FIESTA::ref_species[match(nspp$SPCD, FIESTA::ref_species$SPCD),]$COMMON

species.of.interest <- c("red maple", "sugar maple","red oak", "white oak","chestnut oak", "eastern hemlock", "american beech", "eastern white pine", 
  "black cherry", "balsam fir", "white ash", "red spruce", "yellow poplar", "hickory", "yellow birch", "white cedar", "paper birch")


# first lets filter the trend summary:
FHM_interest <- trend_summary_by_host_pest_year %>% 
  filter(Host %in% species.of.interest) %>%
  group_by(Year, Agent_std, Agent_Type) %>%
  summarise(
    any_increasing = any(any_increasing, na.rm = TRUE),
    any_mortality  = any(any_mortality,  na.rm = TRUE),
    n_states = max(n_states, na.rm = TRUE),   # optional context if present
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      any_mortality & any_increasing ~ "Increasing + Mortality",
      any_mortality                  ~ "Mortality",
      any_increasing                 ~ "Increasing",
      TRUE                           ~ "Mentioned/Other"
    )
  )

# Order pests by how often they show up as increasing/mortality (so the heatmap is readable)
FHM_interest

# alternatively, order by pest types:
pest.type.df <- read.csv("data/northeast_forest_agents_types.csv")

pest_order <- FHM_interest %>%
  group_by(Agent_std) %>%
  summarise(
    years_increasing = sum(any_increasing, na.rm = TRUE),
    years_mortality  = sum(any_mortality,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(years_mortality), desc(years_increasing)) %>%
  pull(Agent_std)

ggplot(FHM_interest, aes(x = Year, y = Agent_std, fill = status)) +
  geom_tile(color = NA) +
  labs(
    title = "Temporal Heatmap of Forest Pests (All Hosts  of interest)",
    x = "Year",
    y = "Forest pest (standardized)",
    fill = "Signal"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )+facet_wrap(~Agent_Type)


allbarplot.spp <- lapply(species.of.interest, function(spp){
  plot_pest_stacked_by_host(
    FHM_roll, spp, host_meta,
    include_species = TRUE,
    include_genus = FALSE,
    include_type = FALSE,
    signal = "increasing",
    weight = "states",
    top_n_pests = 15, 
    fill = c("Agent_std")
  )
  ggsave(paste0(output.dir, "images/top_pests_barplot_", spp, ".png"))
})

allbarplot.spp <- lapply(species.of.interest, function(spp){
  plot_pest_stacked_by_host(
    FHM_roll, spp, host_meta,
    include_species = TRUE,
    include_genus = FALSE,
    include_type = FALSE,
    signal = "increasing",
    weight = "count",
    top_n_pests = 15, 
    fill = c("Agent_type")
  )
  ggsave(paste0(output.dir, "images/top_agent_types_barplot_", spp, ".png"))
})


# next:
# plot up timeseries by Agent_type using the top most common pests/pathogens for each species





# plot mentions over time in NE

# Define your Northeastern set (edit if needed)
northeast_states <- c(
  "Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut",
  "New York","New Jersey","Pennsylvania","Delaware","Maryland","Ohio","West Virginia"
)

heat_df <- FHM_host_state_enriched %>%
  # keep only records with location data
  filter(!is.na(State), State != "") %>%
  # NE only
  filter(State %in% northeast_states) %>%
  # collapse to Year x Agent presence (any mention anywhere in NE that year)
  group_by(Year, Agent_std) %>%
  summarise(present = 1L, .groups = "drop")

# Order agents for readability (most frequently present at top)
agent_order <- heat_df %>%
  count(Agent_std, name = "n_years") %>%
  arrange(desc(n_years)) %>%
  pull(Agent_std)

heat_df <- heat_df %>%
  mutate(Agent_std = factor(Agent_std, levels = agent_order))



ne_host_state_enriched <- FHM_host_state_enriched %>%
  filter(!is.na(State), State != "", State %in% northeast_states)
ne_host_state_enriched %>% filter(Agent_std %in% "Spruce Budworm")%>% View()

heat_df_states <- FHM_host_state_enriched %>%
  filter(!is.na(State), State != "", State %in% northeast_states) %>%
  group_by(Year, Agent_std) %>%
  summarise(n_states = n_distinct(State), .groups = "drop")

agent_order <- heat_df_states %>%
  summarise(n_years = n(), .by = Agent_std) %>%
  arrange(desc(n_years)) %>%
  pull(Agent_std)

heat_df_states <- heat_df_states %>%
  mutate(Agent_std = factor(Agent_std, levels = agent_order))

state.heatmap.all <- ggplot(heat_df_states %>% filter(!is.na(Year)), aes(x = as.character(Year), y = Agent_std, fill = n_states)) +
  geom_tile(color = "lightgrey") +
  labs(
    title = "Forest pests and conditions over time",
    subtitle = "number of NE states with a mention that year",
    x = "Year",
    y = "Agent",
    fill = "# states"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())+
  scale_fill_viridis_c()
state.heatmap.all
ggsave(paste0(output.dir, "images/heat_map_NE_pests_diseases_time_state.png"), plot = state.heatmap.all,
       height = 8, width = 7)

## do this but group by hosts:

plot_agent_heatmap_by_host <- function(df,
                                       host,
                                       host_meta,
                                       fill = c("presence", "states"),
                                       include_species = TRUE,
                                       include_genus   = TRUE,
                                       include_type    = FALSE,
                                       include_alltrees = FALSE,
                                       drop_missing_location = TRUE) {

  library(dplyr)
  library(stringr)
  library(ggplot2)

  fill <- match.arg(fill)

  # --- Host selection ---
  host_set <- get_host_set(
    target_host = host,
    host_meta = host_meta,
    include_species = include_species,
    include_genus = include_genus,
    include_type = include_type,
    include_alltrees = include_alltrees
  )

  d <- df %>%
    mutate(
      Host = str_to_lower(str_squish(Host)),
      State = str_squish(State)
    ) %>%
    filter(Host %in% host_set)

  if (drop_missing_location) {
    d <- d %>% filter(!is.na(State), State != "")
  }

  # --- NE only ---
  northeast_states <- c(
    "Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut",
    "New York","New Jersey","Pennsylvania","Delaware","Maryland","Ohio","West Virginia"
  )

  d <- d %>% filter(State %in% northeast_states)

  # --- Derive clean Agent_group (for ordering) ---
  d <- d %>%
    mutate(
      Agent_group = case_when(

        str_detect(str_to_lower(Agent_Type %||% ""), "abiotic") ~ "Abiotic",
        str_detect(str_to_lower(Agent_Type %||% ""), "disease") ~ "Disease",
        str_detect(str_to_lower(Agent_Type %||% ""), "insect")  ~ "Insect",
        str_detect(str_to_lower(Agent_std), "decline complex") ~ "Decline complex",
        TRUE ~ "Other"
      )
    )

  # --- Collapse to Year × Agent ---
  heat_df <- if (fill == "presence") {
    d %>%
      group_by(Year, Agent_std, Agent_group) %>%
      summarise(value = 1L, .groups = "drop")
  } else {
    d %>%
      group_by(Year, Agent_std, Agent_group) %>%
      summarise(value = n_distinct(State), .groups = "drop")
  }

  # --- Order agents by Agent_group then frequency ---
  agent_order <- heat_df %>%
    group_by(Agent_std, Agent_group) %>%
    summarise(n_years = n_distinct(Year), .groups = "drop") %>%
    arrange(
      factor(Agent_group, levels = c("Insect","Disease","Other","Decline complex","Abiotic")),
      desc(n_years)
    ) %>%
    pull(Agent_std)

  heat_df$Agent_std <- factor(heat_df$Agent_std, levels = agent_order )


  # --- Plot ---
  ggplot(heat_df, aes(x = Year, y = Agent_std, fill = value)) +
    geom_tile() +
    labs(
      title = paste0("Forest pests & conditions over time — ", host),
      subtitle = paste0(
        "Ordered by Agent_Type then frequency | fill = ", fill
      ),
      x = "Year",
      y = "Agent",
      fill = ifelse(fill == "presence", "Present", "# states")
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )+scale_fill_viridis_c()
}

FHM_host_state_enriched %>% select(Agent_std, Agent_Type) %>% unique() %>%
  ungroup()%>%
  group_by(Agent_std) %>%
  summarise(n.rec = n()) %>%
  filter(n.rec > 1)



plot_agent_heatmap_by_host(
  df = FHM_host_state_enriched,
  host = "sugar maple",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = TRUE
)

plot_agent_heatmap_by_host(
  df = FHM_host_state_enriched,
  host = "red maple",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = FALSE
)

plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "oak",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = FALSE
)

plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "eastern white pine",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = FALSE
)

plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "eastern hemlock",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = TRUE
)

plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "ash",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = TRUE
)

plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "american beech",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = TRUE
)


plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "paper birch",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = FALSE
)

plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "yellow birch",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = FALSE
)


plot_agent_heatmap_by_host(
  FHM_host_state_enriched,
  host = "birch",
  host_meta = host_meta,
  fill = "states",
  include_species = TRUE,
  include_genus = TRUE,
  include_type = FALSE,
  include_alltrees = FALSE
)

# Goal: 
# Plot of some common pests and diseases to our species by the type:


# get all of the common species and in NE states:
Region.pest.species <- FHM_host_state_enriched %>% filter(State %in% northeast_states) %>% 
  filter(Host %in% c(species.of.interest, "maple", "oak", 
                     "birch", "poplar", "fir", "spruce", "all trees", "hardwoods", "conifers")) 
  

Region.pest.species %>% select(Year, Agent_Type, Agent_Name, Host, State) %>%
  group_by(Year, Agent_Type, Agent_Name, Host) %>% 
  summarise(nstates = n()) %>%
  arrange(desc(nstates)) %>% View()


# Type == Abiotic:

abiotic.plot.region <- Region.pest.species %>% filter(Agent_Type %in% "Abiotic") %>% 
  # code the drought-delayed mortality after 1988 drought as something different
  mutate(Agent_std = ifelse(Agent_std %in% "Drought" & Year > 1988, "post-Drought mortality", Agent_std)) %>% 
  mutate(Agent_std = ifelse(Agent_std %in% c("Winter injury", "Frost Injury", "Ice Storm Damage") , "Winter/Frost/Ice Storm Injury", Agent_std)) %>% 
  
  group_by(Year,  Agent_std, State) %>%
  summarise(nspecies = n()) %>% 
  ungroup() %>% 
  group_by(Year, Agent_std) %>% 
  summarise(nstates = n(), 
            Nspecies = mean(nspecies))%>%
  ungroup() %>% 
  mutate(Agent_reorder = fct_relevel(Agent_std, 
                                     rev(c("Drought", "post-Drought mortality", 
                                       "Pollution", "Winter/Frost/Ice Storm Injury", "Hurricane Bob"))))|>
  ggplot()+geom_bar(aes(x = Year, y = nstates, fill = Agent_reorder), stat = "identity", color = "grey")+
  theme_bw(base_size  = 12)+
  ylab("Number of states impacted")+
  xlab("Year")+
  xlim(1977, 1992)+
  scale_fill_manual(values = c("Drought"= "#d7191c",
                               "post-Drought mortality"="#fdae61",
                               "Pollution"="#ffffbf",
    "Hurricane Bob"="#abd9e9",
    "Winter/Frost/Ice Storm Injury"="#2c7bb6"), name = " Abiotic Agent")+
  scale_x_continuous(
    breaks = min(FHM_host_state_enriched$Year):max(FHM_host_state_enriched$Year) # Sets breaks for every year
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(filename = paste0(output.dir, "images/abiotic_impacts_by_year_NE.png"), abiotic.plot.region, 
       height = 4, width = 8)



# Type == Disease:
disease.region <- Region.pest.species %>% filter(Agent_Type %in% "Disease") %>% 
  ungroup()%>%
  left_join(., pest.type.df) %>%
  filter(!Class %in% c("defoliating insect", "parasite"))
  
# Decline complexes:---
  
unique.declines <- disease.region %>% 
  filter(Class %in% "decline complex") %>% 
  mutate(Agent_std_decline = ifelse(Agent_std %in% c("Spruce decline complex", "Stillwells Syndrome"), "Spruce Fir decline complex", Agent_std))%>%
  mutate(Host = ifelse(Host %in% "fir", "balsam fir", Host))%>%
  mutate(Host = ifelse(Host %in% "spruce", "red spruce", Host))%>%
  mutate(Host = ifelse(Host %in% c("red", "red oak", "oak", "white oak"), "oak", Host))%>%
  mutate(Host = ifelse(Host %in% c("birch"), "paper birch", Host)) %>%
  mutate(Host = ifelse(Host %in% c("maple"), "sugar maple", Host)) %>%
   ungroup() %>%
 filter(mortality == TRUE | mortality_row == TRUE | increasing == TRUE) %>%
 select(Year, Class, Host, Agent_std_decline, State) %>%
  group_by(Year, Class, Host, Agent_std_decline) %>% 
  summarise(nstate = n())


decline.plot <- unique.declines |>

  ggplot()+#geom_point(aes(x = Year, y = Host, color = Host, size = nstate),  stat = "identity")+
  geom_bar(aes(x = Year, y = nstate, fill = Host),  stat = "identity", color = "grey")+
  
  
  theme_bw(base_size  = 12)+
  ylab("Number of states impacted")+
  xlab("Year")+
  xlim(1977, 1992)+
  scale_fill_manual(values = c("white ash"= "#bababa",
                               "balsam fir"= "#b2df8a",
                               "red spruce"="#b2182b",
                               "sugar maple" =  "#a6cee3",
                              
                               "red maple" = "#e31a1c",
                               "oak"= "#54278f",
                                "red oak"="#8073ac", 
                               "eastern white pine"="#33a02c", 
                               "yellow birch" =   "#fdbf6f",
                               "paper birch"="#fccde5", 
                               "hardwoods" = "#7f3b08"
                               ), name = "Decline complexes")+
  scale_x_continuous(
    breaks = min(FHM_host_state_enriched$Year):max(FHM_host_state_enriched$Year) # Sets breaks for every year
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(filename = paste0(output.dir, "images/decline_complex_impacts_by_year_NE.png"), 
       decline.plot, 
       height = 4, width = 8)

# primary pathogens causing mortality-----
unique.pathogens <- disease.region %>% 
  filter(Class %in% c("fungal pathogen" ,"pathogen")) %>% 
  mutate(Agent_std_decline = ifelse(Agent_std %in% c("Spruce decline complex", "Stillwells Syndrome"), "Spruce Fir decline complex", Agent_std))%>%
  mutate(Host = ifelse(Host %in% "fir", "balsam fir", Host))%>%
  mutate(Host = ifelse(Host %in% "spruce", "red spruce", Host))%>%
  ungroup() %>%
  filter(mortality == TRUE | mortality_row == TRUE | increasing == TRUE) %>%
  #filter(detected == TRUE ) %>%
  select(Year, Class, Host, Agent_std_decline, State, Impact) %>%
  group_by(Year, Class, Host, Agent_std_decline, Impact) %>% 
  summarise(nstate = n())
  
 


pathogen.plot.primary <- unique.pathogens %>% filter(Impact %in% "primary disease") %>% 
  mutate(Host_spp = ifelse(Host %in% c("oak", "red oak"), "northern red oak", 
                           ifelse(Host %in% c("ash", "white ash"), "white ash", Host)))|>

  ggplot()+#geom_point(aes(x = Year, y = Host, color = Host, size = nstate),  stat = "identity")+
  geom_bar(aes(x = Year, y = nstate, fill = Agent_std_decline),  stat = "identity", color = "grey")+
  
  
  theme_bw(base_size  = 12)+
  ylab("Number of states impacted")+
  xlab("Year")+
  xlim(1977, 1992)+
  scale_fill_manual(values = c("Ash Yellows"= "#bababa",
                               "Beech Bark Disease"="#1f78b4",
                               "Oak Wilt"="#8073ac",
                               "White Pine Blister Rust"="#33a02c"), name = "Mortality causeing disease Agent")+
  scale_x_continuous(
    breaks = min(FHM_host_state_enriched$Year):max(FHM_host_state_enriched$Year) # Sets breaks for every year
  ) +
  theme(panel.grid.minor = element_blank())




ggsave(filename = paste0(output.dir, "images/primary_disease_impacts_by_year_NE.png"), 
       pathogen.plot.primary, 
       height = 4, width = 8)

# secondary pathogens, by host species:
secondary.pathogens <- disease.region %>% 
  filter(Class %in% c("fungal pathogen" ,"pathogen")) %>% 
  mutate(Agent_std_decline = ifelse(Agent_std %in% c("Spruce decline complex", "Stillwells Syndrome"), "Spruce Fir decline complex", Agent_std))%>%
  mutate(Host = ifelse(Host %in% "fir", "balsam fir", Host))%>%
  mutate(Host = ifelse(Host %in% "spruce", "red spruce", Host))%>%
  ungroup() %>%
  #filter(mortality == TRUE | mortality_row == TRUE | increasing == TRUE) %>%
  #filter(detected == TRUE ) %>%
  select(Year, Class, Host, Agent_std_decline, State, Impact) %>%
  group_by(Year, Class, Host, Agent_std_decline, Impact) %>% 
  summarise(nstate = n()) %>% 
  filter(Impact %in% "secondary stressor") 


type.of.secondary <-  tibble::tribble(
  ~Agent_std_decline, ~type_secondary,
  "Annosus Root and Butt Rot", "Root rot",
  "Shoestring Root rot",   "Root rot",
  "Armillaria Root Rot", "Root rot",
  "Cytospora Canker", "Canker",        
   "Yellow Birch Canker",  "Canker", 
  "Scleroderris Canker","Canker",
  "Caliciopsis Canker" ,  "Canker",    
  "Hypoxylon Canker","Canker",
  "Maple Canker",   "Canker",
  
   "Balsam Fir Needle Rust" , "needle/leaf disease",  
   
  "Anthracnose" , "needle/leaf disease", 
  
  "Leaf Spot" ,   "needle/leaf disease",             
  "Pine Wood Nematode","Other",
  "Cenangium Twig Blight","Other",
  "Stegansporium", "Canker"
   
  )

pathogen.plot.secondary <- secondary.pathogens%>% 
  mutate(Host = ifelse(Host %in% c("oak", "red oak"), "oak", 
                           ifelse(Host %in% c("ash", "white ash"), "white ash", Host)))%>%
  mutate(Host = ifelse(Host %in% c("hickory"), "hardwoods", Host)) %>% 
  mutate(Host = ifelse(Host %in% c("sugar maple", "red maple", "maple"), "maple", Host)) %>%
  left_join(., type.of.secondary)|>
  
  ggplot()+#geom_point(aes(x = Year, y = Host, color = Host, size = nstate),  stat = "identity")+
  geom_bar(aes(x = Year, fill = Host),  stat = "count")+
  
  
  theme_bw(base_size  = 12)+
  ylab("Number of secondary pathogens reported")+
  xlab("Year")+
  xlim(1977, 1992)+

  scale_fill_manual(values = c("conifers"=  "#1a9850",
                               "balsam fir"= "#b2df8a",
                               #"red spruce"="#b2182b",
                               "sugar maple" =  "#a6cee3",
                               "maple" =  "#a6ced4",
                               "red maple" = "#e31a1c",
                               "oak"= "#54278f",
                               "red oak"="#8073ac", 
                               "eastern white pine"="#33a02c", 
                               "yellow birch" =   "#fdbf6f",
                               "paper birch"="#fccde5", 
                               "hardwoods" = "#7f3b08"
  ), name = "Decline complexes")+
  scale_x_continuous(
    breaks = min(FHM_host_state_enriched$Year):max(FHM_host_state_enriched$Year) # Sets breaks for every year
  ) +
  theme(panel.grid.minor = element_blank())

pathogen.plot.secondary

ggsave(filename = paste0(output.dir, "images/secondary_disease_impacts_by_year_NE.png"), 
       pathogen.plot.secondary, 
       height = 4, width = 8)

# Insects reported----
insects.region <- Region.pest.species %>% filter(Agent_Type %in% "Insect") %>% 
  ungroup()%>%
  left_join(., pest.type.df)

# conifer insects:
insects.conifers <- Region.pest.species %>% 
  filter(Agent_Type %in% "Insect") %>% 
  ungroup()%>%
  left_join(., pest.type.df) %>% 
  filter(Host %in% c("spruce", "eastern hemlock", "fir", "balsam fir", "red spruce", "eastern white pine", "conifers", "cedar"))%>%
  arrange(Agent_std, Year) %>% select(Year, Agent_Type, Agent_std, Host, Remarks, State) %>% 
  group_by(Year, Agent_Type, Agent_std, Host, Remarks) %>% 
  summarise(nstates = n())

unique(insects.conifers$Agent_std)
write.csv(insects.conifers, "data/conifer_insects.csv")

# read through these and selected several for plotting:
conifer.insects.c <- read.csv("data/conifer_insects_checked.csv")


conifer.insects.c %>% filter(Include == TRUE) |> ggplot()+
  geom_bar(aes(x = Year, fill = Agent_std))

conifer.defoliators.species <- conifer.insects.c %>% filter(Include == TRUE) %>%
  select(Start, End, Agent_std, Host, Type_host)%>% 
  distinct()%>%
  mutate(Host = fct_relevel(Host, rev(c("balsam fir", "red spruce", "eastern hemlock", "eastern white pine"))))%>%
  mutate(Agent_std = fct_relevel(Agent_std, 
                                 rev(c("Spruce Budworm", "Spruce Beetle", "White Pine Weevil", "Hemlock Looper", "Hemlock Scale", "Hemlock Woolly Adelgid",
                                      "Conifer Swift Moth", "Conifer Sawflies", "Balsam Woolly Adelgid", "Balsam Twig Aphid", "Balsam Gall Midge"))))|> 
  ggplot()+
  geom_segment(aes(x = Start, xend = End, y = Agent_std, color = Host, linetype  = Type_host),
               position = position_dodge(width = 0.6), 
               #linejoin = "round", 
               lineend = "round", size = 3)+
  theme_bw(base_size = 12)+
  scale_color_manual(values = c(
   "balsam fir"= "#b2df8a", # balsam fir
   "eastern hemlock"="#003f21", # hemlock
   "red spruce"= "#b2182b", # red spruce
    #"#fee090", # n white cedar
   "eastern white pine"= "#33a02c" # eastern white pine
    
  ))+
  #scale_alpha_manual(values = c("secondary" = 0.65, "primary" = 1))+
  scale_linetype_manual(values = c("secondary" = "11", "primary" = "solid"))+
  ylab("Conifer Insect Defoliators")+xlab("Year")+
  scale_x_continuous(
    breaks = min(FHM_host_state_enriched$Year):max(FHM_host_state_enriched$Year) # Sets breaks for every year
  ) 
conifer.defoliators.species

ggsave(filename = paste0(output.dir, "images/conifer_defoliators_impacts_by_year_NE.png"), 
       conifer.defoliators.species, 
       height = 5, width = 8)

# get other defoliators impacting our species of interest:
hardwood.defoliators <- Region.pest.species %>% 
  filter(Agent_Type %in% "Insect") %>% 
  ungroup()%>%
  left_join(., pest.type.df) %>% 
  filter(Host %in% c("maple", "oak", "hardwoods", "poplar", "paper birch", "white oak", "birch", "black cherry", "hickory", "red oak", "sugar maple", "red maple", "american beech", "yellow poplar", "yellow birch"))%>%
  arrange(Agent_std, Year) %>% select(Year, Agent_Type, Agent_std, Host, Remarks, State) %>% 
  group_by(Year, Agent_Type, Agent_std, Host, Remarks) %>% 
  summarise(nstates = n())
unique(hardwood.defoliators$Agent_std)
write.csv(hardwood.defoliators, "data/hardwood_insects.csv")


insects.region %>%
  #filter(increasing == TRUE | mortality == TRUE)%>%
  group_by(Year, Agent_std, Host, Impact) %>%
  filter(Impact %in% "defoliation")%>% 
  summarise(nstate = n()) |>

ggplot()+#geom_point(aes(x = Year, y = Host, color = Host, size = nstate),  stat = "identity")+
  geom_bar(aes(x = Year, fill = Agent_std),  stat = "count", color = "grey")+
  #geom_tile(aes(x = Year, y = Host, fill = nstate), color = "grey")+
  
  
  theme_bw(base_size  = 12)+
  #ylab("Number of states impacted")+
  #xlab("Year")+
  #xlim(1977, 1992)+
  # scale_fill_manual(values = c("Ash Yellows"= "#bababa",
  #                              "Beech Bark Disease"="#1f78b4",
  #                              "Oak Wilt"="#8073ac",
  #                              "White Pine Blister Rust"="#33a02c"), name = "Mortality causeing disease Agent")+
 # scale_x_continuous(
#    breaks = min(FHM_host_state_enriched$Year):max(FHM_host_state_enriched$Year) # Sets breaks for every year
 # ) +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90,hjust = 1 ))+facet_wrap(~Host)

# for each species get 

# Next steps:

# 1. Expand out from 1975-1997
# 2. arrange heatmaps by defoliators, chronic diseases, etc
# 3.
# 4. Save full dataset for each pest & evaluate increasing vs decreasing

