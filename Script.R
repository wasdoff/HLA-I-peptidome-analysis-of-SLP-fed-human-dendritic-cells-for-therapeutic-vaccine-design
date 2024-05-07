library(tidyverse)
library(glue)
library(conflicted)
conflicts_prefer(dplyr::filter)

# Prepare input data ----

## Read data ----

psm_ <- read_csv(
  "DIA search 2023-06-28/dia_db.psms.csv",
  guess_max = 1e6,
  na = c("", "NA", "0"),
  show_col_types = FALSE
) |>
  rownames_to_column("Index") |>
  glimpse()

read_fasta <- function(path) {
  tibble(
    fasta_text = read_lines(path),
    sequence_block = cumsum(str_detect(fasta_text, ">"))
  ) |>
    group_by(sequence_block) |>
    summarise(
      Header = fasta_text[1] |> str_extract("(?<=^...\\|)[^ ]+"),
      Seq = fasta_text[2:n()] |> paste(collapse = "")
    ) |>
    select(-sequence_block) |>
    deframe()
}

db <- read_fasta("homo_sapiens_02-Mar-2022_incl_additional_12SLP_4-HepB-RefProteomes_validated.fasta")

psm_2 <- psm_ |>
  select(Peptide, Accession) |>
  mutate(
    Stripped = str_remove_all(Peptide, "\\([^\\)]+\\)"),
    ManualAccession = local({
      missing_acc <- which(is.na(Accession))
      out <- rep(NA_character_, times = n())
      out[missing_acc] <- map(Stripped[missing_acc], function(x) {
        
        x_il <- x |> str_replace_all("[IL]", "L")
        db_il <- db |> str_replace_all("[IL]", "L")
        match <- str_detect(db_il, fixed(x_il))
        names(db[match])
        
      }, .progress = TRUE)
      out
    })
  ) |>
  glimpse()

psm_3 <- psm_2 |>
  group_by(Stripped) |>
  mutate(
    AccessionIndex = c(ManualAccession[1], str_split(Accession[1], fixed(";"))) |> 
      unlist() |>
      discard(is.na) |>
      map_int(\(x) which(x == names(db))) |>
      list()
  ) |>
  ungroup() |>
  glimpse()

psm_3 |>
  mutate(Missing = is.na(Accession)) |>
  select(Missing, AccessionIndex) |>
  unchop(AccessionIndex) |>
  (\(x) {
    patchwork::wrap_plots(
      ggplot(x) +
        geom_histogram(
          aes(
            fill = Missing,
            x = AccessionIndex
          ),
          alpha = 0.5,
          binwidth = 1000, 
          boundary = 0
        ) +
        theme_bw(),
      ggplot(x) +
        geom_density(
          aes(
            fill = Missing,
            x = AccessionIndex
          ),
          alpha = 0.5
        ) +
        theme_bw()
    )
  })()

psm <- psm_ |> mutate(
  Accession = if_else(
    !is.na(Accession),
    Accession,
    map_chr(psm_2$ManualAccession, \(x) paste(x, collapse = ";"))
  )
)

## Add sample annotation ----

annot <- read_csv(
  "Raw file annotation DIA.csv",
  show_col_types = FALSE
) |>
  print(n = Inf)

annot_split <- local({
  temp <- annot |>
    mutate(Raw = str_split(Raw, "; "), Name = str_remove(Name, "^Area ")) |>
    unnest(Raw)
  temp |>
    select(-Raw) |>
    split(as_factor(temp$Raw))
}) |>
  glimpse()

psm_annot <- psm |>
  group_by(`Source File`) |>
  mutate(
    Annotation = rlang::inject(switch(`Source File`[1], !!!annot_split)),
    .after = `Source File`
  ) |>
  ungroup() |>
  unpack(Annotation) |>
  filter(Relevant == "Yes") |>
  select(-Relevant) |>
  glimpse()

## Deal with missing accessions by overwriting them with my own fasta-matched
psm_annot <- psm_annot |>
  mutate(
    Sequence = str_remove_all(Peptide, "\\([^\\)]*\\)"),
    Old.Accession = Accession
  ) |>
  mutate(
    .by = Sequence,
    Accession = Sequence[1] |> 
      fixed() |> 
      str_which(db, pattern = _) |> 
      (\(x) names(db[x]))() |> 
      paste(collapse = ";")
  ) |>
  glimpse()

## Dealing with equivalent PSMs and summing precursors ----
# Two cases of duplicates:
#
# 1. Same precusor (sequence, modifcations and charge) reported more than once
#   ==> we pick the PSM with the highest Area (or first if all Areas are NA)
#
# 2. PSM is reported twice but with a different sequence because I/L
#   ==> we keep these, but flag them in the next section
#
# After doing this we will sum all precursors to their peptide sequence, carrying over the I/L flag
is_first_max <- function(x) {
  # x is a numeric vector
  # Returns a length(x) sized lgl with 1 TRUE for first occurence of max(x)
  # If all of x is NA, then the first element is TRUE
  if (all(is.na(x))) {
    lgl <- rep(TRUE, times = length(x))
    lgl[-1] <- NA
    return(lgl)
  }
  lgl <- x == max(x, na.rm = TRUE)
  first_true  <- which(lgl)[1]
  lgl[-first_true] <- FALSE
  return(lgl)
}

slp <- c(
  SLP1  = "HYFQTRHYLHTLWKAGILYKRETTR",
  SLP2  = "TSFPWLLGCAANWILRGTSFVYVPS",
  SLP3  = "SVVRRAFPHCLAFSYMDDVVLGAKS",
  SLP4  = "LSAMSTTDLEAYFKDCLFKDWEELG",
  SLP5  = "HLSLRGLPVCAFSSAGPCALRFTSA",
  SLP6  = "HHIRIPRTPARVTGGVFLVDKNPHN",
  SLP7  = "AARLCCQLDPARDVLCLRPVGAESR",
  SLP8  = "RKLHLYSHPIILGFRKIPMGVGLSP",
  SLP9  = "ARQRPGLCQVFADATPTGWGLAIGH",
  SLP10 = "SPSVPSHLPDRVHFASPLHVAWRPP",
  SLP11 = "ASSSSSCLHQSAVRKAAYSHLSTSK",
  SLP12 = "GFAAPFTQCGYPALMPLYACIQAKQA"
)

slp_I2L <- str_replace_all(slp, "I", "L")

input <- psm_annot |>
  ### Remove some columns ----
  select(
    Index,
    Peptide,
    Sequence,
    Length,
    z,
    Area,
    `Source File`:Accession
  ) |>
  ### Deal with case 1 duplicates ----
  mutate(
    .by = c(Peptide, `Source File`, z),
    first_max = is_first_max(Area),
    Index = glue(
      "{main_index}={secondary_indices}",
      main_index = Index[first_max],
      secondary_indices = Index[!first_max] |> paste(collapse = ";")
    ),
    Index = ifelse(
      str_detect(Index, "=$"),
      str_remove(Index, "=$"),
      glue("<{Index}>")
    )
  ) |>
  filter(first_max) |>
  select(-first_max) |>
  ### Sum precursors ----
  summarise(
    .by = c(Sequence, `Source File`, Length, `Source File`:Accession),
    Index = paste(Index, collapse = "+"),
    Index = ifelse(
      str_detect(Index, "\\+"),
      glue("[{Index}]"),
      Index
    ),
    `Summed Area` = sum(Area, na.rm = TRUE)
  ) |>
  ### Calculate geometric mean between replicates ----
  summarise(
    .by = c(Sequence, Name:Accession),
    Index = paste(Index, collapse = "&"),
    `Geo. mean Summed Area` = exp(mean(log(`Summed Area`[`Summed Area` > 0]), na.rm = TRUE)),
    `Summed Area` = set_names(`Summed Area`, `Source File`) |> list()
  ) |>
  relocate(Index) |>
  ### Adding I2L column ----
  mutate(
    I2L = str_replace_all(Sequence, "I", "L"),
    SLP = case_when(
      str_detect(slp_I2L[01], I2L) ~ "SLP1",
      str_detect(slp_I2L[02], I2L) ~ "SLP2",
      str_detect(slp_I2L[03], I2L) ~ "SLP3",
      str_detect(slp_I2L[04], I2L) ~ "SLP4",
      str_detect(slp_I2L[05], I2L) ~ "SLP5",
      str_detect(slp_I2L[06], I2L) ~ "SLP6",
      str_detect(slp_I2L[07], I2L) ~ "SLP7",
      str_detect(slp_I2L[08], I2L) ~ "SLP8",
      str_detect(slp_I2L[09], I2L) ~ "SLP9",
      str_detect(slp_I2L[10], I2L) ~ "SLP10",
      str_detect(slp_I2L[11], I2L) ~ "SLP11",
      str_detect(slp_I2L[12], I2L) ~ "SLP12",
      .default = "Not SLP"
    ),
    .after = Sequence
  )

table(input$SLP)

# Inventorize and flag I/L duplicates ----

## Finding the I/L duplicates ----

duplicates_overall <- input |>
  select(Sequence, I2L) |>
  distinct() |>
  chop(Sequence) |>
  mutate(
    N = lengths(Sequence),
    Sequence = map_chr(Sequence, \(x) paste(x, collapse = ";"))
  ) |>
  print()

duplicates_perAdj <- input |>
  select(Sequence, I2L, Adjuvant) |>
  distinct() |>
  chop(Sequence) |>
  mutate(
    N = lengths(Sequence),
    Sequence = map_chr(Sequence, \(x) paste(x, collapse = ";"))
  ) |>
  print()

duplicates_perName <- input |>
  select(Sequence, I2L, Name, Adjuvant) |>
  distinct() |>
  chop(Sequence) |>
  mutate(
    N = lengths(Sequence),
    Sequence = map_chr(Sequence, \(x) paste(x, collapse = ";"))
  ) |>
  print()

## Flagging the I/L duplicates ----

input <- input |>
  mutate(
    `Potential IL duplicate` = I2L %in% duplicates_overall$I2L[duplicates_overall$N > 1],
    .after = I2L
  )

# Looking at length distributions ----

# NOTE: two version (1) keep I/L duplicates (2) drop

lengthdist_overall <- input |>
  select(Sequence) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(.by = Length, N = n()) |>
  arrange(Length) |>
  print()

lengthdist_overall_noILdups <- input |>
  filter(!`Potential IL duplicate`) |>
  select(Sequence) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(.by = Length, N = n()) |>
  arrange(Length) |>
  print()

lengthdist_perAdj <- input |>
  select(Sequence, Adjuvant) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(.by = c(Length, Adjuvant), N = n()) |>
  arrange(Length) |>
  print()

lengthdist_perAdj_noILdups <- input |>
  filter(!`Potential IL duplicate`) |>
  select(Sequence, Adjuvant) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(.by = c(Length, Adjuvant), N = n()) |>
  arrange(Length) |>
  print()

lengthdist_perName <- input |>
  select(Sequence, Name, Adjuvant) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(.by = c(Length, Name, Adjuvant), N = n()) |>
  arrange(Length) |>
  print()

lengthdist_perName_noILdups <- input |>
  filter(!`Potential IL duplicate`) |>
  select(Sequence, Name, Adjuvant) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(.by = c(Length, Name, Adjuvant), N = n()) |>
  arrange(Length) |>
  print()

# AA dist ----

# NOTE: only while keeping I/L duplicates
# Recreate previous stacked bar graph using adjusted color palette (see email somewhere)

aa_levels <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

my_palette <- c(
  rgb(041, 097, 117, maxColorValue = 255),
  rgb(068, 160, 194, maxColorValue = 255),
  rgb(091, 127, 140, maxColorValue = 255),
  rgb(158, 196, 204, maxColorValue = 255),
  rgb(023, 055, 066, maxColorValue = 255),
  rgb(044, 080, 117, maxColorValue = 255),
  rgb(074, 134, 194, maxColorValue = 255),
  rgb(085, 105, 125, maxColorValue = 255),
  rgb(171, 184, 196, maxColorValue = 255),
  rgb(025, 046, 066, maxColorValue = 255),
  rgb(036, 093, 094, maxColorValue = 255),
  rgb(068, 200, 175, maxColorValue = 255),
  rgb(067, 142, 143, maxColorValue = 255),
  rgb(161, 191, 191, maxColorValue = 255),
  rgb(017, 045, 044, maxColorValue = 255),
  rgb(049, 109, 117, maxColorValue = 255),
  rgb(081, 181, 194, maxColorValue = 255),
  rgb(099, 138, 143, maxColorValue = 255),
  rgb(192, 204, 206, maxColorValue = 255),
  rgb(028, 062, 066, maxColorValue = 255)
)

uniq_pep <- input |>
  select(Sequence, Name) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  filter(Length |> between(8, 11))

aa_freq_counts <- uniq_pep |>
  split(uniq_pep$Length) |>
  map(function(tbl) {
    tbl |>
      mutate(aa = str_split(Sequence, "") |> map(function(y) factor(y, aa_levels))) |>
      unnest_wider(aa, names_sep = "_") |>
      reframe(
        across(contains("aa_"),  function(column) {
          table(column) |> list()
        }),
        .by = c(Name, Length)
      ) |>
      mutate(
        AA = list(aa_levels),
        .before = aa_1,
        .by = Name
      ) |>
      unnest(AA) |>
      rowwise() |>
      mutate(
        across(contains("aa_"), function(column) pluck(column, AA))
      ) |>
      ungroup() |>
      pivot_longer(
        contains("aa_"),
        names_to = "Position",
        values_to = "Count"
      )
  }) |>
  bind_rows() |>
  mutate(Position = Position |> str_remove("aa_") |> as.integer()) |>
  print()

aa_freq_plots <- aa_freq_counts |>
  split(aa_freq_counts$Name) |>
  imap(function(tbl, name) {
    tbl |>
      mutate(
        Fraction = Count / sum(Count),
        .by = c(Length, Position)
      ) |>
      ggplot() +
      geom_col(
        aes(
          x = Position,
          y = Fraction,
          fill = AA
        ),
        width = 0.95
      ) +
      scale_x_continuous(breaks = 1:100, expand = expansion()) +
      scale_y_continuous(labels = scales::label_percent(), expand = expansion()) +
      scale_fill_discrete(type = my_palette) +
      facet_wrap(vars(Length), scales = "free_x") +
      ggtitle(name) +
      theme_bw() +
      theme(
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "cm")
      )
  })

# Overlappings IDs between Adjuvants ----

adjuvant_id_overlaps <- input |>
  select(Sequence, Adjuvant) |>
  distinct() |>
  mutate(Length = nchar(Sequence)) |>
  summarise(
    .by = c(Sequence, Length),
    Adjuvant = sort(Adjuvant) |> paste(collapse = "&")
  ) |>
  summarise(
    .by = c(Adjuvant, Length),
    N = n()
  )

adjuvant_id_overlaps <- adjuvant_id_overlaps |>
  mutate(Length = str_pad(Length, width = 2, pad = 0)) |>
  bind_rows(tibble(
    Adjuvant = c("AV&PIC", "AV&PIC", "AV", "AV", "PIC", "PIC"),
    Length = c("total", "total 8--11") |> rep(3),
    N = c(
      adjuvant_id_overlaps |> filter(Adjuvant == "AV&PIC") |> pull(N) |> sum(),
      adjuvant_id_overlaps |> filter(Adjuvant == "AV&PIC", Length |> between(8, 11)) |> pull(N) |> sum(),
      adjuvant_id_overlaps |> filter(Adjuvant == "AV") |> pull(N) |> sum(),
      adjuvant_id_overlaps |> filter(Adjuvant == "AV", Length |> between(8, 11)) |> pull(N) |> sum(),
      adjuvant_id_overlaps |> filter(Adjuvant == "PIC") |> pull(N) |> sum(),
      adjuvant_id_overlaps |> filter(Adjuvant == "PIC", Length |> between(8, 11)) |> pull(N) |> sum()
    )
  )) |>
  arrange(str_width(Adjuvant)) |>
  pivot_wider(names_from = Adjuvant, values_from = N) |>
  arrange(Length) |>
  mutate(Total = AV + PIC + `AV&PIC`) |>
  print()

# Count peptides per protein / gene ----

# Variants of analysis:
# 1 Use first reported acc
# 2 Use first reported gene
# 3 Use most reported acc in overall data
# 4 Use most reported gene in overall data
# All: overall + split by adjuvant

# 1
pep_per_prot_first_acc <- input |>
  select(Adjuvant, Name, Sequence, Accession) |>
  mutate(Accession = str_split_i(Accession, ";", 1) |> str_extract("^[^\\|]+")) |>
  distinct() |> 
  mutate(
    .by = Sequence,
    `Count of Sequence` = n()
  ) |>
  mutate(
    .by = Accession,
    `Count of Accession` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Sequence),
    `Count of same Adjuvant and Sequence` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Accession),
    `Count of same Adjuvant and Accession` = n()
  ) |>
  mutate(
    .by = Sequence,
    `Sequence is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Accession,
    `Accession is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  glimpse()

# 2
pep_per_prot_first_gene <- input |>
  select(Adjuvant, Name, Sequence, Accession) |>
  mutate(
    Gene = str_split_i(Accession, ";", 1) |> str_extract("[^\\|]+$"),
    Accession = NULL
  ) |>
  distinct() |>
  mutate(
    .by = Sequence,
    `Count of Sequence` = n()
  ) |>
  mutate(
    .by = Gene,
    `Count of Gene` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Sequence),
    `Count of same Adjuvant and Sequence` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Gene),
    `Count of same Adjuvant and Gene` = n()
  ) |>
  mutate(
    .by = Sequence,
    `Sequence is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Gene,
    `Gene is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  glimpse()

# 1+2
pep_per_prot_first_accgene <- input |>
  select(Adjuvant, Name, Sequence, Accession) |>
  mutate(
    Gene = str_split_i(Accession, ";", 1) |> str_extract("[^\\|]+$"),
    Accession = str_split_i(Accession, ";", 1) |> str_extract("^[^\\|]+")
  ) |>
  distinct() |>
  mutate(
    .by = Sequence,
    `Count of Sequence` = n()
  ) |>
  mutate(
    .by = Gene,
    `Count of Gene` = n()
  ) |>
  mutate(
    .by = Accession,
    `Count of Accession` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Sequence),
    `Count of same Adjuvant and Sequence` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Gene),
    `Count of same Adjuvant and Gene` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Accession),
    `Count of same Adjuvant and Accession` = n()
  ) |>
  mutate(
    .by = Sequence,
    `Sequence is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Gene,
    `Gene is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Accession,
    `Accession is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  glimpse()

# 3
ranked_acc <- input |>
  select(Sequence, Accession) |>
  distinct() |>
  pull(Accession) |>
  str_split_i(";", 1) |> 
  str_extract("^[^\\|]+") |>
  table() |>
  enframe("Accession", "N") |>
  mutate(
    N = as.integer(N),
    Rank = rank(-N, ties.method = "first")
  ) |>
  arrange(Rank) |>
  glimpse()

pep_per_prot_toprank_acc <- input |>
  select(Adjuvant, Name, Sequence, Accession) |>
  mutate(Accession = str_split(Accession, ";") |> map_chr(\(x) {
    x <- str_extract(x, "^[^\\|]+")
    ranked_acc |> 
      filter(Accession %in% x) |>
      slice(1) |>
      pull(Accession)
  }, .progress = TRUE)) |>
  distinct() |>
  mutate(
    .by = Sequence,
    `Count of Sequence` = n()
  ) |>
  mutate(
    .by = Accession,
    `Count of Accession` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Sequence),
    `Count of same Adjuvant and Sequence` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Accession),
    `Count of same Adjuvant and Accession` = n()
  ) |>
  mutate(
    .by = Sequence,
    `Sequence is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Accession,
    `Accession is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  glimpse()

# 4
ranked_gene <- input |>
  select(Sequence, Accession) |>
  distinct() |>
  pull(Accession) |>
  str_split_i(";", 1) |> 
  str_extract("[^\\|]+$") |>
  table() |>
  enframe("Gene", "N") |>
  mutate(
    N = as.integer(N),
    Rank = rank(-N, ties.method = "first")
  ) |>
  arrange(Rank) |>
  glimpse()

pep_per_prot_toprank_gene <- input |>
  select(Adjuvant, Name, Sequence, Accession) |>
  mutate(
    Gene = str_split(Accession, ";") |> map_chr(\(x) {
      x <- str_extract(x, "[^\\|]+$")
      ranked_gene |> 
        filter(Gene %in% x) |>
        slice(1) |>
        pull(Gene)
    }, .progress = TRUE),
    Accession = NULL
  ) |>
  distinct() |>
  mutate(
    .by = Sequence,
    `Count of Sequence` = n()
  ) |>
  mutate(
    .by = Gene,
    `Count of Gene` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Sequence),
    `Count of same Adjuvant and Sequence` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Gene),
    `Count of same Adjuvant and Gene` = n()
  ) |>
  mutate(
    .by = Sequence,
    `Sequence is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Gene,
    `Gene is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  glimpse()

# 3+4
pep_per_prot_toprank_accgene <- input |>
  select(Adjuvant, Name, Sequence, Accession) |>
  mutate(
    Gene = str_split(Accession, ";") |> map_chr(\(x) {
      x <- str_extract(x, "[^\\|]+$")
      ranked_gene |> 
        filter(Gene %in% x) |>
        slice(1) |>
        pull(Gene)
    }, .progress = TRUE),
    Accession = str_split(Accession, ";") |> map_chr(\(x) {
      x <- str_extract(x, "^[^\\|]+")
      ranked_acc |> 
        filter(Accession %in% x) |>
        slice(1) |>
        pull(Accession)
    }, .progress = TRUE)
  ) |>
  distinct() |>
  mutate(
    .by = Sequence,
    `Count of Sequence` = n()
  ) |>
  mutate(
    .by = Gene,
    `Count of Gene` = n()
  ) |>
  mutate(
    .by = Accession,
    `Count of Accession` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Sequence),
    `Count of same Adjuvant and Sequence` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Gene),
    `Count of same Adjuvant and Gene` = n()
  ) |>
  mutate(
    .by = c(Adjuvant, Accession),
    `Count of same Adjuvant and Accession` = n()
  ) |>
  mutate(
    .by = Sequence,
    `Sequence is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Gene,
    `Gene is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  mutate(
    .by = Accession,
    `Accession is unique to` = {
      adj <- unique(Adjuvant)
      if (length(adj) > 1) NA_character_ else adj
    }
  ) |>
  glimpse()

# Exporting ----

writexl::write_xlsx(
  x = list(
    input = input,
    duplicates_overall = duplicates_overall,
    duplicates_perAdj = duplicates_perAdj,
    duplicates_perName = duplicates_perName,
    lengthdist_overall = lengthdist_overall,
    lengthdist_overall_noILdups = lengthdist_overall_noILdups,
    lengthdist_perAdj = lengthdist_perAdj,
    lengthdist_perAdj_noILdups = lengthdist_perAdj_noILdups,
    lengthdist_perName = lengthdist_perName,
    lengthdist_perName_noILdups = lengthdist_perName_noILdups,
    aa_freq_counts = aa_freq_counts,
    adjuvant_id_overlaps = adjuvant_id_overlaps,
    pep_per_prot_first_acc = pep_per_prot_first_acc,
    pep_per_prot_first_gene = pep_per_prot_first_gene,
    pep_per_prot_first_accgene = pep_per_prot_first_accgene,
    pep_per_prot_toprank_acc = pep_per_prot_toprank_acc,
    pep_per_prot_toprank_gene = pep_per_prot_toprank_gene,
    pep_per_prot_toprank_accgene = pep_per_prot_toprank_accgene
  ),
  path = "R_Export.xlsx"
)

for (i in seq_along(aa_freq_plots)) {
  ggsave(
    filename = paste0(names(aa_freq_plots)[i], ".pdf"),
    plot = aa_freq_plots[[i]]
  )
}
