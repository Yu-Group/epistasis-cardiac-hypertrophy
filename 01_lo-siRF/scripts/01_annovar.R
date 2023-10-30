rm(list = ls())
library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)

source(file.path("..", "functions", "utils.R"), chdir = TRUE)

# references:
# https://www.siweizhang.com/2018/10/30/How-to-do-Annotation/
# bim format: https://www.cog-genomics.org/plink/1.9/formats#bim
# avinput format: https://annovar.openbioinformatics.org/en/latest/user-guide/input/

# set file paths
BIM_DATA_PATH <- file.path(
  "..", "data", paste0("ukbb_wbr_imp_morgan_chr", 1:22)
)
ANNOVAR_DIR <- file.path("..", "software", "annovar")
OUT_PATH <- file.path("..", "results", "annovar", "ukbb_wbr_imp_morgan_chr")

if (!dir.exists(dirname(OUT_PATH))) {
  dir.create(dirname(OUT_PATH), recursive = TRUE)
}

# get the reference genome
genome <- BSgenome.Hsapiens.UCSC.hg19
chrs <- seqnames(Hsapiens)[1:22]
# seqlengths(genome)
# seqnames(genome)
# Hsapiens$chr2

for (chr in 1:22) {
  print(chr)

  # format plink data for annovar input
  plink.dat <- read.table(paste0(BIM_DATA_PATH, chr, ".bim"),
    header = FALSE, sep = "\t", stringsAsFactors = FALSE
  )
  ref.seq <- getSeq(Hsapiens,
    names = chrs[chr],
    start = 1, end = max(plink.dat[, 4])
  )
  av.dat <- plink.dat %>%
    select(CHR = V1, Start = V4, End = V4, Minor = V5, Major = V6) %>%
    mutate(
      Ref = strsplit(as.character(ref.seq[Start]), split = "")[[1]],
      Alt = ifelse(Ref == Major, Minor,
        ifelse(Ref == Minor, Major, NA)
      )
    )
  av.dat[is.na(av.dat$Alt), c("Ref", "Alt")] <- av.dat[
    is.na(av.dat$Alt),
    c("Major", "Minor")
  ]
  write.table(av.dat %>%
    select(CHR, Start, End, Ref, Alt),
  paste0(OUT_PATH, chr, ".avinput"),
  row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
  )

  # run annovar
  av <- file.path(ANNOVAR_DIR, "table_annovar.pl")
  db <- file.path(ANNOVAR_DIR, "humandb")
  system(paste0(
    "perl ", av, " ", paste0(OUT_PATH, chr, ".avinput"), " ", db,
    " -buildver hg19",
    " -outfile ", paste0(OUT_PATH, chr, "_avresults"),
    " -remove ",
    " -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a",
    " -operation g,r,f,f,f -nastring . -csvout -polish"
  ))
}

# merge results into a single data frame
snps_ls <- list()
for (chr in 1:22) {
  av_out <- fread(paste0(OUT_PATH, chr, "_avresults.hg19_multianno.csv"),
    header = TRUE
  )
  snps <- fread(paste0(BIM_DATA_PATH, chr, ".bim"))
  snps_ls[[chr]] <- cbind(
    snps[, -3],
    av_out %>% select(Gene.refGene, Func.refGene, ExonicFunc.refGene)
  ) %>%
    setNames(c(
      "Chr", "rsID", "Pos", "Alt", "Ref",
      "Gene", "Gene Function", "Exonic Function"
    ))
}

snps_df <- map_dfr(snps_ls, ~.x) %>%
  distinct(rsID, .keep_all = TRUE) %>%
  mutate(Name = renameSNP(rsID, Chr)) %>%
  select(Name, everything())
# saveRDS(snps_df, paste0(OUT_PATH, "snp2gene_df.rds"))
write.csv(snps_df, file.path(dirname(OUT_PATH), "snp2gene_df.csv"),
  row.names = FALSE
)
