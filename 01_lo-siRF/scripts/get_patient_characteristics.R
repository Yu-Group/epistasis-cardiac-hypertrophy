library(magrittr)
rm(list = ls())

DATA_DIR <- file.path("..", "data")
RESULTS_DIR <- file.path("..", "results")
TABLES_SUPP_DIR <- file.path(RESULTS_DIR, "tables_supplementary")
CLIN_DATA_PATH <- file.path(DATA_DIR, "ukbb_clinical_data.tab")
CMRI_DATA_PATH <- file.path(DATA_DIR, "table_ventricular_volume_with_indexing.csv")
GENO_IDS_PATH <- file.path(DATA_DIR, "geno_ids.rds")

#### Read in raw data ####
clin_data <- data.table::fread(CLIN_DATA_PATH, nrows = 5)

#### Clean data ####
clin_data_filtered <- clin_data %>%
  dplyr::select(
    f.eid, f.31.0.0, f.34.0.0,
    tidyselect::starts_with("f.41202."),
    tidyselect::starts_with("f.41204."),
    tidyselect::starts_with("f.20002."),
    tidyselect::starts_with("f.6150."),
    tidyselect::starts_with("f.2443."),
    tidyselect::starts_with("f.2473."),
    tidyselect::starts_with("f.6177."),
    tidyselect::starts_with("f.6153."),
    tidyselect::starts_with("f.20003.")
  )
colnames(clin_data_filtered)
dim(clin_data_filtered)

keep_idx <- which(colnames(clin_data) %in% colnames(clin_data_filtered))
clin_data <- data.table::fread(CLIN_DATA_PATH, select = keep_idx)

# merge clinical data with cmri data
cmri_data <- data.table::fread(CMRI_DATA_PATH)
geno_ids <- readRDS(GENO_IDS_PATH)
clin_data <- dplyr::left_join(
  cmri_data, clin_data,
  by = c("id" = "f.eid")
) %>%
  dplyr::filter(id %in% !!geno_ids)

clin_data_orig <- clin_data

# 6150	Vascular/heart problems diagnosed by doctor	Medical conditions
lvl.100605 <- c(-7, -3, 1, 2, 3, 4)
lbl.100605 <- c(
  "None of the above", "Prefer not to answer", "Heart attack",
  "Angina", "Stroke", "High blood pressure"
)
clin_data_orig %>%
  dplyr::select(tidyselect::starts_with("f.6150.")) %>%
  unlist() %>%
  table()
clin_data <- clin_data %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::starts_with("f.6150."),
      ~ ordered(.x, levels = lvl.100605, labels = lbl.100605)
    )
  )

# 2443	Diabetes diagnosed by doctor	Medical conditions
lvl.100349 <- c(-3, -1, 0, 1)
lbl.100349 <- c("Prefer not to answer", "Do not know", "No", "Yes")
clin_data_orig %>%
  dplyr::select(tidyselect::starts_with("f.2443.")) %>%
  unlist() %>%
  table()
clin_data <- clin_data %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::starts_with("f.2443."),
      ~ ordered(.x, levels = lvl.100349, labels = lbl.100349)
    )
  )

# 2473	Other serious medical condition/disability diagnosed by doctor	Medical conditions
lvl.100603 <- c(-3, -1, 0, 1)
lbl.100603 <- c("Prefer not to answer", "Do not know", "No", "Yes - you will be asked about this later by an interviewer")
clin_data_orig %>%
  dplyr::select(tidyselect::starts_with("f.2473.")) %>%
  unlist() %>%
  table()
clin_data <- clin_data %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::starts_with("f.2443."),
      ~ ordered(.x, levels = lvl.100603, labels = lbl.100603)
    )
  )

# 6177	Medication for cholesterol, blood pressure or diabetes	Medication
lvl.100625 <- c(-7, -3, -1, 1, 2, 3)
lbl.100625 <- c("None of the above", "Prefer not to answer", "Do not know", "Cholesterol lowering medication", "Blood pressure medication", "Insulin")
clin_data_orig %>%
  dplyr::select(tidyselect::starts_with("f.6177.")) %>%
  unlist() %>%
  table()
clin_data <- clin_data %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::starts_with("f.6177."),
      ~ ordered(.x, levels = lvl.100625, labels = lbl.100625)
    )
  )

# 6153	Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones	Medication
lvl.100626 <- c(-7, -3, -1, 1, 2, 3, 4, 5)
lbl.100626 <- c("None of the above", "Prefer not to answer", "Do not know", "Cholesterol lowering medication", "Blood pressure medication", "Insulin", "Hormone replacement therapy", "Oral contraceptive pill or minipill")
clin_data_orig %>%
  dplyr::select(tidyselect::starts_with("f.6153.")) %>%
  unlist() %>%
  table()
clin_data <- clin_data %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::starts_with("f.6153."),
      ~ ordered(.x, levels = lvl.100626, labels = lbl.100626)
    )
  )

get_prevalence <- function(df, x, prefix) {
  df %>%
    dplyr::select(tidyselect::starts_with(sprintf("f.%s.", prefix))) %>%
    dplyr::mutate(
      dplyr::across(tidyselect::everything(), as.character)
    ) %>%
    apply(., MARGIN = 1, FUN = function(.x) any(.x %in% x)) %>%
    mean(na.rm = TRUE)
}

all_codes <- clin_data %>%
  dplyr::select(
    tidyselect::starts_with("f.41202."),
    tidyselect::starts_with("f.41204.")
  ) %>%
  unlist() %>%
  unique()
all_codes <- all_codes[!is.na(all_codes)]

htn_icd10 <- purrr::map(
  paste0("I", 10:16),
  ~ all_codes[startsWith(all_codes, .x)]
) %>%
  do.call(c, .)
as_icd10 <- all_codes[startsWith(all_codes, "I35")]
diabetes_icd10 <- all_codes[startsWith(all_codes, "E11")]
hf_icd10 <- all_codes[startsWith(all_codes, "I50")]

patient_table <- clin_data %>%
  dplyr::summarise(
    `Mean iLVM` = mean(`iLVM (g/m2)`, na.rm = TRUE),
    `SD iLVM` = sd(`iLVM (g/m2)`, na.rm = TRUE),
    `Mean Age` = mean(`age (y)`, na.rm = TRUE),
    `SD Age` = sd(`age (y)`, na.rm = TRUE),
    `Mean Weight` = mean(`weight (kg)`, na.rm = TRUE),
    `SD Weight` = sd(`weight (kg)`, na.rm = TRUE),
    `Mean Height` = mean(`height (cm)`, na.rm = TRUE),
    `SD Height` = sd(`height (cm)`, na.rm = TRUE),
    Male = sum(f.31.0.0 == 1),
    Female = sum(f.31.0.0 == 0)
  ) %>%
  dplyr::mutate(
    `HTN Diagnosed By Doctor` = get_prevalence(
      clin_data, "High blood pressure", 6150
    ),
    `Diabetes Diagnosed By Doctor` = get_prevalence(
      clin_data, "Yes", 2443
    ),
    `Cholesterol Lowering Medication` = get_prevalence(
      clin_data, "Cholesterol lowering medication", 6177
    ),
    `Cholesterol Lowering Medication2` = get_prevalence(
      clin_data, "Cholesterol lowering medication", 6153
    ),
    `Blood Pressure Medication` = get_prevalence(
      clin_data, "Blood pressure medication", 6177
    ),
    `Blood Pressure Medication2` = get_prevalence(
      clin_data, "Blood pressure medication", 6153
    ),
    `Insulin` = get_prevalence(
      clin_data, "Insulin", 6177
    ),
    `Insulin2` = get_prevalence(
      clin_data, "Insulin", 6153
    ),
    `HTN ICD10` = get_prevalence(
      clin_data, htn_icd10, c(41202, 41204)
    ),
    `AS ICD10` = get_prevalence(
      clin_data, as_icd10, c(41202, 41204)
    ),
    `Diabetes ICD10` = get_prevalence(
      clin_data, diabetes_icd10, c(41202, 41204)
    ),
    `HF ICD10` = get_prevalence(
      clin_data, hf_icd10, c(41202, 41204)
    ),
    `Self-diagnosed HTN` = get_prevalence(
      clin_data, "1065", 20002
    ),
    `Self-diagnosed AS` = get_prevalence(
      clin_data, "1490", 20002
    ),
    `Self-diagnosed HF` = get_prevalence(
      clin_data, "1076", 20002
    ),
    `Self-diagnosed Diabetes` = get_prevalence(
      clin_data, "1223", 20002
    ),
    HTN = get_prevalence(
      clin_data, c(htn_icd10, "1065"), c(41202, 41204, 20002)
    ),
    AS = get_prevalence(
      clin_data, c(as_icd10, "1490"), c(41202, 41204, 20002)
    ),
    HF = get_prevalence(
      clin_data, c(hf_icd10, "1076"), c(41202, 41204, 20002)
    ),
    Diabetes = get_prevalence(
      clin_data, c(diabetes_icd10, "1223"), c(41202, 41204, 20002)
    )
  )

format_percentage <- function(x) {
  paste0(formatC(x * 100, digits = 1, format = "f"), "%")
}

clin_data_cleaned <- clin_data %>%
  dplyr::filter(!is.na(`iLVM (g/m2)`)) %>%
  dplyr::rename(
    "Sex" = "f.31.0.0",
    "Age (y)" = "age (y)",
    "Height (cm)" = "height (cm)",
    "Weight (kg)" = "weight (kg)",
  ) %>%
  dplyr::mutate(
    Sex = ifelse(Sex == 0, "Female", "Male")
  )

male_diagnoses <- clin_data_cleaned %>%
  dplyr::filter(Sex == "Male") %>%
  dplyr::summarise(
    Sex = "Male",
    HTN = get_prevalence(
      ., c(htn_icd10, "1065", "High blood pressure"), c(41202, 41204, 20002, 6150)
    ),
    AS = get_prevalence(
      ., c(as_icd10, "1490"), c(41202, 41204, 20002)
    ),
    HF = get_prevalence(
      ., c(hf_icd10, "1076"), c(41202, 41204, 20002)
    ),
    `Type II Diabetes` = get_prevalence(
      ., c(diabetes_icd10, "1223"), c(41202, 41204, 20002)
    ),
    `Blood Pressure Medication` = get_prevalence(
      ., "Blood pressure medication", c(6153, 6177)
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(-Sex, format_percentage)
  )

female_diagnoses <- clin_data_cleaned %>%
  dplyr::filter(Sex == "Female") %>%
  dplyr::summarise(
    Sex = "Female",
    HTN = get_prevalence(
      ., c(htn_icd10, "1065", "High blood pressure"), c(41202, 41204, 20002, 6150)
    ),
    AS = get_prevalence(
      ., c(as_icd10, "1490"), c(41202, 41204, 20002)
    ),
    HF = get_prevalence(
      ., c(hf_icd10, "1076"), c(41202, 41204, 20002)
    ),
    `Type II Diabetes` = get_prevalence(
      ., c(diabetes_icd10, "1223"), c(41202, 41204, 20002)
    ),
    `Blood Pressure Medication` = get_prevalence(
      ., "Blood pressure medication", c(6153, 6177)
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(-Sex, format_percentage)
  )

patient_table <- clin_data_cleaned %>%
  dplyr::group_by(Sex) %>%
  dplyr::summarise(
    N = dplyr::n(),
    dplyr::across(
      c(`iLVM (g/m2)`, `LVM (g)`, `Age (y)`, `Height (cm)`, `Weight (kg)`),
      ~ sprintf(
        "%s (%s)",
        formatC(mean(.x, na.rm = TRUE), digits = 1, format = "f"),
        formatC(sd(.x, na.rm = TRUE), digits = 1, format = "f")
      )
    )
  )

patient_table_all <- dplyr::left_join(
  patient_table,
  dplyr::bind_rows(
    male_diagnoses, female_diagnoses
  ),
  by = "Sex"
) %>%
  dplyr::rename(
    `Hypertensive Diseases` = HTN,
    `Aortic Stenosis` = AS,
    `Heart Failure` = HF,
    `LVMi (g/m2)` = `iLVM (g/m2)`
  )

vthemes::pretty_DT(
  patient_table_all,
  rownames = FALSE, options = list(dom = "t", ordering = F)
)

tab <- patient_table_all %>%
  tibble::column_to_rownames("Sex") %>%
  t() %>%
  as.data.frame()
write.csv(tab, file.path(TABLES_SUPP_DIR, "pop_characteristics.csv"))

tab %>%
  vthemes::pretty_DT(
    options = list(
      dom = "t", ordering = FALSE,
      columnDefs = list(list(className = "dt-center", targets = "_all"))
    )
  )

tab %>%
  tibble::rownames_to_column() %>%
  vthemes::pretty_kable(
    align = "ccc", col.names = c(" ", "Female", "Male"), full_width = FALSE
  )
