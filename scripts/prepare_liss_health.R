# Prepare the LISS mood and physical-health panel for the empirical RI-CLPM.
#
# Input:  the LISS Core Study Health module, one SPSS file per wave, in
#         ../dataset/LISS/ (files chYYw_EN_*.sav from
#         https://www.dataarchive.lissdata.nl), next to this project.
# Output: ../dataset/liss_health_panel.csv (read by scripts/empirical_clpm.R),
#         one row per person and wave:
#           id, wave, year, followed by the 11 variables below,
#         with only persons who answered all 11 items in every selected wave.
#
# Run from the project root:  Rscript scripts/prepare_liss_health.R
#
# Variables (item numbers within wave chYYw), all coded so that higher = worse:
#   Mental Health Inventory-5, how often in the past month the respondent felt
#   (1 = never ... 6 = continuously):
#     anxious        011  very anxious
#     down           012  so down that nothing could cheer them up
#     calm_r         013  calm and peaceful (reversed)
#     depressed      014  depressed and gloomy
#     happy_r        015  happy (reversed)
#   Physical health:
#     poor_health    004  general health, 1 = poor ... 5 = excellent (reversed)
#     hinder_daily   020  physical health or emotional problems hindered daily
#                         activities in the past month (1 = not at all ... 5 = very much)
#     hinder_social  021  idem, social activities
#   Recurrent complaints (0 = no, 1 = yes):
#     pain           070  back, knee, hip or other joint pain
#     fatigue        076  fatigue
#     sleep          077  sleeping problems
#
# Waves. Zainal and Newman (2021, Journal of Affective Disorders) analysed the
# MHI-5 in the LISS Health module from 2008 to 2018 with an RI-CLPM; there is
# no Health wave in 2014, so their window has 10 waves with one two-year gap.
# Our RI-CLPM uses one lag matrix for all transitions, so we take the
# equally spaced first part of that window: the six annual waves 2008-2013.
# The Health module is fielded in November-December of each year.

library(foreign)

waves <- data.frame(
  code = c("08b", "09c", "10d", "11e", "12f", "13g"),
  year = 2008:2013)

items <- data.frame(
  number  = c("011", "012", "013", "014", "015", "004", "020", "021", "070", "076", "077"),
  name    = c("anxious", "down", "calm_r", "depressed", "happy_r", "poor_health",
              "hinder_daily", "hinder_social", "pain", "fatigue", "sleep"),
  lowest  = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
  highest = c(6, 6, 6, 6, 6, 5, 5, 5, 1, 1, 1),
  reverse = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))

# ---- Read the selected waves ----------------------------------------------------

read_wave <- function(code, year, wave_number) {
  file <- list.files("../dataset/LISS", paste0("^ch", code, "_EN_.*\\.sav$"), full.names = TRUE)
  stopifnot(length(file) == 1)
  raw <- suppressWarnings(read.spss(file, to.data.frame = TRUE,
                                    use.value.labels = FALSE))
  columns <- paste0("ch", code, items$number)
  stopifnot(columns %in% names(raw))
  responses <- as.matrix(raw[, columns])
  for (j in seq_len(nrow(items))) {
    stopifnot(all(responses[, j] %in% c(items$lowest[j]:items$highest[j], NA)))
    if (items$reverse[j])
      responses[, j] <- items$lowest[j] + items$highest[j] - responses[, j]
  }
  colnames(responses) <- items$name
  data.frame(id = raw$nomem_encr, wave = wave_number, year = year, responses)
}

panel <- do.call(rbind, Map(read_wave, waves$code, waves$year, seq_len(nrow(waves))))
rownames(panel) <- NULL
stopifnot(!anyDuplicated(panel[c("id", "wave")]))

# ---- Balanced panel: complete responses in every wave -----------------------------

complete <- panel[stats::complete.cases(panel[items$name]), ]
waves_per_person <- table(complete$id)
keep <- as.numeric(names(waves_per_person)[waves_per_person == nrow(waves)])
panel <- complete[complete$id %in% keep, ]
panel <- panel[order(panel$id, panel$wave), ]

# ---- Write and describe ------------------------------------------------------------

write.csv(panel, "../dataset/liss_health_panel.csv", row.names = FALSE)

cat(sprintf("Waves %s (T = %d); %d persons with complete responses in every wave.\n",
            paste(waves$year, collapse = ", "), nrow(waves), length(keep)))
cat("\nMeans by year (higher = worse; complaints are proportions):\n")
print(round(aggregate(panel[items$name], panel["year"], mean), 2), row.names = FALSE)
cat("\nCorrelation of each variable with itself one year later (first two waves):\n")
first <- panel[panel$wave == 1, ]; second <- panel[panel$wave == 2, ]
stopifnot(identical(first$id, second$id))
print(round(diag(cor(first[items$name], second[items$name])), 2))
