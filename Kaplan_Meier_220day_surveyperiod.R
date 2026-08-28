# Kaplan-Meier survival analysis by treatment through 220 days 

library(readxl)
library(dplyr)
library(survival)

# -----------------------------------------------------------------------------
# 1. Select and read the data
# -----------------------------------------------------------------------------

survey_days <- 220

# Opens a window for selecting the Excel workbook.
data_path <- file.choose()

fawns_raw <- read_excel(
  path = data_path,
  .name_repair = "unique"
)

# -----------------------------------------------------------------------------
# 2. Prepare the survival data
# -----------------------------------------------------------------------------

km_data <- fawns_raw %>%
  mutate(
    # Convert treatment codes into treatment labels.
    treatment = dplyr::recode(
      tolower(trimws(as.character(treatment))),
      "1" = "T1",
      "t1" = "T1",
      "2" = "T2",
      "t2" = "T2",
      "3" = "T3",
      "t3" = "T3",
      "4" = "T4",
      "t4" = "T4",
      "c" = "SWSP",
      "swsp" = "SWSP",
    ),
    
    status = tolower(
      trimws(as.character(Status))
    ),
    
    collar = tolower(
      trimws(as.character(collar))
    ),
    
    time_to_death = as.numeric(`time to death`)
  ) %>%
  filter(
    collar == "yes",
    status %in% c("dead", "alive", "censor")
  ) %>%
  mutate(
    # Alive fawns are censored at day 220.
    # Earlier collar losses remain censored on their observed day.
    # Deaths after day 220 are censored at day 220.
    analysis_time = pmin(
      coalesce(time_to_death, survey_days),
      survey_days
    ),
    
    event = as.integer(
      status == "dead" &
        !is.na(time_to_death) &
        time_to_death <= survey_days
    ),
    
    treatment = factor(
      treatment,
      levels = c("T1", "T2", "T3", "T4", "SWSP")
    )
  )

# -----------------------------------------------------------------------------
# 3. Check sample sizes and events
# -----------------------------------------------------------------------------

sample_summary <- km_data %>%
  group_by(treatment) %>%
  summarise(
    n = n(),
    deaths = sum(event),
    censored = sum(event == 0),
    .groups = "drop"
  )

print(sample_summary)

# -----------------------------------------------------------------------------
# 4. Fit the Kaplan-Meier curves
# -----------------------------------------------------------------------------

km_fit <- survfit(
  Surv(
    time = analysis_time,
    event = event
  ) ~ treatment,
  data = km_data
)

# Report survival estimates at day 220.
print(
  summary(
    km_fit,
    times = survey_days,
    extend = TRUE
  )
)

# -----------------------------------------------------------------------------
# 5. Conduct the overall log-rank test
# -----------------------------------------------------------------------------

log_rank <- survdiff(
  Surv(
    time = analysis_time,
    event = event
  ) ~ treatment,
  data = km_data
)

log_rank_df <- length(log_rank$n) - 1

log_rank_p <- pchisq(
  log_rank$chisq,
  df = log_rank_df,
  lower.tail = FALSE
)

print(log_rank)

cat(
  "Overall log-rank p-value:",
  format.pval(log_rank_p, digits = 3),
  "\n"
)

# -----------------------------------------------------------------------------
# 6. Display the Kaplan-Meier curves
# -----------------------------------------------------------------------------

treatment_colors <- c(
  "#D55E00", # T1
  "#0072B2", # T2
  "#009E73", # T3
  "#CC79A7", # T4
  "#333333"  # SWSP
)

plot(
  km_fit,
  col = treatment_colors,
  lwd = 2,
  mark.time = TRUE,
  conf.int = FALSE,
  xlim = c(0, survey_days),
  ylim = c(0, 1),
  xaxt = "n",
  xlab = "Days since release",
  ylab = "Survival probability",
  main = paste0(
    survey_days,
    "-day juvenile white-tailed deer survival"
  )
)

axis(
  side = 1,
  at = seq(
    from = 0,
    to = survey_days,
    by = 20
  )
)

legend(
  "bottomleft",
  legend = c("T1", "T2", "T3", "T4", "SWSP"),
  col = treatment_colors,
  lwd = 2,
  title = "Treatment",
  bty = "n"
)