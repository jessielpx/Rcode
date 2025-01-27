install.packages("survival")
install.packages("survminer")
install.packages("readxl")
library(survival)
library(survminer)
library(readxl)

# The name of the column by which the analysis is done
index <- "VHL"

# Read in the data and exclude non-localized sample
tumors <- read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin')
tumors <- tumors[tumors$LocalisedRCC == "1", ]

# Filtering
#tumors <- tumors[tumors$VHL == "1", ]
#tumors <- tumors[tumors$Sex == "Male", ]

# Define the 5-year CSS
max_time <- 5
tumors$CSSYears_5yr <- pmin(tumors$CSSYears, max_time)  # Cap the CSS time at 5 years
tumors$CSSEvent_5yr <- ifelse(tumors$CSSYears > max_time, 0, tumors$CSSEvent)  # Mark events beyond 5 years as censored

# Define the survival object & Fit the Kaplan-Meier survival model
tumorKM_5yr <- survfit(Surv(tumors$CSSYears_5yr, tumors$CSSEvent_5yr) ~ tumors[[index]], data = tumors, type = "kaplan-meier")
# Same but stratified by Stage
tumorKM_5yr <- survfit(Surv(tumors$CSSYears_5yr, tumors$CSSEvent_5yr) ~ tumors[[index]] + strata(Stage) , data = tumors, type = "kaplan-meier")



# Perform a log-rank test
logrank_test_5yr <- survdiff(Surv(tumors$CSSYears_5yr, tumors$CSSEvent_5yr) ~ tumors[[index]], data = tumors)
## KG - 26.01.2025 - Additionally, you could check if a cox-PH model is significant
coxph_5yr <- coxph(Surv(tumors$CSSYears_5yr, tumors$CSSEvent_5yr) ~ tumors[[index]], data = tumors)
# Same but stratified by Stage
coxph_5yr <- coxph(Surv(tumors$CSSYears_5yr, tumors$CSSEvent_5yr) ~ tumors[[index]] + strata(Stage), data = tumors)
summary(coxph_5yr)
## KG - 26.01.2025 - for the coxph, if necessary, you may need to refactor to sort the categories (this is an example, may not be necessary and would need to be changed to fit your categories)
tumors$SETD2_TP53 <- factor(tumors[[index]], levels = c("None", "SETD2", "TP53", "SETD2_TP53"))




# Calculate the p-value for the log-rank test
logrank_pval_5yr <- 1 - pchisq(logrank_test_5yr$chisq, df = length(logrank_test_5yr$n) - 1)

# Print the log-rank p-value for reference
print(paste("5-year log-rank p-value:", round(logrank_pval_5yr, 3)))

# Plot the Kaplan-Meier curve with the p-value and risk table
ggsurvplot(
  tumorKM_5yr,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval_5yr, 3)),
  pval.coord = c(3, 1),  # Add the p-value to the plot
  risk.table = TRUE,     # Enable the risk table at the bottom
  legend.labs = c("VHLwt-I", "VHLwt-II", "VHLwt-III", "VHLmut-I","VHLmut-II","VHLmut-III"),
  legend = c(0.1, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("#99FFCC", "#33CC66", "#006633", "#FFCC99","#FF9933","#993300"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)



ggsurvplot(
  tumorKM_5yr,
  conf.int = FALSE,
  pval = paste("p =", signif(logrank_pval_5yr, 3)),
  pval.coord = c(3, 0.95),# Add the p-value to the plot
  risk.table = TRUE,                            # Enable the risk table at the bottom
  legend.labs = c("VHLwt", "VHLmut"),
  legend = c(0.2, 0.2),
  break.time.by = 1,
  legend.title = "Subtype",
  censor.shape = "+",
  censor.size = 5,
  palette = c("#33CC66", "#FF9933"),
  xlab = "Survival (years)",
  ylab = "Proportion surviving"
)

