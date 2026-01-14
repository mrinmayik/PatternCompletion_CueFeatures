# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/PatternCompletion_CueFeatures/")
source("Initialise.R")

############## Read in behavioural data and organise it ############

# Initialize empty dataframe to store all participants' data
behavioural_data <- c()

# Loop through all participants and read their task event files
for(part in unique(participants$Participant)){
  # Read TSV file containing task events for each participant
  part_data <- read.table(paste0(data_path, part, "/", part, "_Data.tsv"),
                          header=TRUE,
                          sep = "\t")
  # Combine with existing data
  behavioural_data <- bind_rows(behavioural_data, part_data)
}

# Filter to only include Test phase trials
test_data <- behavioural_data %>% 
  filter(Phase == "Test")

# Prepare test data for analysis
# Create accuracy variable and set factor levels for proper statistical modeling
test_data <- test_data %>% 
  # Determine what the correct response should be based on trial type
  # "Old" trials should be responded to as "Old", all others as "New"
  mutate(CorrectResponse = ifelse(TrialType == "Old",
                                  "Old", 
                                  "New"),
         # Calculate accuracy: compare participant's response to correct response
         Accuracy = CorrectResponse == Response,
         # Convert variables to factors with specified levels for proper statistical analysis
         Group = factor(Group,
                        levels = factor_labels$Group$levels),
         TrialType = factor(TrialType,
                            levels = factor_labels$TrialType$levels),
         Scrambling = factor(Scrambling,
                             levels = factor_labels$Scrambling$levels))

# Count total number of trials per condition combination
# This is used later to calculate percentages and identify missing data
# Groups by: Participant, Group (YA/OA), TrialType (Old/Similar/New), 
#            Degradation level (40%/60%/85%/95%), 
#            and Scrambling condition (Intact/Scrambled)
total_trials_GrpDegTtScr <- test_data %>% 
  group_by(Participant, Group, TrialType, 
           Degradation, Scrambling) %>% 
  summarise(TotalTrials = length(Scrambling))

############## Analyse Test Data ############

########## Trials without a response ########

# Identify trials that should be excluded from analysis
# Bad trials are those with missing RT (no response) or RT <= 100ms (likely anticipatory responses)
bad_trials <- test_data %>% 
  filter(is.na(RT) | RT <= 100)

# Calculate percentage of bad trials per condition
# This helps assess data quality and whether exclusion rates differ across conditions
badtrials_GrpDegTtScr <- bad_trials %>% 
  group_by(Participant, Group, TrialType, 
           Degradation, Scrambling) %>%
  summarise(BadTrials = length(Scrambling)) %>% 
  ungroup() %>% 
  # Join with total trial counts to calculate percentages
  full_join(total_trials_GrpDegTtScr,
            by = c("Participant", "Group", "TrialType", 
                   "Degradation", "Scrambling")) %>% 
  # Replace NA with 0 for conditions with no bad trials
  mutate(BadTrials = ifelse(is.na(BadTrials),
                            0,
                            BadTrials),
         # Calculate percentage of bad trials
         PercBadTrials = (BadTrials/TotalTrials)*100)

# ANOVA to test if bad trial rates differ across conditions
# Tests whether exclusion rates vary by Group (between-subjects), 
# TrialType, Degradation, and Scrambling (within-subjects factors)
badtrials_GrpDegTtScr_anova <- anova_test(data = badtrials_GrpDegTtScr,
                                          dv = PercBadTrials,
                                          wid = Participant,
                                          between = Group,
                                          within = c(TrialType, Degradation, Scrambling),
                                          type = 3)

# Probe the Group x Degradation interaction in the proportion of bad trials
# Aggregate bad trials by Group and Degradation only (collapsing across TrialType and Scrambling)
badtrials_GrpDeg <- bad_trials %>% 
  group_by(Participant, Group, Degradation) %>%
  summarise(BadTrials = length(Scrambling)) %>% 
  ungroup() %>% 
  # Join with total trials aggregated at same level
  full_join(test_data %>% 
              group_by(Participant, Group, Degradation) %>% 
              summarise(TotalTrials = length(Scrambling)),
            by = c("Participant", "Group", "Degradation")) %>% 
  mutate(BadTrials = ifelse(is.na(BadTrials),
                            0,
                            BadTrials),
         PercBadTrials = (BadTrials/TotalTrials)*100)

# Post-hoc paired t-tests comparing bad trial rates across degradation levels
# Performed separately for each group to test if exclusion rates change with degradation
badtrials_GrpDeg_posthoc <- badtrials_GrpDeg %>% 
  group_by(Group) %>% 
  t_test(formula = PercBadTrials ~ Degradation,
         paired = T)

########## Accuracy ########

# Prepare accuracy data for analysis
acc_data_GrpDegTtScr <- test_data %>% 
  # Exclude bad trials (missing RT or RT <= 100ms) before calculating accuracy
  filter(!(is.na(RT) | RT <= 100)) %>% 
  # Calculate number of correct and incorrect trials per condition
  group_by(Participant, Group, TrialType, 
           Degradation, Scrambling) %>%
  summarise(NumCorr = sum(Accuracy),
            NumIncorr = sum(!Accuracy)) %>% 
  # Join with total trial counts (including bad trials) to calculate percentage accuracy
  # This ensures accuracy is calculated relative to all trials, not just valid ones
  full_join(total_trials_GrpDegTtScr,
            by = c("Participant", "Group", 
                   "TrialType", "Degradation", 
                   "Scrambling")) %>% 
  # Calculate percentage accuracy
  mutate(PercAcc = (NumCorr/TotalTrials)*100)

# Set up contrast coding for Scrambling factor
contrasts(acc_data_GrpDegTtScr$Scrambling) <- rbind(-1, 1)
colnames(contrasts(acc_data_GrpDegTtScr$Scrambling)) <- 
  levels(acc_data_GrpDegTtScr$Scrambling)[2]

# Fit linear mixed-effects model for accuracy
# Full model includes all main effects and interactions:
# - Degradation (continuous): image degradation level
# - TrialType (categorical): Old, Similar, or New
# - Scrambling (categorical): Intact vs Scrambled
# - Group (categorical): YA (Young Adults) vs OA (Older Adults)
# Random intercept for Participant accounts for individual differences
acc_GrpDegTtScr_lmer <- lmerTest::lmer(
  formula = PercAcc ~ Degradation * TrialType * Scrambling * Group +
    (1|Participant), 
  data = acc_data_GrpDegTtScr)

# Perform backward stepwise model selection
# Removes non-significant terms to find the most parsimonious model
acc_GrpDegTtScr_lmerstep <- lmerTest::step(acc_GrpDegTtScr_lmer, 
                                           direction = "backward")
# Extract the final model after stepwise selection
acc_GrpDegTtScr_lmerfinal <- lmerTest::get_model(acc_GrpDegTtScr_lmerstep)

##### Group x Degradation x Trial Type interaction #####

# Calculate estimated marginal means (EMMs) for Group x TrialType x Degradation
# EMMs provide predicted accuracy values at specific degradation levels (40, 70, 85, 95)
# These are averaged across Scrambling conditions (marginalized)
acc_GrpDegTt_emmeans <- emmeans::emmeans(
  acc_GrpDegTtScr_lmerfinal, 
  ~ Group * TrialType * Degradation,
  at = list(Degradation = c(40, 70, 85, 95)))

# Test if accuracy differs from chance (50%) for each condition
# Two-sided test (side = 0) tests both above and below chance
acc_GrpDegTt_chance <- emmeans::test(acc_GrpDegTt_emmeans, 
                                     null = 50, 
                                     side = 0)
# Apply FDR correction for multiple comparisons
acc_GrpDegTt_chance <- correct_p_vals(acc_GrpDegTt_chance, 
                                      "p.value")

# Calculate linear trends (slopes) for Degradation effect
# Tests how accuracy changes as a function of degradation level
# Computed separately for each Group x TrialType combination
acc_GrpDegTt_emtrends <- emmeans::emtrends(
  acc_GrpDegTtScr_lmerfinal, 
  ~ Group * TrialType, var = "Degradation")
# Test if slopes differ from zero (i.e., is there a degradation effect?)
acc_GrpDegTt_trend <- emmeans::test(acc_GrpDegTt_emtrends, 
                                    null = 0, 
                                    side = 0)
# Apply FDR correction for multiple comparisons
acc_GrpDegTt_trend <- correct_p_vals(acc_GrpDegTt_trend, 
                                     "p.value")

# Calculate median degrees of freedom for effect size calculations
# Used as a single value to approximate df across all comparisons
median_df <- median(acc_GrpDegTt_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

# Perform pairwise contrasts to compare groups within each condition
# Combines contrast results with effect sizes for comprehensive reporting
acc_GrpDegTt_contrast <- bind_rows(
  # Pairwise comparisons between groups (YA vs OA) at each TrialType x Degradation combination
  full_join(
    # Contrast
    emmeans::contrast(acc_GrpDegTt_emmeans, 
                      method = "pairwise", 
                      by = c("TrialType", "Degradation")) %>% 
      as_tibble() %>% 
      rename(Within1 = TrialType,
             Within2 = Degradation),
    # Cohen's d effect size
    emmeans::eff_size(acc_GrpDegTt_emmeans,
                      by = c("TrialType", "Degradation"),
                      sigma = sigma(lmerTest::get_model(acc_GrpDegTtScr_lmerstep)),
                      edf = median_df) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = TrialType,
             Within2 = Degradation),
    by = c("contrast", "Within1", "Within2")),
  # Pairwise comparisons of degradation slopes between groups (YA vs OA) within each TrialType
  emmeans::contrast(acc_GrpDegTt_emtrends, 
                    method = "pairwise", 
                    by = "TrialType") %>% 
    as_tibble() %>% 
    rename(Within1 = TrialType) %>% 
    mutate(Within2 = NA)) %>% 
  as.data.frame() %>% 
  # Apply FDR correction to all p-values
  do(correct_p_vals(., "p.value"))

# Create annotation text for plot showing degradation slopes and significance
# Formats slope coefficients and p-values for display in the figure
# Uses plotmath syntax for mathematical notation in ggplot
acc_GrpDegTt_annot <- acc_GrpDegTt_trend %>% 
  mutate(
    # Format slope to 2 decimal places
    slope_str = sprintf("%.2f", Degradation.trend),
    # Format p-values: < .001, = .001, or exact value
    p_str = case_when(
      round(`p.value_fdrcorrected`, 3) < .001 ~ "<.001",
      round(`p.value_fdrcorrected`, 3) == .001 ~ "=.001",
      round(`p.value_fdrcorrected`, 3) > .001 ~ sprintf("%.3f", `p.value_fdrcorrected`)),
    # Remove leading zero from p-values (e.g., "0.05" -> ".05")
    p_str = gsub("^0\\.", ".", p_str),
    # Create annotation string with slope, p-value, and significance marker
    # Different formatting depending on whether p-value uses <, =, or exact value
    annot =
      case_when(
        grepl("^<", p_str) ~
          paste0(
            Group, ": italic(beta)==", slope_str,
            '~","~italic(p)<"', sub("^<", "", p_str), '"',
            ifelse(sig_fdrcorrected == "*", '~"*"', '')),
        grepl("^=", p_str) ~
          paste0(
            Group, ": italic(beta)==", slope_str,
            '~","~italic(p)=="', sub("^=", "", p_str), '"',
            ifelse(sig_fdrcorrected == "*", '~"*"', '')),
        TRUE ~
          paste0(
            Group, ": italic(beta)==", slope_str,
            '*","~italic(p)=="', p_str, '"',
            ifelse(sig_fdrcorrected == "*", '~"*"', '')))) %>%
  # Set annotation position: left side of plot, bottom, with vertical offset for each group
  mutate(x_pos = -Inf,
         y_pos = acc_GrpDegTt_emmeans %>% 
           as_tibble() %>% 
           pull(emmean) %>% 
           min(),
         vjust = ifelse(Group == "YA", 1.3, 2.6),
         hjust = -0.05)

# Create line plot showing accuracy as a function of degradation
# Separate panels for each TrialType, with separate lines for each Group
acc_GrpDegTt_line <- 
  ggplot(acc_GrpDegTt_emmeans %>% 
           as_tibble(), 
         aes(x = Degradation, 
             y = emmean,
             colour = Group)) +
  # Create separate columns for each TrialType (Old, Similar, New)
  facet_grid(cols = vars(TrialType),
             scales = "free") +
  geom_line(size = 1.5,
            aes(group = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 2, color = "black") +
  geom_point(size = 3, aes(group = Group, colour = Group)) +
  # Add annotation text showing slope coefficients and p-values
  geom_text(
    data = acc_GrpDegTt_annot,
    aes(x = x_pos, y = y_pos, label = annot, hjust = hjust, 
        vjust = vjust, colour = Group),
    parse = TRUE,
    size = 4,
    show.legend = FALSE) +
  # Set axis labels and legend titles
  labs(x = factor_labels$Degradation$axis_label, 
       y = "Estimated Accuracy", 
       fill = factor_labels$Group$axis_label, 
       color = factor_labels$Group$axis_label) +
  # Apply color scheme for groups
  scale_fill_manual(values = factor_labels$Group$colours) +
  scale_color_manual(values = factor_labels$Group$colours) +
  # Set y-axis breaks and expansion
  scale_y_continuous(breaks = c(25, 50, 75, 100), 
                     expand = expansion(mult = c(0, 0.1))) +
  # Set x-axis breaks at actual degradation levels
  scale_x_continuous(breaks = c(40, 70, 85, 95),
                     expand = expansion(mult = c(0.1, 0.1))) +
  # Set y-axis limits (allows slight overflow for annotations)
  coord_cartesian(ylim = c(25, 102)) +
  # Add horizontal line at chance level (50%)
  geom_hline(yintercept = 50, linetype = "dashed") + 
  # Apply custom themes for consistent formatting
  y_axis_theme + x_axis_theme +
  blank_bg_theme + legend_theme + paper_facet_theme


##### Group x Trial Type x Scrambling interaction #####

# Calculate estimated marginal means for Group x TrialType x Scrambling
# These are averaged across Degradation levels (marginalized)
# Shows the effect of scrambling condition separately for each group and trial type
acc_GrpTtScr_emmeans <- emmeans::emmeans(
  acc_GrpDegTtScr_lmerfinal, 
  ~ Group * TrialType * Scrambling)
# Test if accuracy differs from chance (50%) for each condition combination
acc_GrpTtScr_chance <- emmeans::test(acc_GrpTtScr_emmeans, 
                                     null = 50, 
                                     side = 0)
# Apply FDR correction for multiple comparisons
acc_GrpTtScr_chance <- correct_p_vals(acc_GrpTtScr_chance, 
                                      "p.value")

# Reuse median degrees of freedom from previous analysis for effect size calculations
median_df <- median(acc_GrpDegTt_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

# Perform pairwise contrasts for Group x TrialType x Scrambling analysis
# Combines multiple types of comparisons with effect sizes
acc_GrpTtScr_contrast <- bind_rows(
  # Compare groups (YA vs OA) within each TrialType x Scrambling combination
  full_join(
    # Contrast
    emmeans::contrast(acc_GrpTtScr_emmeans, 
                      method = "pairwise", 
                      by = c("TrialType", "Scrambling")) %>% 
      as_tibble() %>% 
      rename(Within1 = TrialType,
             Within2 = Scrambling),
    # Cohen's d effect size
    emmeans::eff_size(acc_GrpTtScr_emmeans, 
                      sigma = sigma(lmerTest::get_model(acc_GrpDegTtScr_lmerstep)), 
                      by = c("TrialType", "Scrambling"),
                      edf = median_df) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = TrialType,
             Within2 = Scrambling),
    by = c("contrast", "Within1", "Within2")),
  # Compare scrambling conditions (Intact vs Scrambled) within each Group x TrialType combination
  # Tests if scrambling affects accuracy differently for each group and trial type
  full_join(
    # Contrast
    emmeans::contrast(acc_GrpTtScr_emmeans, method = "pairwise", by = c("TrialType", "Group")) %>%
      as_tibble() %>%
      rename(Within1 = TrialType,
             Within2 = Group),
    # Cohen's d effect size
    emmeans::eff_size(acc_GrpTtScr_emmeans,
                      sigma = sigma(lmerTest::get_model(acc_GrpDegTtScr_lmerstep)), 
                      by = c("TrialType", "Group"),
                      edf = median_df) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = TrialType,
             Within2 = Group),
    by = c("contrast", "Within1", "Within2"))) %>% 
  as.data.frame() %>% 
  # Apply FDR correction to all p-values
  do(correct_p_vals(., "p.value"))

# Create line plot showing accuracy as a function of scrambling condition
# Separate panels for each TrialType, with separate lines for each Group
acc_GDTt.predline <- 
  ggplot(acc_GrpTtScr_emmeans %>% 
           as_tibble(), 
         aes(x = Scrambling, 
             y = emmean,
             colour = Group)) +
  # Create separate columns for each TrialType
  facet_grid(cols = vars(TrialType),
             scales = "free") +
  geom_line(size = 1.5,
            aes(group = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.08, color = "black") +
  geom_point(size = 3, aes(group = Group, colour = Group)) +
  # Set axis labels and legend titles
  labs(x = factor_labels$Scrambling$axis_label, 
       y = "Estimated Accuracy", 
       fill = factor_labels$Group$axis_label, 
       color = factor_labels$Group$axis_label) +
  # Apply color scheme for groups
  scale_fill_manual(values = factor_labels$Group$colours) +
  scale_color_manual(values = factor_labels$Group$colours) +
  # Set y-axis breaks and expansion
  scale_y_continuous(breaks = c(25, 50, 75, 100), 
                     expand = expansion(mult = c(0, 0.1))) +
  # Set y-axis limits
  coord_cartesian(ylim = c(25, 102)) +
  # Add horizontal line at chance level (50%)
  geom_hline(yintercept = 50, linetype = "dashed") + 
  # Apply custom themes for consistent formatting
  y_axis_theme + x_axis_theme +
  blank_bg_theme + legend_theme + paper_facet_theme

# Combine both accuracy plots into a single figure
# Top panel: Degradation effect, Bottom panel: Scrambling effect
# Share common legend at bottom
ggarrange(acc_GrpDegTt_line,
          acc_GDTt.predline +
            theme(axis.title.x = element_blank()),
          nrow = 2, ncol = 1, 
          common.legend = TRUE, 
          legend = "bottom")

########## Reaction Time ##########

# Prepare reaction time data for analysis
# Calculate mean RT for each condition combination per participant
rt_data_GrpDegTtScr <- test_data %>%
  # Exclude bad trials (missing RT or RT <= 100ms) from RT analysis
  filter(!(is.na(RT) | RT <= 100)) %>% 
  # Calculate mean RT for each condition combination
  group_by(Participant, Group, TrialType, Degradation, 
           Scrambling) %>%
  summarise(MeanRT = mean(RT))

# Set up contrast coding for Scrambling factor (same as for accuracy)
contrasts(rt_data_GrpDegTtScr$Scrambling) <- rbind(-1, 1)
colnames(contrasts(rt_data_GrpDegTtScr$Scrambling)) <- 
  levels(rt_data_GrpDegTtScr$Scrambling)[2]

# Fit linear mixed-effects model for reaction time
# Full model includes all main effects and interactions:
# - Degradation (continuous): image degradation level
# - TrialType (categorical): Old, Similar, or New
# - Scrambling (categorical): Intact vs Scrambled
# - Group (categorical): YA (Young Adults) vs OA (Older Adults)
# Random intercept for Participant accounts for individual differences
rt_GrpDegTtScr_lmer <- lmerTest::lmer(
  formula = MeanRT ~ Degradation * TrialType * Scrambling * Group +
    (1|Participant), 
  data = rt_data_GrpDegTtScr)

# Perform backward stepwise model selection
# Removes non-significant terms to find the most parsimonious model
rt_GrpDegTtScr_lmerstep <- lmerTest::step(rt_GrpDegTtScr_lmer, 
                                          direction = "backward")
# Extract the final model after stepwise selection
rt_GrpDegTtScr_lmerfinal <- lmerTest::get_model(rt_GrpDegTtScr_lmerstep)

##### Group x Degradation x Trial Type x Scrambling #####

# Calculate estimated marginal means for the four-way interaction
# EMMs provide predicted RT values at specific degradation levels (40, 70, 85, 95)
# Shows RT for each combination of Group, Scrambling, Degradation, and TrialType
rt_GrpDegTtScr_emmeans <- emmeans::emmeans(rt_GrpDegTtScr_lmerfinal, 
                                           ~ Group * Scrambling * Degradation * TrialType,
                                           at = list(Degradation = c(40, 70, 85, 95)))

# Calculate linear trends (slopes) for Degradation effect on RT
# Tests how RT changes as a function of degradation level
# Computed separately for each Group x Scrambling x TrialType combination
rt_GrpDegTtScr_emtrends <- emmeans::emtrends(rt_GrpDegTtScr_lmerfinal, 
                                             ~ Group * Scrambling * TrialType,
                                             var = "Degradation")
# Test if slopes differ from zero (i.e., is there a degradation effect on RT?)
rt_GrpDegTtScr_trend <- emmeans::test(rt_GrpDegTtScr_emtrends, 
                                      null = 0, 
                                      side = 0)
# Apply FDR correction for multiple comparisons
rt_GrpDegTtScr_trend <- correct_p_vals(rt_GrpDegTtScr_trend, 
                                       "p.value")

# Calculate median degrees of freedom for effect size calculations
median_df <- median(rt_GrpDegTtScr_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

# Perform pairwise contrasts for RT analysis
# Combines multiple types of comparisons with effect sizes
rt_GrpDegTtScr_contrast <- bind_rows(
  # Compare groups (YA vs OA) within each Degradation x TrialType x Scrambling combination
  full_join(
    # Contrast
    emmeans::contrast(rt_GrpDegTtScr_emmeans, 
                      method = "pairwise", 
                      by = c("Degradation", "TrialType", "Scrambling")) %>% 
      as_tibble() %>% 
      rename(Within1 = Degradation,
             Within2 = TrialType,
             Within3 = Scrambling),
    # Cohen's d effect size
    emmeans::eff_size(rt_GrpDegTtScr_emmeans,
                      by = c("Degradation", "TrialType", "Scrambling"),
                      edf = median_df,
                      sigma = sigma(lmerTest::get_model(rt_GrpDegTtScr_lmerstep))) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = Degradation,
             Within2 = TrialType,
             Within3 = Scrambling),
    by = c("contrast", "Within1", "Within2", "Within3")),
  # Compare scrambling conditions (Intact vs Scrambled) within each Degradation x TrialType x Group combination
  # Tests if scrambling affects RT differently across conditions
  full_join(
    # Contrast
    emmeans::contrast(rt_GrpDegTtScr_emmeans, method = "pairwise", 
                      by = c("Degradation", "TrialType", "Group")) %>% 
      as_tibble() %>% 
      rename(Within1 = Degradation,
             Within2 = TrialType,
             Within3 = Group),
    # Cohen's d effect size
    emmeans::eff_size(rt_GrpDegTtScr_emmeans,
                      by = c("Degradation", "TrialType", "Group"),
                      edf = median_df,
                      sigma = sigma(lmerTest::get_model(rt_GrpDegTtScr_lmerstep))) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = Degradation,
             Within2 = TrialType,
             Within3 = Group),
    by = c("contrast", "Within1", "Within2", "Within3")),
  # Compare degradation slopes between groups within each TrialType x Scrambling combination
  # Tests if groups differ in how RT changes with degradation
  # Contrast
  emmeans::contrast(rt_GrpDegTtScr_emtrends, 
                    method = "pairwise", 
                    by = c("TrialType", "Scrambling")) %>% 
    as_tibble() %>% 
    rename(Within2 = TrialType,
           Within3 = Scrambling) %>% 
    mutate(Within1 = 1)) %>% 
  as.data.frame() %>% 
  # Apply FDR correction to all p-values
  do(correct_p_vals(., "p.value"))

# Calculate y-axis range for RT plot
# Determines appropriate limits based on the data range
y_range <- range(rt_GrpDegTtScr_emmeans %>% 
                   as_tibble() %>% 
                   pull(emmean), 
                 na.rm = TRUE)
# Add padding above and below the data range
y_range <- c(y_range[1] - 100, y_range[2] + 100)
# Generate evenly spaced breaks for y-axis
y_breaks <- pretty(y_range, n = 3)

# Create annotation text for RT plot showing degradation slopes and significance
# Formats slope coefficients and p-values for display in the figure
rt_GrpDegTtScr_annot <- rt_GrpDegTtScr_trend %>% 
  mutate(
    # Format slope to 2 decimal places
    slope_str = sprintf("%.2f", Degradation.trend),
    # Format p-values: < .001, = .001, or exact value
    p_str = case_when(
      round(`p.value_fdrcorrected`, 3) < .001 ~ "<.001",
      round(`p.value_fdrcorrected`, 3) == .001 ~ "=.001",
      round(`p.value_fdrcorrected`, 3) > .001 ~ sprintf("%.3f", `p.value_fdrcorrected`)
    ),
    # Remove leading zero from p-values (e.g., "0.05" -> ".05")
    p_str = gsub("^0\\.", ".", p_str),
    # Create annotation string with slope, p-value, and significance marker
    annot =
      case_when(
        grepl("^<", p_str) ~
          paste0(
            Group, ": italic(beta)==", slope_str,
            '~","~italic(p)<"', sub("^<", "", p_str), '"',
            ifelse(sig_fdrcorrected == "*", '~"*"', '')
          ),
        grepl("^=", p_str) ~
          paste0(
            Group, ": italic(beta)==", slope_str,
            '~","~italic(p)=="', sub("^=", "", p_str), '"',
            ifelse(sig_fdrcorrected == "*", '~"*"', '')
          ),
        TRUE ~
          paste0(
            Group, ": italic(beta)==", slope_str,
            '*","~italic(p)=="', p_str, '"',
            ifelse(sig_fdrcorrected == "*", '~"*"', '')))) %>%
  # Set annotation position: left side of plot, top, with vertical offset for each group
  mutate(x_pos = -Inf,
         y_pos = y_range[2],  # both bottom
         vjust = ifelse(Group == "YA", 0.5, 1.8),    # stack vertically
         hjust = -0.05)

# Get unique trialtypes (for potential future use in separate plots)
trialtypes <- unique(rt_GrpDegTtScr_emmeans %>% 
                       as_tibble() %>% 
                       pull(TrialType))

# Create line plot showing RT as a function of degradation
# Separate panels for each TrialType (rows) and Scrambling condition (columns)
# Separate lines for each Group
rt_GrpDegTtScr_line <- 
  ggplot(rt_GrpDegTtScr_emmeans %>% 
           as_tibble(), 
         aes(x = Degradation, 
             y = emmean,
             colour = Group)) +
  # Create grid: rows = TrialType, columns = Scrambling
  facet_grid(rows = vars(TrialType),
             cols = vars(Scrambling)) +
  geom_line(size = 1.5, aes(group = Group,
                            colour = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 2, color = "black") +
  geom_point(size = 2, aes(group = Group, colour = Group)) +
  # Add annotation text showing slope coefficients and p-values
  geom_text(
    data = rt_GrpDegTtScr_annot,
    aes(x = x_pos, y = y_pos, label = annot, hjust = hjust, 
        vjust = vjust, colour = Group),
    parse = TRUE,
    size = 3,
    show.legend = FALSE) +
  # Set axis labels and legend titles
  labs(x = factor_labels$Percent$axis_label, 
       y = "Estimated RT", 
       fill = factor_labels$Group$axis_label, 
       color = factor_labels$Group$axis_label) +
  # Apply color scheme for groups
  scale_fill_manual(values = factor_labels$Group$colours) +
  scale_color_manual(values = factor_labels$Group$colours) +
  # Set y-axis limits and breaks (calculated above based on data range)
  scale_y_continuous(limits = y_range,
                     breaks = y_breaks,
                     expand = expansion(mult = c(0, 0.1))) +
  # Set x-axis breaks at actual degradation levels
  scale_x_continuous(breaks = c(40, 70, 85, 95),
                     expand = expansion(mult = c(0.1, 0.1))) +
  # Apply custom themes for consistent formatting
  y_axis_theme + x_axis_theme +
  blank_bg_theme + legend_theme + paper_facet_theme
# Display plot
rt_GrpDegTtScr_line
