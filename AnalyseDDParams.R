# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/PatternCompletion_CueFeatures/")
source("Initialise.R")

############## Read in drift diffusion data and organise it ############

# Initialize empty dataframe to store all participants' data
driftdiffusion_data <- c()

# Loop through all participants and read their drift diffusion parameter files
for(part in unique(participants$Participant)){
  # Read TSV file containing drift diffusion parameters for each participant
  part_data <- read.table(paste0(data_path, part, "/", part, "_DriftDiffusionParameters.tsv"),
                          header=TRUE,
                          sep = "\t")
  # Combine with existing data
  driftdiffusion_data <- bind_rows(driftdiffusion_data, part_data)
}

# Prepare test data for analysis
# Set factor levels for proper statistical modeling
driftdiffusion_data <- driftdiffusion_data %>% 
  # Convert variables to factors with specified levels for proper statistical analysis
  mutate(Group = factor(Group,
                       levels = factor_labels$Group$levels),
         TrialType = factor(TrialType,
                           levels = factor_labels$TrialType$levels),
         Scrambling = factor(Scrambling,
                            levels = factor_labels$Scrambling$levels))

############## Analyse Drift Rate (delta) ############

# Set up contrast coding for Scrambling factor
contrasts(driftdiffusion_data$Scrambling) <- rbind(-1, 1)
colnames(contrasts(driftdiffusion_data$Scrambling)) <- 
  levels(driftdiffusion_data$Scrambling)[2]

# Fit linear mixed-effects model for delta
# Full model includes all main effects and interactions:
# - Degradation (continuous): image degradation level
# - TrialType (categorical): Old, Similar, or New
# - Scrambling (categorical): Intact vs Scrambled
# - Group (categorical): YA (Young Adults) vs OA (Older Adults)
# Random intercept for Participant accounts for individual differences
delta_GrpDegTtScr_lmer <- lmerTest::lmer(
  formula = delta ~ Degradation * TrialType * Scrambling * Group +
    (1|Participant), 
  data = driftdiffusion_data)

# Perform backward stepwise model selection
# Removes non-significant terms to find the most parsimonious model
delta_GrpDegTtScr_lmerstep <- lmerTest::step(delta_GrpDegTtScr_lmer, 
                                             direction = "backward")
# Extract the final model after stepwise selection
delta_GrpDegTtScr_lmerfinal <- lmerTest::get_model(delta_GrpDegTtScr_lmerstep)

##### Group x Degradation x Trial Type x Scrambling #####

# Calculate estimated marginal means for the four-way interaction
# EMMs provide predicted delta values at specific degradation levels (40, 70, 85, 95)
# Shows delta for each combination of Group, Scrambling, Degradation, and TrialType
delta_GrpDegTtScr_emmeans <- emmeans::emmeans(delta_GrpDegTtScr_lmerfinal, 
                                              ~ Group * Scrambling * Degradation * TrialType,
                                              at = list(Degradation = c(40, 70, 85, 95)))

# Calculate linear trends (slopes) for Degradation effect on delta
# Tests how delta changes as a function of degradation level
# Computed separately for each Group x Scrambling x TrialType combination
delta_GrpDegTtScr_emtrends <- emmeans::emtrends(delta_GrpDegTtScr_lmerfinal, 
                                                ~ Group * Scrambling * TrialType,
                                                var = "Degradation")
# Test if slopes differ from zero (i.e., is there a degradation effect on delta?)
delta_GrpDegTtScr_trend <- emmeans::test(delta_GrpDegTtScr_emtrends, 
                                         null = 0, 
                                         side = 0)
# Apply FDR correction for multiple comparisons
delta_GrpDegTtScr_trend <- correct_p_vals(delta_GrpDegTtScr_trend, 
                                          "p.value")

# Calculate median degrees of freedom for effect size calculations
median_df <- median(delta_GrpDegTtScr_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

# Perform pairwise contrasts for delta analysis
# Combines multiple types of comparisons with effect sizes
delta_GrpDegTtScr_contrast <- bind_rows(
  # Compare groups (YA vs OA) within each Degradation x TrialType x Scrambling combination
  full_join(
    # Contrast
    emmeans::contrast(delta_GrpDegTtScr_emmeans, 
                      method = "pairwise", 
                      by = c("Degradation", "TrialType", "Scrambling")) %>% 
      as_tibble(),
    # Cohen's d effect size
    emmeans::eff_size(delta_GrpDegTtScr_emmeans,
                      by = c("Degradation", "TrialType", "Scrambling"),
                      edf = median_df,
                      sigma = sigma(lmerTest::get_model(delta_GrpDegTtScr_lmerstep))) %>% 
      as_tibble() %>% 
      select(-c(df, SE)),
    by = c("contrast", "Degradation", "TrialType", "Scrambling")) %>% 
    rename(Within1 = Degradation,
           Within2 = TrialType,
           Within3 = Scrambling),
  # Compare scrambling conditions (Intact vs Scrambled) within each Degradation x TrialType x Group combination
  # Tests if scrambling affects delta differently across conditions
  full_join(
    # Contrast
    emmeans::contrast(delta_GrpDegTtScr_emmeans, method = "pairwise", 
                      by = c("Degradation", "TrialType", "Group")) %>% 
      as_tibble(),
    # Cohen's d effect size
    emmeans::eff_size(delta_GrpDegTtScr_emmeans,
                      by = c("Degradation", "TrialType", "Group"),
                      edf = median_df,
                      sigma = sigma(lmerTest::get_model(delta_GrpDegTtScr_lmerstep))) %>% 
      as_tibble() %>% 
      select(-c(df, SE)),
    by = c("contrast", "Degradation", "TrialType", "Group")) %>% 
    rename(Within1 = Degradation,
           Within2 = TrialType,
           Within3 = Group),
  # Compare degradation slopes between groups within each TrialType x Scrambling combination
  # Tests if groups differ in how delta changes with degradation
  # Contrast
  emmeans::contrast(delta_GrpDegTtScr_emtrends, 
                    method = "pairwise", 
                    by = c("TrialType", "Scrambling")) %>% 
    as_tibble() %>% 
    rename(Within2 = TrialType,
           Within3 = Scrambling) %>% 
    mutate(Within1 = 1)) %>% 
  as.data.frame() %>% 
  # Apply FDR correction to all p-values
  do(correct_p_vals(., "p.value"))

##### Plot #####

# Create annotation text for delta plot showing degradation slopes and significance
# Formats slope coefficients and p-values for display in the figure
delta_GrpDegTtScr_annot <- delta_GrpDegTtScr_trend %>% 
  mutate(
    # Format slope to 3 decimal places
    slope_str = case_when(
      round(abs(Degradation.trend), 3) < .001 ~ "<.001",
      round(abs(Degradation.trend), 3) >= .001 ~ sprintf("%.3f", Degradation.trend)
    ),
    # Create annotation string for slope
    slope_str_m = case_when(
      grepl("^<", slope_str) ~ paste0('<"', sub("^<", "", slope_str), '"'),
      grepl("^=", slope_str) ~ paste0('=="', sub("^=", "", slope_str), '"'),
      TRUE                   ~ paste0('=="', slope_str, '"')
    ),
    # Format p-values: < .001, = .001, or exact value
    p_str = case_when(
      round(p.value_fdrcorrected, 3) < .001 ~ "<.001",
      round(p.value_fdrcorrected, 3) == .001 ~ "=.001",
      TRUE                                  ~ sprintf("%.3f", p.value_fdrcorrected)
    ),
    # Remove leading zero from p-values (e.g., "0.05" -> ".05")
    p_str = gsub("^0\\.", ".", p_str),
    # Create annotation string for p-value
    p_str_m = case_when(
      grepl("^<", p_str) ~ paste0('<"', sub("^<", "", p_str), '"'),
      grepl("^=", p_str) ~ paste0('=="', sub("^=", "", p_str), '"'),
      TRUE               ~ paste0('=="', p_str, '"')
    ),
    # Concatenate and add significance marker
    annot = paste0(
      Group,
      ": italic(beta)", slope_str_m,
      '~","~italic(p)', p_str_m,
      ifelse(sig_fdrcorrected == "*", '~"*"', '')))

# Calculate consistent extent across all trialtypes
# Because the y-ranges of drift rate differ vastly between trial types,
# the patterns are obscured if the same range is maintained.
# We can use the actual values corresponding to each trial type but
# use a common extent (ymax-ymin) to ensure that patterns are visible
# while maintining the values

# Find the range that ensures zero is visible for all trialtypes
all_ranges <- delta_GrpDegTtScr_emmeans %>% 
  as_tibble() %>% 
  group_by(TrialType) %>% 
  summarise(
    min_val = min(lower.CL, na.rm = TRUE),
    max_val = max(upper.CL, na.rm = TRUE),
    # Ensure zero is visible by extending range if needed
    range_min = min(min_val, -0.1),
    range_max = max(max_val, 0.2),
    range_width = range_max - range_min,
    .groups = 'drop'
  )

# Use the maximum range width across all trialtypes
extent <- max(all_ranges$range_width)

# Create line plot showing delta as a function of degradation
# Separate panels for each TrialType (rows) and Scrambling condition (columns)
# Since it is not possible to give y-axis limits separately for each facet,
# we will construct three separate plots that will be stitched together
delta_GrpDegTtScr_line <- lapply(factor_labels$TrialType$levels, function(trialtype_name) {
  # Filter data for current trialtype
  plot_data <- delta_GrpDegTtScr_emmeans %>% 
    as_tibble() %>% 
    filter(TrialType == trialtype_name)
  
  # Calculate range for this trialtype with consistent extent
  # Ensure zero is visible in the range
  trialtype_range <- c(min(plot_data$lower.CL, na.rm = TRUE),
                       max(plot_data$upper.CL))
  trialtype_min <- min(trialtype_range[1], -0.1)
  trialtype_max <- max(trialtype_range[2], 0.1)
  
  # Center the range around the trialtype's data with consistent extent
  trialtype_center <- (trialtype_min + trialtype_max) / 2
  y_range <- c(trialtype_center - extent/2, trialtype_center + extent/2)
  y_breaks <- pretty(y_range, n = 3)
  
  # Set annotation position: left side of plot, top, with vertical offset for each group
  sigtrend_pos <- delta_GrpDegTtScr_annot %>% 
    filter(TrialType == trialtype_name) %>% 
    mutate(x_pos = 37.5,
           y_pos = case_when(trialtype_name == "New" ~ 0.3,
                             trialtype_name == "Old" ~ Inf, 
                             trialtype_name == "Similar" ~ y_range[2]),
           vjust = case_when(trialtype_name == "New" ~ ifelse(Group == "YA", 0.5, 1.8),
                             trialtype_name == "Old" ~ ifelse(Group == "YA", 1.7, 3),
                             trialtype_name == "Similar" ~ ifelse(Group == "YA", -0.7, 0.6)),
           hjust = -0.05)   
  
  # Create plot
  delta.predline <- ggplot(plot_data, 
                           aes(x = Degradation, 
                               y = emmean,
                               colour = Group)) +
    facet_grid(cols = vars(Scrambling),
               rows = vars(TrialType)) +
    geom_line(size = 1.5, aes(group = Group,
                              colour = Group)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 2, color = "black") +
    geom_point(size = 2, aes(group = Group, colour = Group)) +
    geom_text(
      data = sigtrend_pos,
      aes(x = x_pos, y = y_pos, label = annot, hjust = hjust, 
          vjust = vjust, colour = Group),
      parse = TRUE,
      size = 3,
      show.legend = FALSE
    ) +
    # Add customizations
    labs(x = factor_labels$Percent$axis_label, 
         y = "Estimated Drift Rate", 
         fill = factor_labels$Group$axis_label, 
         color = factor_labels$Group$axis_label) +
    scale_fill_manual(values = factor_labels$Group$colours) +
    scale_color_manual(values = factor_labels$Group$colours) +
    # scale_x_discrete(expand = expansion(add = c(0.7, 0.2))) +  # More space on left, some on right
    scale_y_continuous(limits = y_range,
                       breaks = y_breaks,
                       expand = expansion(mult = c(0.05, 0.25)),
                       labels = function(x) sprintf("%4.1f", x)) +
    scale_x_continuous(breaks = c(40, 70, 85, 95)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    y_axis_theme + x_axis_theme +
    blank_bg_theme + legend_theme + paper_facet_theme
  
  return(delta.predline)
})

ggarrange(delta_GrpDegTtScr_line[1][[1]] +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_blank()),
          delta_GrpDegTtScr_line[2][[1]] +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.text.x.top = element_blank()),
          delta_GrpDegTtScr_line[3][[1]] + 
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.text.x.top = element_blank()),
          nrow = 3, ncol = 1, 
          common.legend = TRUE, 
          legend = "bottom")

############## Analyse Non-decision Time (tau) ############

# Fit linear mixed-effects model for tau
# Full model includes all main effects and interactions:
# - Degradation (continuous): image degradation level
# - TrialType (categorical): Old, Similar, or New
# - Scrambling (categorical): Intact vs Scrambled
# - Group (categorical): YA (Young Adults) vs OA (Older Adults)
# Random intercept for Participant accounts for individual differences
tau_GrpDegTtScr_lmer <- lmerTest::lmer(
  formula = tau ~ Degradation * TrialType * Scrambling * Group +
    (1|Participant), 
  data = driftdiffusion_data)

# Perform backward stepwise model selection
# Removes non-significant terms to find the most parsimonious model
tau_GrpDegTtScr_lmerstep <- lmerTest::step(tau_GrpDegTtScr_lmer, 
                                           direction = "backward")
# Extract the final model after stepwise selection
tau_GrpDegTtScr_lmerfinal <- lmerTest::get_model(tau_GrpDegTtScr_lmerstep)

##### Group x Degradation x Scrambling #####

# Calculate estimated marginal means for the three-way interaction
# EMMs provide predicted tau values at specific degradation levels (40, 70, 85, 95)
# Shows tau for each combination of Group, Scrambling, and Degradation
# Marginalized across TrialType
tau_GrpDegScr_emmeans <- emmeans::emmeans(tau_GrpDegTtScr_lmerfinal, 
                                          ~ Group * Scrambling * Degradation,
                                          at = list(Degradation = c(40, 70, 85, 95)))

# Calculate linear trends (slopes) for Degradation effect on tau
# Tests how tau changes as a function of degradation level
# Computed separately for each Group x Scrambling combination
# Marginalized across TrialType
tau_GrpDegScr_emtrends <- emmeans::emtrends(tau_GrpDegTtScr_lmerfinal, 
                                            ~ Group * Scrambling,
                                            var = "Degradation")
# Test if slopes differ from zero (i.e., is there a degradation effect on tau?)
tau_GrpDegScr_trend <- emmeans::test(tau_GrpDegScr_emtrends, 
                                     null = 0, 
                                     side = 0)
# Apply FDR correction for multiple comparisons
tau_GrpDegScr_trend <- correct_p_vals(tau_GrpDegScr_trend, 
                                      "p.value")

# Calculate median degrees of freedom for effect size calculations
median_df <- median(tau_GrpDegScr_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

# Perform pairwise contrasts for tau analysis
# Combines multiple types of comparisons with effect sizes
tau_GrpDegScr_contrast <- bind_rows(
  # Compare groups (YA vs OA) within each Degradation x Scrambling combination
  full_join(
    # Contrast
    emmeans::contrast(tau_GrpDegScr_emmeans, 
                      method = "pairwise", 
                      by = c("Degradation", "Scrambling")) %>% 
      as_tibble(),
    # Cohen's d effect size
    emmeans::eff_size(tau_GrpDegScr_emmeans,
                      by = c("Degradation", "Scrambling"),
                      edf = median_df,
                      sigma = sigma(lmerTest::get_model(tau_GrpDegTtScr_lmerstep))) %>% 
      as_tibble() %>% 
      select(-c(df, SE)),
    by = c("contrast", "Degradation", "Scrambling")) %>% 
    rename(Within1 = Degradation,
           Within2 = Scrambling),
  # Compare scrambling conditions (Intact vs Scrambled) within each Degradation x Group combination
  # Tests if scrambling affects tau differently across conditions
  full_join(
    # Contrast
    emmeans::contrast(tau_GrpDegScr_emmeans, method = "pairwise", 
                      by = c("Degradation", "Group")) %>% 
      as_tibble(),
    # Cohen's d effect size
    emmeans::eff_size(tau_GrpDegScr_emmeans,
                      by = c("Degradation", "Group"),
                      edf = median_df,
                      sigma = sigma(lmerTest::get_model(tau_GrpDegTtScr_lmerstep))) %>% 
      as_tibble() %>% 
      select(-c(df, SE)),
    by = c("contrast", "Degradation", "Group")) %>% 
    rename(Within1 = Degradation,
           Within2 = Group),
  # Compare degradation slopes between groups within each Scrambling combination
  # Tests if groups differ in how tau changes with degradation
  # Contrast
  emmeans::contrast(tau_GrpDegScr_emtrends, 
                    method = "pairwise", 
                    by = "Scrambling") %>% 
    as_tibble() %>% 
    mutate(Within1 = 1) %>% 
    rename(Within2 = Scrambling)) %>% 
  as.data.frame() %>% 
  # Apply FDR correction to all p-values
  do(correct_p_vals(., "p.value"))

##### Plot #####

# Calculate y-axis range for tau plot
# Determines appropriate limits based on the data range
y_range <- c(min(tau_GrpDegScr_emmeans %>% 
                   as_tibble() %>% 
                   pull(lower.CL), na.rm = TRUE),
             max(tau_GrpDegScr_emmeans %>% 
                   as_tibble() %>% 
                   pull(upper.CL), na.rm = TRUE))
y_range <- c(y_range[1] - 0.05, y_range[2] + 0.05)
# Generate evenly spaced breaks for y-axis
y_breaks <- pretty(y_range, n = 3)

# Create annotation text for tau plot showing degradation slopes and significance
# Formats slope coefficients and p-values for display in the figure
tau_GrpDegScr_annot <- tau_GrpDegScr_trend %>% 
  mutate(
    # Format slope to 3 decimal places
    slope_str = case_when(
      round(abs(Degradation.trend), 3) < .001 ~ "<.001",
      round(abs(Degradation.trend), 3) >= .001 ~ sprintf("%.3f", Degradation.trend)
    ),
    # Create annotation string for slope
    slope_str_m = case_when(
      grepl("^<", slope_str) ~ paste0('<"', sub("^<", "", slope_str), '"'),
      grepl("^=", slope_str) ~ paste0('=="', sub("^=", "", slope_str), '"'),
      TRUE                   ~ paste0('=="', slope_str, '"')
    ),
    # Format p-values: < .001, = .001, or exact value
    p_str = case_when(
      round(p.value_fdrcorrected, 3) < .001 ~ "<.001",
      round(p.value_fdrcorrected, 3) == .001 ~ "=.001",
      TRUE                                  ~ sprintf("%.3f", p.value_fdrcorrected)
    ),
    # Remove leading zero from p-values (e.g., "0.05" -> ".05")
    p_str = gsub("^0\\.", ".", p_str),
    # Create annotation string for p-value
    p_str_m = case_when(
      grepl("^<", p_str) ~ paste0('<"', sub("^<", "", p_str), '"'),
      grepl("^=", p_str) ~ paste0('=="', sub("^=", "", p_str), '"'),
      TRUE               ~ paste0('=="', p_str, '"')
    ),
    # Concatenate and add significance marker
    annot = paste0(
      Group,
      ": italic(beta)", slope_str_m,
      '~","~italic(p)', p_str_m,
      ifelse(sig_fdrcorrected == "*", '~"*"', ''))) %>%
  # Set annotation position: left side of plot, top, with vertical offset for each group
  mutate(x_pos = -Inf,
         y_pos = y_range[2],  # both bottom
         vjust = ifelse(Group == "YA", 0.5, 1.8),    # stack vertically
         hjust = -0.05)

# Create line plot showing tau as a function of degradation
# Two columns for Scrambling condition, one row
# Separate lines for each Group
tau_GrpDegScr_line <- 
  ggplot(tau_GrpDegScr_emmeans %>% 
           as_tibble(), 
         aes(x = Degradation, 
             y = emmean,
             colour = Group)) +
  # Create grid: columns = Scrambling
  facet_grid(cols = vars(Scrambling)) +
  geom_line(size = 1.5, aes(group = Group,
                            colour = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 2, color = "black") +
  geom_point(size = 2, aes(group = Group, colour = Group)) +
  # Add annotation text showing slope coefficients and p-values
  geom_text(
    data = tau_GrpDegScr_annot,
    aes(x = x_pos, y = y_pos, label = annot, hjust = hjust, 
        vjust = vjust, colour = Group),
    parse = TRUE,
    size = 3,
    show.legend = FALSE) +
  # Set axis labels and legend titles
  labs(x = factor_labels$Degradation$axis_label, 
       y = "Estimated Tau", 
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
tau_GrpDegScr_line

############## Analyse Bias (beta) and Boundary Separation (alpha) ##############

alphabeta_data_Grp <- driftdiffusion_data %>%
  select(-c("delta", "tau")) %>%
  pivot_longer(cols = c("alpha", "beta"),
               names_to = "Parameter",
               values_to = "ParamVal") %>% 
  group_by(Participant, Parameter, Group) %>%
  summarise(ParamVal = mean(ParamVal))

alphabeta_Grp_ttest <- full_join(
  alphabeta_data_Grp %>% 
    group_by(Parameter) %>% 
    rstatix::t_test(formula = ParamVal ~ Group,
                    paired = FALSE),
  alphabeta_data_Grp %>% 
    group_by(Parameter) %>% 
    rstatix::cohens_d(formula = ParamVal ~ Group,
                      paired = FALSE),
  by = join_by(Parameter, .y., group1, group2, n1, n2)
)

beta_Grp_oneway <- full_join(
  alphabeta_data_Grp %>% 
    filter(Parameter == "beta") %>% 
    group_by(Group) %>% 
    rstatix::t_test(formula = ParamVal ~ 1,
                    mu = 0.5),
  alphabeta_data_Grp %>% 
    filter(Parameter == "beta") %>% 
    group_by(Group) %>% 
    rstatix::cohens_d(formula = ParamVal ~ 1,
                      mu = 0.5),
  by = join_by(Group, .y., group1, group2, n)
) %>% 
  as.data.frame() %>% 
  do(correct_p_vals(., "p"))
