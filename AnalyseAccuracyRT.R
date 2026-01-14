# Set working directory and source initialization file
setwd("~/GitDir/CodeWithPapers/PatternCompletion_CueFeatures/")
source("Initialise.R")

#Setup common properties for plots
MeanProperties <- list("Line" = list("Size"  = 1.5),
                       "Errorbar" = list("Width" = 0.2,
                                         "Size" =  0.7),
                       "Point" = list("Size" = 4,
                                      "Shape" = 21))

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

test_data <- behavioural_data %>% 
  filter(Phase == "Test")
test_data <- test_data %>% 
  mutate(CorrectResponse = ifelse(TrialType == "Old",
                                  "Old", 
                                  "New"),
         Accuracy = CorrectResponse == Response,
         Group = factor(Group,
                        levels = factor_labels$Group$levels),
         TrialType = factor(TrialType,
                            levels = factor_labels$TrialType$levels),
         Scrambling = factor(Scrambling,
                             levels = factor_labels$Scrambling$levels))

total_trials_GrpDegTtScr <- test_data %>% 
  group_by(Participant, Group, TrialType, 
           Degradation, Scrambling) %>% 
  summarise(TotalTrials = length(Scrambling))

############## Analyse Test Data ############

########## Trials without a response ########

bad_trials <- test_data %>% 
  filter(is.na(RT) | RT <= 100)

badtrials_GrpDegTtScr <- bad_trials %>% 
  group_by(Participant, Group, TrialType, 
           Degradation, Scrambling) %>%
  summarise(BadTrials = length(Scrambling)) %>% 
  ungroup() %>% 
  full_join(total_trials_GrpDegTtScr,
            by = c("Participant", "Group", "TrialType", 
                   "Degradation", "Scrambling")) %>% 
  mutate(BadTrials = ifelse(is.na(BadTrials),
                            0,
                            BadTrials),
         PercBadTrials = (BadTrials/TotalTrials)*100)

badtrials_GrpDegTtScr_anova <- anova_test(data = badtrials_GrpDegTtScr,
                                          dv = PercBadTrials,
                                          wid = Participant,
                                          between = Group,
                                          within = c(TrialType, Degradation, Scrambling),
                                          type = 3)

badtrials_GrpDeg <- bad_trials %>% 
  group_by(Participant, Group, Degradation) %>%
  summarise(BadTrials = length(Scrambling)) %>% 
  ungroup() %>% 
  full_join(test_data %>% 
              group_by(Participant, Group, Degradation) %>% 
              summarise(TotalTrials = length(Scrambling)),
            by = c("Participant", "Group", "Degradation")) %>% 
  mutate(BadTrials = ifelse(is.na(BadTrials),
                            0,
                            BadTrials),
         PercBadTrials = (BadTrials/TotalTrials)*100)

badtrials_GrpDeg_posthoc <- badtrials_GrpDeg %>% 
  group_by(Group) %>% 
  t_test(formula = PercBadTrials ~ Degradation,
         paired = T)

########## Accuracy ########

acc_data_GrpDegTtScr <- test_data %>% 
  filter(!(is.na(RT) | RT <= 100)) %>% 
  group_by(Participant, Group, TrialType, 
           Degradation, Scrambling) %>%
  summarise(NumCorr = sum(Accuracy),
            NumIncorr = sum(!Accuracy)) %>% 
  full_join(total_trials_GrpDegTtScr,
            by = c("Participant", "Group", 
                   "TrialType", "Degradation", 
                   "Scrambling")) %>% 
  mutate(PercAcc = (NumCorr/TotalTrials)*100)

contrasts(acc_data_GrpDegTtScr$Scrambling) <- rbind(-1, 1)
colnames(contrasts(acc_data_GrpDegTtScr$Scrambling)) <- 
  levels(acc_data_GrpDegTtScr$Scrambling)[2]

acc_GrpDegTtScr_lmer <- lmerTest::lmer(
  formula = PercAcc ~ Degradation * TrialType * Scrambling * Group +
    (1|Participant), 
  data = acc_data_GrpDegTtScr)

acc_GrpDegTtScr_lmerstep <- lmerTest::step(acc_GrpDegTtScr_lmer, 
                                           direction = "backward")
acc_GrpDegTtScr_lmerfinal <- lmerTest::get_model(acc_GrpDegTtScr_lmerstep)

##### Group x Degradation x Trial Type interaction #####

acc_GrpDegTt_emmeans <- emmeans::emmeans(
  acc_GrpDegTtScr_lmerfinal, 
  ~ Group * TrialType * Degradation,
  at = list(Degradation = c(40, 70, 85, 95)))
acc_GrpDegTt_chance <- emmeans::test(acc_GrpDegTt_emmeans, 
                                     null = 50, 
                                     side = 0)
acc_GrpDegTt_chance <- correct_p_vals(acc_GrpDegTt_chance, 
                                      "p.value")

acc_GrpDegTt_emtrends <- emmeans::emtrends(
  acc_GrpDegTtScr_lmerfinal, 
  ~ Group * TrialType, var = "Degradation")
acc_GrpDegTt_trend <- emmeans::test(acc_GrpDegTt_emtrends, 
                                    null = 0, 
                                    side = 0)
acc_GrpDegTt_trend <- correct_p_vals(acc_GrpDegTt_trend, 
                                     "p.value")

median_df <- median(acc_GrpDegTt_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

acc_GrpDegTt_contrast <- bind_rows(
  full_join(
    emmeans::contrast(acc_GrpDegTt_emmeans, 
                      method = "pairwise", 
                      by = c("TrialType", "Degradation")) %>% 
      as_tibble() %>% 
      rename(Within1 = TrialType,
             Within2 = Degradation),
    emmeans::eff_size(acc_GrpDegTt_emmeans,
                      by = c("TrialType", "Degradation"),
                      sigma = sigma(lmerTest::get_model(acc_GrpDegTtScr_lmerstep)),
                      edf = median_df) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = TrialType,
             Within2 = Degradation),
    by = c("contrast", "Within1", "Within2")),
  emmeans::contrast(acc_GrpDegTt_emtrends, 
                    method = "pairwise", 
                    by = "TrialType") %>% 
    as_tibble() %>% 
    rename(Within1 = TrialType) %>% 
    mutate(Within2 = NA)) %>% 
  as.data.frame() %>% 
  do(correct_p_vals(., "p.value"))

acc_GrpDegTt_annot <- acc_GrpDegTt_trend %>% 
  mutate(
    slope_str = sprintf("%.2f", Degradation.trend),
    p_str = case_when(
      round(`p.value_fdrcorrected`, 3) < .001 ~ "<.001",
      round(`p.value_fdrcorrected`, 3) == .001 ~ "=.001",
      round(`p.value_fdrcorrected`, 3) > .001 ~ sprintf("%.3f", `p.value_fdrcorrected`)),
    p_str = gsub("^0\\.", ".", p_str),
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
  mutate(x_pos = -Inf,
         y_pos = acc_GrpDegTt_emmeans %>% 
           as_tibble() %>% 
           pull(emmean) %>% 
           min(),
         vjust = ifelse(Group == "YA", 1.3, 2.6),
         hjust = -0.05)

acc_GrpDegTt_line <- 
  ggplot(acc_GrpDegTt_emmeans %>% 
           as_tibble(), 
         aes(x = Degradation, 
             y = emmean,
             colour = Group)) +
  facet_grid(cols = vars(TrialType),
             scales = "free") +
  geom_line(size = 1.5,
            aes(group = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 2, color = "black") +
  geom_point(size = 3, aes(group = Group, colour = Group)) +
  geom_text(
    data = acc_GrpDegTt_annot,
    aes(x = x_pos, y = y_pos, label = annot, hjust = hjust, 
        vjust = vjust, colour = Group),
    parse = TRUE,
    size = 4,
    show.legend = FALSE) +
  labs(x = factor_labels$Degradation$axis_label, 
       y = "Estimated Accuracy", 
       fill = factor_labels$Group$axis_label, 
       color = factor_labels$Group$axis_label) +
  scale_fill_manual(values = factor_labels$Group$colours) +
  scale_color_manual(values = factor_labels$Group$colours) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), 
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = c(40, 70, 85, 95),
                     expand = expansion(mult = c(0.1, 0.1))) +
  coord_cartesian(ylim = c(25, 102)) +
  geom_hline(yintercept = 50, linetype = "dashed") + 
  y_axis_theme + x_axis_theme +
  blank_bg_theme + legend_theme + paper_facet_theme +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20))

##### Group x Trial Type x Scrambling interaction #####

acc_GrpTtScr_emmeans <- emmeans::emmeans(
  acc_GrpDegTtScr_lmerfinal, 
  ~ Group * TrialType * Scrambling)
acc_GrpTtScr_chance <- emmeans::test(acc_GrpTtScr_emmeans, 
                                     null = 50, 
                                     side = 0)
acc_GrpTtScr_chance <- correct_p_vals(acc_GrpTtScr_chance, 
                                      "p.value")

median_df <- median(acc_GrpDegTt_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))

acc_GrpTtScr_contrast <- bind_rows(
  full_join(
    emmeans::contrast(acc_GrpTtScr_emmeans, 
                      method = "pairwise", 
                      by = c("TrialType", "Scrambling")) %>% 
      as_tibble() %>% 
      rename(Within1 = TrialType,
             Within2 = Scrambling),
    emmeans::eff_size(acc_GrpTtScr_emmeans, 
                      sigma = sigma(lmerTest::get_model(acc_GrpDegTtScr_lmerstep)), 
                      by = c("TrialType", "Scrambling"),
                      edf = median_df) %>% 
      as_tibble() %>% 
      select(-c(df, SE)) %>% 
      rename(Within1 = TrialType,
             Within2 = Scrambling),
    by = c("contrast", "Within1", "Within2")),
  # Compare displacements within each Group
  full_join(
    emmeans::contrast(acc_GrpTtScr_emmeans, method = "pairwise", by = c("TrialType", "Group")) %>%
      as_tibble() %>%
      rename(Within1 = TrialType,
             Within2 = Group),
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
  do(correct_p_vals(., "p.value"))

acc_GDTt.predline <- 
  ggplot(acc_GrpTtScr_emmeans %>% 
           as_tibble(), 
         aes(x = Scrambling, 
             y = emmean,
             colour = Group)) +
  facet_grid(cols = vars(TrialType),
             scales = "free") +
  geom_line(size = 1.5,
            aes(group = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.08, color = "black") +
  geom_point(size = 3, aes(group = Group, colour = Group)) +
  labs(x = factor_labels$Scrambling$axis_label, 
       y = "Estimated Accuracy", 
       fill = factor_labels$Group$axis_label, 
       color = factor_labels$Group$axis_label) +
  scale_fill_manual(values = factor_labels$Group$colours) +
  scale_color_manual(values = factor_labels$Group$colours) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), 
                     expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(ylim = c(25, 102)) +
  geom_hline(yintercept = 50, linetype = "dashed") + 
  y_axis_theme + x_axis_theme +
  blank_bg_theme + legend_theme + paper_facet_theme +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

ggarrange(acc_GrpDegTt_line,
          acc_GDTt.predline +
            theme(axis.title.x = element_blank()),
          nrow = 2, ncol = 1, 
          common.legend = TRUE, 
          legend = "bottom")

########## Reaction Time ##########

rt_data_GrpDegTtScr <- test_data %>%
  filter(!(is.na(RT) | RT <= 100)) %>% 
  group_by(Participant, Group, TrialType, Degradation, 
           Scrambling) %>%
  summarise(MeanRT=mean(RT))

contrasts(rt_data_GrpDegTtScr$Scrambling) <- rbind(-1, 1)
colnames(contrasts(rt_data_GrpDegTtScr$Scrambling)) <- 
  levels(rt_data_GrpDegTtScr$Scrambling)[2]

rt_GrpDegTtScr_lmer <- lmerTest::lmer(
  formula = MeanRT ~ Degradation * TrialType * Scrambling * Group +
    (1|Participant), 
  data = rt_data_GrpDegTtScr)

rt_GrpDegTtScr_lmerstep <- lmerTest::step(rt_GrpDegTtScr_lmer, 
                                          direction = "backward")
rt_GrpDegTtScr_lmerfinal <- lmerTest::get_model(rt_GrpDegTtScr_lmerstep)

##### Group x Degradation x Trial Type x Scrambling #####

rt_GrpDegTtScr_emmeans <- emmeans::emmeans(rt_GrpDegTtScr_lmerfinal, 
                                           ~ Group * Scrambling * Degradation * TrialType,
                                           at = list(Degradation = c(40, 70, 85, 95)))

rt_GrpDegTtScr_emtrends <- emmeans::emtrends(rt_GrpDegTtScr_lmerfinal, 
                                             ~ Group * Scrambling * TrialType,
                                             var = "Degradation")
rt_GrpDegTtScr_trend <- emmeans::test(rt_GrpDegTtScr_emtrends, 
                                      null = 0, 
                                      side = 0)
rt_GrpDegTtScr_trend <- correct_p_vals(rt_GrpDegTtScr_trend, 
                                       "p.value")

median_df <- median(rt_GrpDegTtScr_emmeans %>% 
                      as_tibble() %>% 
                      pull(df))
rt_GrpDegTtScr_contrast <- bind_rows(
  full_join(
    emmeans::contrast(rt_GrpDegTtScr_emmeans, 
                      method = "pairwise", 
                      by = c("Degradation", "TrialType", "Scrambling")) %>% 
      as_tibble() %>% 
      rename(Within1 = Degradation,
             Within2 = TrialType,
             Within3 = Scrambling),
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
  full_join(
    emmeans::contrast(rt_GrpDegTtScr_emmeans, method = "pairwise", 
                      by = c("Degradation", "TrialType", "Group")) %>% 
      as_tibble() %>% 
      rename(Within1 = Degradation,
             Within2 = TrialType,
             Within3 = Group),
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
  emmeans::contrast(rt_GrpDegTtScr_emtrends, 
                    method = "pairwise", 
                    by = c("TrialType", "Scrambling")) %>% 
    as_tibble() %>% 
    rename(Within2 = TrialType,
           Within3 = Scrambling) %>% 
    mutate(Within1 = 1)) %>% 
  as.data.frame() %>% 
  do(correct_p_vals(., "p.value"))

y_range <- range(rt_GrpDegTtScr_emmeans %>% 
                   as_tibble() %>% 
                   pull(emmean), 
                 na.rm = TRUE)
y_range <- c(y_range[1] - 100, y_range[2] + 100)
y_breaks <- pretty(y_range, n = 3)
rt_GrpDegTtScr_annot <- rt_GrpDegTtScr_trend %>% 
  mutate(
    slope_str = sprintf("%.2f", Degradation.trend),
    p_str = case_when(
      round(`p.value_fdrcorrected`, 3) < .001 ~ "<.001",
      round(`p.value_fdrcorrected`, 3) == .001 ~ "=.001",
      round(`p.value_fdrcorrected`, 3) > .001 ~ sprintf("%.3f", `p.value_fdrcorrected`)
    ),
    p_str = gsub("^0\\.", ".", p_str),
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
  mutate(x_pos = -Inf,
         y_pos = y_range[2],  # both bottom
         vjust = ifelse(Group == "YA", 0.5, 1.8),    # stack vertically
         hjust = -0.05)

# Get unique trialtypes
trialtypes <- unique(rt_GrpDegTtScr_emmeans %>% 
                       as_tibble() %>% 
                       pull(TrialType))

# Create and save plots for each trialtype
rt_GrpDegTtScr_line <- 
  ggplot(rt_GrpDegTtScr_emmeans %>% 
           as_tibble(), 
         aes(x = Degradation, 
             y = emmean,
             colour = Group)) +
  geom_line(size = 1.5, aes(group = Group,
                            colour = Group)) +
  facet_grid(rows = vars(TrialType),
             cols = vars(Scrambling)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 2, color = "black") +
  geom_point(size = 2, aes(group = Group, colour = Group)) +
  geom_text(
    data = rt_GrpDegTtScr_annot,
    aes(x = x_pos, y = y_pos, label = annot, hjust = hjust, 
        vjust = vjust, colour = Group),
    parse = TRUE,
    size = 3,
    show.legend = FALSE) +
  # Add customizations
  labs(x = factor_labels$Percent$axis_label, 
       y = "Estimated RT", 
       fill = factor_labels$Group$axis_label, 
       color = factor_labels$Group$axis_label) +
  scale_fill_manual(values = factor_labels$Group$colours) +
  scale_color_manual(values = factor_labels$Group$colours) +
  scale_y_continuous(limits = y_range,
                     breaks = y_breaks,
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = c(40, 70, 85, 95),
                     expand = expansion(mult = c(0.1, 0.1))) +
  y_axis_theme + x_axis_theme +
  blank_bg_theme + legend_theme + paper_facet_theme +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

