library(RWiener)

# Function to simulate accuracy and RT data based on DDM parameters.
#
# Simulates responses using the Wiener diffusion model (rwiener) for each
# experimental condition. The number of trials simulated per condition matches
# the comparison dataset, allowing direct comparison between model predictions
# and empirical data.
#
# Some DDM parameters may be fixed (collapsed) across certain experimental
# factors. For example, if "alpha" does not vary by Degradation, the same
# alpha value is used across all Degradation levels. The function automatically
# detects which parameters vary by which factors by checking for "collapsed"
# values in mean_params_df.
#
# Args:
#   mean_params_df: A dataframe containing parameter values for the DDM.
#     Required columns: Parameter (alpha, tau, beta, delta), Group,
#     Degradation, Scrambling, TrialType, and ParamVal (the parameter
#     values). When a parameter is fixed across a factor, that factor's
#     value should be "collapsed" in this dataframe.
#   comparison_data: The empirical data to match trial counts against.
#     Required columns: Group, Degradation, Scrambling, TrialType,
#     resp (response), and q (RT).
#
# Returns:
#   A list containing:
#     - simulations: A dataframe with both empirical ("Data") and simulated
#       ("Model") data in long format, including accuracy coding.


# Helper function to look up a single DDM parameter value for a condition.
#
# Each parameter may or may not vary by each experimental factor (Group,
# Degradation, Scrambling, TrialType). If a parameter is collapsed across
# a factor, the lookup uses "collapsed" instead of the actual condition
# level.
#
# Args:
#   param_name: Which DDM parameter to look up ("alpha", "tau", "beta", or
#     "delta").
#   param_variation: Dataframe indicating which factors each parameter varies
#     by (output of the variation detection step in simulate_DDM).
#   mean_params_df: The full parameter dataframe to look up values from.
#   grp, deg, scr, tt: The current condition's factor levels.
#
# Returns:
#   The parameter value (numeric scalar) for the given parameter and condition.
lookup_ddm_param <- function(param_name, param_variation,
                             mean_params_df,
                             grp, deg, scr, tt) {

  # Get the variation flags for this specific parameter
  variation <- param_variation %>%
    filter(Parameter == param_name)

  # For each factor, decide whether to look up the actual condition level
  # or "collapsed" (if the parameter doesn't vary by that factor)
  lookup_group <- ifelse(variation$varies_by_group,
                         as.character(grp), "collapsed")
  lookup_degradation <- ifelse(variation$varies_by_degradation,
                               as.character(deg), "collapsed")
  lookup_scrambling <- ifelse(variation$varies_by_scrambling,
                              as.character(scr), "collapsed")
  lookup_trialtype <- ifelse(variation$varies_by_trialtype,
                             as.character(tt), "collapsed")

  # Filter to the matching row and extract the parameter value
  mean_params_df %>%
    filter(Parameter == param_name,
           Group == lookup_group,
           Degradation == lookup_degradation,
           Scrambling == lookup_scrambling,
           TrialType == lookup_trialtype) %>%
    pull(ParamVal)
}


simulate_DDM <- function(mean_params_df,
                         comparison_data) {
  set.seed(1)

  # ---- Step 1: Determine which factors each parameter varies by ----
  # For each DDM parameter, check whether it has different values across
  # levels of each experimental factor, or is "collapsed" (fixed) across
  # that factor. This is used later to look up the correct parameter value
  # for each simulated condition.
  param_variation <- mean_params_df %>%
    group_by(Parameter) %>%
    summarise(
      varies_by_group = !any(Group == "collapsed"),
      varies_by_degradation = n_distinct(Degradation) > 1 &&
        !any(Degradation == "collapsed"),
      varies_by_scrambling = n_distinct(Scrambling) > 1 &&
        !any(Scrambling == "collapsed"),
      varies_by_trialtype = n_distinct(TrialType) > 1 &&
        !any(TrialType == "collapsed")
    )

  # ---- Step 2: Identify all unique conditions to simulate ----
  # Extract every unique combination of experimental factors from the
  # comparison data. Each row defines one condition to simulate.
  conditions <- comparison_data %>%
    distinct(Group, Degradation, Scrambling, TrialType)

  # Accumulator for simulation results across all conditions
  all_simulations <- c()

  # ---- Step 3: Simulate data for each condition ----
  for (row_idx in seq_len(nrow(conditions))) {

    # Extract the current condition's factor levels
    grp <- conditions$Group[row_idx]
    deg <- conditions$Degradation[row_idx]
    scr <- conditions$Scrambling[row_idx]
    tt  <- conditions$TrialType[row_idx]

    # Get the empirical trials for this condition.
    # The number of rows determines how many trials to simulate.
    current_data <- comparison_data %>%
      filter(Group == grp,
             Degradation == deg,
             Scrambling == scr,
             TrialType == tt)

    # ---- Step 3a: Look up DDM parameters for this condition ----
    # Each parameter is looked up individually because they may vary by
    # different subsets of factors (e.g., alpha varies by Group only,
    # but delta varies by Group and TrialType).
    alpha_val <- lookup_ddm_param("alpha", param_variation,
                                  mean_params_df,
                                  grp, deg, scr, tt)
    tau_val   <- lookup_ddm_param("tau", param_variation,
                                  mean_params_df,
                                  grp, deg, scr, tt)
    beta_val  <- lookup_ddm_param("beta", param_variation,
                                  mean_params_df,
                                  grp, deg, scr, tt)
    delta_val <- lookup_ddm_param("delta", param_variation,
                                  mean_params_df,
                                  grp, deg, scr, tt)

    # ---- Step 3b: Simulate responses ----
    # Use the Wiener diffusion model to generate simulated responses
    # and RTs, with the same number of trials as the empirical data.
    simulated <- rwiener(
      n = nrow(current_data),
      alpha = alpha_val,
      tau   = tau_val,
      beta  = beta_val,
      delta = delta_val
    )

    # ---- Step 3c: Combine empirical and simulated data ----
    # Rename columns to distinguish empirical ("Data_") from simulated
    # ("Model_") responses and RTs before combining.
    empirical_renamed <- current_data %>%
      dplyr::rename(Empirical_Response = Response,
                    Empirical_RT = RT)

    simulated_renamed <- simulated %>%
      dplyr::rename(Simulated_Response = resp,
                    Simulated_RT = q) %>% 
      mutate(Simulated_RT = Simulated_RT*1000,
             Simulated_Response = as.character(
               factor(Simulated_Response,
                      levels = c("lower", "upper"),
                      labels = c("Old", "New"))))
    
    combined <- cbind(empirical_renamed, simulated_renamed)

    # ---- Step 3d: Pivot to long format and compute accuracy ----
    # Reshape so that each row represents either a "Data" (empirical) or
    # "Model" (simulated) observation, with columns resp and RT.
    combined_long <- combined %>%
      pivot_longer(cols = starts_with(c("Empirical", "Simulated")),
                   names_to = c("EmpSim", ".value"),
                   names_sep = "_")

    # Code accuracy: recode trial type to the expected Wiener response
    # boundary (Old = "lower", Similar/New = "upper"), then check if the
    # actual response matches the expected boundary.
    combined_long <- combined_long %>%
      mutate(
        Acc = ifelse(CorrectResponse == Response, 1, 0)
      )

    all_simulations <- bind_rows(all_simulations, combined_long)
  }

  return(list(simulations = all_simulations))
}
