# Load hBayesDM package for drift diffusion modeling
# Install if not already installed
if(!("hBayesDM" %in% rownames(installed.packages()))){
  install.packages("hBayesDM")
}
library(hBayesDM)

# Run Drift Diffusion Model (DDM) analysis for all experimental conditions.
#
# This function fits a DDM using hBayesDM::choiceRT_ddm for each combination
# of Group x Degradation x Scrambling x TrialType. Model outputs and diagnostic
# plots (trace and posterior density) are saved to disk.
#
# Requires variables from Initialise.R in the global environment:
#   participants, factor_labels, derivatives_path, and plotting themes.
#
# Args:
#   dataset_path: Path to the data directory containing participant folders.
#     Each participant folder should contain a {Participant}_TrialData.tsv file.
#   output_folder_modifier: String used as the output folder name under
#     derivatives_path (e.g., "MaximalModel" creates derivatives/MaximalModel/).
#
# Returns:
#   A dataframe tracking the status of each condition (Completed/Skipped/Error).
run_ddm <- function(test_data, output_folder_modifier) {

  set.seed(1)

  ############## Read in behavioural data and organise it ############
  # Prepare test data for DDM analysis
  # Exclude trials with missing RT or RT <= 100ms (anticipatory responses)
  # Convert variables to appropriate formats for choiceRT_ddm
  ddm_data <- test_data %>%
    filter(!(is.na(RT) | RT <= 100)) %>%
    # Determine what the correct response should be based on trial type
    # "Old" trials should be responded to as "Old", all others as "New"
    mutate(CorrectResponse = ifelse(TrialType == "Old",
                                    "Old",
                                    "New"),
           # Convert variables to factors with specified levels
           Group = factor(Group,
                          levels = factor_labels$Group$levels),
           TrialType = factor(TrialType,
                              levels = factor_labels$TrialType$levels),
           Scrambling = factor(Scrambling,
                               levels = factor_labels$Scrambling$levels))

  ############## Run DDM for each condition ############

  # Define all condition levels for the 2 x 4 x 2 x 3 design
  groups <- unique(ddm_data$Group)              # YA, OA
  degradations <- unique(ddm_data$Degradation)  # 40, 70, 85, 95
  scramblings <- unique(ddm_data$Scrambling)    # Intact, Scrambled
  trialtypes <- unique(ddm_data$TrialType)      # Old, Similar, New

  # Create output directory for DDM model results if it doesn't exist
  ddm_output_path <- paste0(derivatives_path, output_folder_modifier, "/Models/")
  if (!dir.exists(ddm_output_path)) {
    dir.create(ddm_output_path, recursive = TRUE)
  }

  # Create output directory for diagnostic plots if it doesn't exist
  diagnostic_plot_path <- paste0(derivatives_path, output_folder_modifier, "/DiagnosticPlots/")
  if (!dir.exists(diagnostic_plot_path)) {
    dir.create(diagnostic_plot_path, recursive = TRUE)
  }

  # Initialize a dataframe to track which conditions have been run
  condition_tracker <- expand.grid(
    Group = groups,
    Degradation = degradations,
    Scrambling = scramblings,
    TrialType = trialtypes,
    stringsAsFactors = FALSE
  ) %>%
    mutate(ConditionID = row_number(),
           Status = "Pending")

  # Save condition tracker for reference
  write.csv(condition_tracker,
            file = paste0(derivatives_path, output_folder_modifier, "/ConditionTracker.csv"),
            row.names = FALSE)

  # Loop through all conditions and run DDM
  for(i in 1:nrow(condition_tracker)){

    # Get current condition parameters
    curr_group <- condition_tracker$Group[i]
    curr_deg <- condition_tracker$Degradation[i]
    curr_scr <- condition_tracker$Scrambling[i]
    curr_tt <- condition_tracker$TrialType[i]
    curr_id <- condition_tracker$ConditionID[i]

    # Print progress
    cat(sprintf("\n========== Condition %d of %d ==========\n", i, nrow(condition_tracker)))
    cat(sprintf("Group: %s, Degradation: %s%%, Scrambling: %s, TrialType: %s\n",
                curr_group, curr_deg, curr_scr, curr_tt))

    # Filter data for current condition
    curr_data <- ddm_data %>%
      filter(Group == curr_group,
             Degradation == curr_deg,
             Scrambling == curr_scr,
             TrialType == curr_tt)

    # Check if there is sufficient data for this condition
    n_trials <- nrow(curr_data)
    n_participants <- length(unique(curr_data$Participant))

    cat(sprintf("Number of trials: %d, Number of participants: %d\n",
                n_trials, n_participants))

    if(n_trials < 10 || n_participants < 2){
      cat("Insufficient data for this condition. Skipping.\n")
      condition_tracker$Status[i] <- "Skipped - Insufficient data"
      next
    }

    # Prepare data for choiceRT_ddm
    # Required columns: subjID, choice, RT
    # choice: 1 = lower boundary (incorrect), 2 = upper boundary (correct)
    # RT: response time in seconds
    ddm_input <- curr_data %>%
      mutate(
        # Create subject ID as character
        subjID = as.character(Participant),
        # Convert RT from milliseconds to seconds
        RT = RT / 1000,
        # Create choice variable: 2 = correct response, 1 = incorrect response
        choice = ifelse(Response == "Old", 1, 2)
      ) %>%
      as.data.frame()

    # Run DDM using choiceRT_ddm
    # Parameters:
    #   niter: total number of iterations per chain
    #   nwarmup: number of warmup iterations (discarded)
    #   nchain: number of MCMC chains
    #   ncore: number of CPU cores to use
    tryCatch({
      ddm_model <- choiceRT_ddm(
        data = ddm_input,
        ncore = 6,
        niter = 15000,
        nwarmup = 2000,
        nchain = 6,
        max_treedepth = 20
      )

      # Save the model output
      output_filename <- sprintf("DDM_%s_%s_%s_%s.RData",
                                 curr_group, curr_deg, curr_scr, curr_tt)
      save(ddm_model, file = paste0(ddm_output_path, output_filename))
      cat(sprintf("Model saved to: %s\n", output_filename))

      # Save diagnostic plots
      # Create base filename for plots
      plot_basename <- sprintf("DDM_%s_%s_%s_%s",
                               curr_group, curr_deg, curr_scr, curr_tt)

      # Setup data to plot diagnostic plots
      ddm_draws <- lapply(c("alpha", "beta", "delta", "tau"),
                          function(param_name){
                            data.frame(Parameter = param_name,
                                       Draw = 1:13000,
                                       ParamVal = ddm_model$parVals[[paste0("mu_", param_name)]]) %>%
                              group_by(Draw) %>%
                              mutate(Chain = 1:6)
                          })
      ddm_draws <- do.call("rbind", ddm_draws) %>%
        mutate(Chain = factor(Chain),
               Parameter = factor(Parameter,
                                  levels = c("alpha", "beta", "delta", "tau"),
                                  labels = c("Boundary Separation",
                                             "Bias",
                                             "Drift Rate",
                                             "Non-decision Time")))

      # Make and save trace plot (shows MCMC chain convergence)
      trace_plot <- ggplot(ddm_draws,
                           aes(x = Draw,
                               y = ParamVal,
                               color = Chain)) +
        geom_line() +
        facet_wrap(~Parameter, nrow = 2, scales = "free_y") +
        # Apply custom themes for consistent formatting
        y_axis_theme + x_axis_theme +
        blank_bg_theme + legend_theme + paper_facet_theme
      # Save
      trace_filename <- paste0(diagnostic_plot_path, plot_basename, "_trace.png")
      png(trace_filename, width = 1200, height = 800)
      print(trace_plot)
      dev.off()
      cat(sprintf("Trace plot saved to: %s\n", trace_filename))


      # Save posterior density plot (shows parameter distributions)
      postdens_plot <- ggdensity(data = ddm_draws,
                                 x = "ParamVal",
                                 color = "Chain",
                                 size = 1) +
        facet_wrap(~Parameter, nrow = 2, scales = "free") +
        x_axis_theme + y_axis_theme + blank_bg_theme +
        paper_facet_theme +
        theme(axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 16),
              strip.text.x = element_text(size = 30))

      # Save
      density_filename <- paste0(diagnostic_plot_path, plot_basename, "_density.png")
      png(density_filename, width = 1200, height = 800)
      print(postdens_plot)
      dev.off()
      cat(sprintf("Density plot saved to: %s\n", density_filename))

      condition_tracker$Status[i] <- "Completed"

    }, error = function(e){
      cat(sprintf("Error running DDM: %s\n", e$message))
      condition_tracker$Status[i] <- paste("Error:", e$message)
    })

    # Update condition tracker file after each iteration
    write.csv(condition_tracker,
              file = paste0(derivatives_path, output_folder_modifier, "/ConditionTracker.csv"),
              row.names = FALSE)
  }

  # Print final summary
  cat("\n========== DDM Analysis Complete ==========\n")
  cat(sprintf("Completed: %d\n", sum(condition_tracker$Status == "Completed")))
  cat(sprintf("Skipped: %d\n", sum(grepl("Skipped", condition_tracker$Status))))
  cat(sprintf("Errors: %d\n", sum(grepl("Error", condition_tracker$Status))))
  cat(sprintf("\nResults saved to: %s\n", ddm_output_path))

  return(condition_tracker)
}
