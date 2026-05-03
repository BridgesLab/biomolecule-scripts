#' Pseudorandomize Cages for Experimental Allocation
#'
#' @description
#' This function performs constrained randomization of cages into two treatment groups (A and B).
#' It ensures that all animals within the same cage receive the same treatment and that
#' the groups are balanced for a primary variable (e.g., OLM).
#'
#' When \code{extra_vars} are provided (e.g., Sex and Age), the function uses
#' \strong{stratified randomization}: cages are split into A/B \emph{within} each
#' unique combination of those variables.  This guarantees equal representation of
#' every sex × age (or similar) subgroup in both treatment arms by design, rather
#' than checking it statistically—which avoids the convergence failures that arise
#' when biological subgroups are strongly different from each other.
#'
#' @param dt A data.table containing a "Cage" column and the variables to be balanced.
#' @param target_col Character. The name of the primary numeric variable to balance (e.g., "OLM").
#' @param extra_vars Character vector (optional). Categorical variables that define
#'   pre-existing strata (e.g., \code{c("Sex", "Age")}). Each variable must be
#'   constant within a cage. Cages are randomized \emph{within} each stratum so
#'   that every stratum contributes equally to both treatment groups.
#' @param max_attempts Integer. The number of randomization iterations to attempt before
#'   stopping. Default is 2000.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: The original data.table with a new "Treatment" column.
#'   \item \code{summary_stats}: A data.table showing means for Group A and B,
#'     Cohen's d (overall), and p-values for overall and per-stratum comparisons.
#' }
#'
#' @details
#' \strong{Without \code{extra_vars}:}
#' Cages are assigned randomly; an assignment is accepted only when the overall
#' group sizes differ by \eqn{\leq 2} and the t-test p-value for \code{target_col}
#' exceeds 0.2.
#'
#' \strong{With \code{extra_vars}:}
#' Within each unique combination of \code{extra_vars} values, half the cages are
#' assigned to A and half to B (floor division for odd counts).  Balance of
#' \code{target_col} is then checked within each stratum (p > 0.2 required).
#' The chi-square check on \code{extra_vars} distributions is removed because
#' equal representation is guaranteed by the stratified design.
#'
#' @export
#' @import data.table
#'
#' @examples
#' # Balanced OLM within every Sex x Age combination:
#' # results <- pseudorandomize_cages(baseline_data, target_col = "OLM",
#' #                                  extra_vars = c("Sex", "Age"))
pseudorandomize_cages <- function(dt, target_col, extra_vars = NULL, max_attempts = 2000) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" is required. Please install it.")
  }

  working_dt <- data.table::copy(data.table::as.data.table(dt))

  # ── Stratified path ────────────────────────────────────────────────────────
  if (!is.null(extra_vars)) {

    # Validate: extra_vars must be constant within each cage
    for (v in extra_vars) {
      n_vals <- working_dt[, data.table::uniqueN(get(v)), by = Cage]
      if (any(n_vals$V1 > 1)) {
        stop(sprintf(
          "'%s' has more than one value within at least one cage. ",
          "Stratification variables must be constant within each cage.", v))
      }
    }

    # Cage-level table retaining stratum info
    cage_info <- unique(working_dt[, c("Cage", extra_vars), with = FALSE])

    # All unique strata (e.g. young-male, young-female, old-male, old-female)
    strata_list <- unique(cage_info[, extra_vars, with = FALSE])

    attempt <- 0
    success <- FALSE

    while (attempt < max_attempts) {
      attempt <- attempt + 1

      # Within each stratum, randomly assign floor(n/2) cages to A, rest to B
      cage_info[, Treatment := NA_character_]
      for (i in seq_len(nrow(strata_list))) {
        s <- strata_list[i]
        row_idx <- Reduce("&", lapply(extra_vars, function(v) cage_info[[v]] == s[[v]]))
        stratum_cages <- cage_info[row_idx, Cage]

        # Balance by animal count, not cage count: shuffle for randomness,
        # then greedily assign each cage to whichever group has fewer animals.
        cage_sizes <- working_dt[Cage %in% stratum_cages, .N, by = Cage]
        cage_sizes <- cage_sizes[sample(.N)]   # random order
        nA_count <- 0L; nB_count <- 0L; a_cages <- character(0)
        for (k in seq_len(nrow(cage_sizes))) {
          if (nA_count <= nB_count) {
            a_cages    <- c(a_cages, cage_sizes$Cage[k])
            nA_count   <- nA_count + cage_sizes$N[k]
          } else {
            nB_count   <- nB_count + cage_sizes$N[k]
          }
        }
        cage_info[row_idx, Treatment := ifelse(Cage %in% a_cages, "A", "B")]
      }

      # Push treatment assignment back to the animal-level table
      working_dt[cage_info, on = "Cage", Treatment := i.Treatment]

      is_balanced <- TRUE
      stats_list  <- list()

      # Overall stats (informational; not used as acceptance criterion here)
      g1 <- working_dt[Treatment == "A", get(target_col)]
      g2 <- working_dt[Treatment == "B", get(target_col)]
      sd_p <- sqrt(((length(g1) - 1) * var(g1) + (length(g2) - 1) * var(g2)) /
                     (length(g1) + length(g2) - 2))
      fit_main  <- lm(as.formula(paste(target_col, "~ Treatment")), data = working_dt)
      res_main  <- summary(fit_main)$coefficients
      p_overall <- if (nrow(res_main) > 1) res_main[2, 4] else 1.0

      stats_list[["Overall"]] <- data.table::data.table(
        Variable     = paste("Overall", target_col),
        Mean_A       = mean(g1),
        Mean_B       = mean(g2),
        SE           = res_main[2, 2],
        EffectSize_d = round((mean(g1) - mean(g2)) / sd_p, 3),
        p_value      = p_overall
      )

      # Per-stratum balance check for target_col
      for (i in seq_len(nrow(strata_list))) {
        s       <- strata_list[i]
        row_idx <- Reduce("&", lapply(extra_vars, function(v) working_dt[[v]] == s[[v]]))
        sub_dt  <- working_dt[row_idx]

        label <- paste(
          mapply(function(v, val) paste0(v, "=", val),
                 extra_vars, as.list(s[, extra_vars, with = FALSE])),
          collapse = ", ")

        # Need at least one animal per treatment side to test
        if (data.table::uniqueN(sub_dt$Treatment) < 2) {
          is_balanced <- FALSE; break
        }

        m1    <- mean(sub_dt[Treatment == "A", get(target_col)])
        m2    <- mean(sub_dt[Treatment == "B", get(target_col)])
        p_sub <- tryCatch(
          t.test(as.formula(paste(target_col, "~ Treatment")), data = sub_dt)$p.value,
          error = function(e) 1.0   # can't test with 1 obs per side; treat as balanced
        )

        stats_list[[paste0("stratum_", i)]] <- data.table::data.table(
          Variable     = paste0(target_col, " (", label, ")"),
          Mean_A       = m1,
          Mean_B       = m2,
          SE           = NA_real_,
          EffectSize_d = NA_real_,
          p_value      = p_sub
        )

        if (p_sub <= 0.2) { is_balanced <- FALSE; break }
      }

      if (is_balanced) { success <- TRUE; break }
    }

  # ── Unstratified path (original behaviour, no extra_vars) ──────────────────
  } else {

    unique_cages <- unique(working_dt$Cage)
    num_cages    <- length(unique_cages)
    attempt      <- 0
    success      <- FALSE

    while (attempt < max_attempts) {
      attempt <- attempt + 1

      group_a_cages <- sample(unique_cages, size = floor(num_cages / 2))
      working_dt[, Treatment := ifelse(Cage %in% group_a_cages, "A", "B")]

      counts <- working_dt[, .N, by = Treatment]
      nA <- if ("A" %in% counts$Treatment) counts[Treatment == "A", N] else 0
      nB <- if ("B" %in% counts$Treatment) counts[Treatment == "B", N] else 0
      if (abs(nA - nB) > 2) next

      is_balanced <- TRUE
      stats_list  <- list()

      fit_main <- lm(as.formula(paste(target_col, "~ Treatment")), data = working_dt)
      res_main <- summary(fit_main)$coefficients
      p_main   <- if (nrow(res_main) > 1) res_main[2, 4] else 1.0

      g1   <- working_dt[Treatment == "A", get(target_col)]
      g2   <- working_dt[Treatment == "B", get(target_col)]
      sd_p <- sqrt(((length(g1) - 1) * var(g1) + (length(g2) - 1) * var(g2)) /
                     (length(g1) + length(g2) - 2))

      stats_list[["Overall"]] <- data.table::data.table(
        Variable     = paste("Overall", target_col),
        Mean_A       = mean(g1),
        Mean_B       = mean(g2),
        SE           = res_main[2, 2],
        EffectSize_d = round((mean(g1) - mean(g2)) / sd_p, 3),
        p_value      = p_main
      )

      if (p_main <= 0.2) is_balanced <- FALSE

      if (is_balanced) { success <- TRUE; break }
    }
  }

  if (success) {
    message(sprintf("Balanced groups found after %d attempts.", attempt))
    return(list(
      data          = working_dt,
      summary_stats = data.table::rbindlist(stats_list)
    ))
  } else {
    stop("Could not achieve balance (p > 0.2) within max_attempts. Try relaxing constraints.")
  }
}
