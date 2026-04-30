#' Pseudorandomize Cages for Experimental Allocation
#'
#' @description
#' This function performs constrained randomization of cages into two treatment groups (A and B).
#' It ensures that all animals within the same cage receive the same treatment and that 
#' the groups are balanced based on row counts and statistical parity (p > 0.2) for 
#' a primary variable (e.g., OLM) and optional metadata (e.g., Age, Sex).
#'
#' @param dt A data.table containing a "Cage" column and the variables to be balanced.
#' @param target_col Character. The name of the primary numeric variable to balance (e.g., "OLM").
#' @param extra_vars Character vector (optional). Additional categorical or numeric variables 
#'   to balance and report (e.g., c("Age", "Sex")).
#' @param max_attempts Integer. The number of randomization iterations to attempt before 
#'   stopping. Default is 2000.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: The original data.table with a new "Treatment" column.
#'   \item \code{summary_stats}: A data.table showing means for Group A and B, 
#'   standard error, Cohen's d (for numeric), and p-values for overall and subgroup comparisons.
#' }
#' 
#' @details 
#' The function uses a "while" loop to iterate through random assignments. An assignment 
#' is accepted only if:
#' 1. The difference in total N between groups is <= 2.
#' 2. The p-value for the primary variable and all subgroup comparisons is > 0.2.
#' 3. Categorical distributions (e.g., Sex ratio) are balanced (Chi-sq p > 0.2).
#'
#' @export
#' @import data.table
#'
#' @examples
#' # my_results <- pseudorandomize_cages(my_data, target_col = "OLM", extra_vars = c("Age", "Sex"))
pseudorandomize_cages <- function(dt, target_col, extra_vars = NULL, max_attempts = 2000) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" is required. Please install it.")
  }
  
  # 1. Setup - Preserve all original metadata
  working_dt <- data.table::copy(data.table::as.data.table(dt))
  unique_cages <- unique(working_dt$Cage)
  num_cages <- length(unique_cages)
  
  attempt <- 0
  success <- FALSE
  
  while (attempt < max_attempts) {
    attempt <- attempt + 1
    
    # 2. Randomly assign Cages (forcing cage-mates together)
    group_a_cages <- sample(unique_cages, size = floor(num_cages / 2))
    working_dt[, Treatment := ifelse(Cage %in% group_a_cages, "A", "B")]
    
    # 3. Check Row Count Constraint (+/- 2)
    counts <- working_dt[, .N, by = Treatment]
    nA <- if("A" %in% counts$Treatment) counts[Treatment == "A", N] else 0
    nB <- if("B" %in% counts$Treatment) counts[Treatment == "B", N] else 0
    
    if (abs(nA - nB) > 2) next
    
    # 4. Comprehensive Statistical Checks
    is_balanced <- TRUE
    stats_list <- list()
    
    # --- Primary Check: Overall Target (e.g. OLM) ---
    fit_main <- lm(as.formula(paste(target_col, "~ Treatment")), data = working_dt)
    res_main <- summary(fit_main)$coefficients
    p_main <- if(nrow(res_main) > 1) res_main[2, 4] else 1.0
    
    g1 <- working_dt[Treatment == "A", get(target_col)]
    g2 <- working_dt[Treatment == "B", get(target_col)]
    
    # Pooled SD for Cohen's d
    sd_p <- sqrt(((length(g1)-1)*var(g1) + (length(g2)-1)*var(g2))/(length(g1)+length(g2)-2))
    eff_size <- (mean(g1) - mean(g2)) / sd_p
    
    stats_list[["Overall"]] <- data.table::data.table(
      Variable = paste("Overall", target_col),
      Mean_A = mean(g1),
      Mean_B = mean(g2),
      SE = res_main[2, 2],
      EffectSize_d = round(eff_size, 3),
      p_value = p_main
    )
    
    # Fail if overall p-value is too low
    if (p_main <= 0.2) is_balanced <- FALSE
    
    # --- Subgroup Checks: OLM split by Age, Sex, etc. ---
    if (is_balanced && !is.null(extra_vars)) {
      for (v in extra_vars) {
        
        # A) Check distribution of the factor itself (Chi-sq) 
        # Prevents e.g. all Males in Group A
        tab <- table(working_dt[[v]], working_dt$Treatment)
        if (suppressWarnings(chisq.test(tab)$p.value) <= 0.2) {
          is_balanced <- FALSE; break
        }
        
        # B) Check Target Mean within each Subgroup Level
        levels_v <- unique(working_dt[[v]])
        for (lvl in levels_v) {
          sub_dt <- working_dt[get(v) == lvl]
          
          # Skip if a subgroup doesn't have both treatments represented
          if (uniqueN(sub_dt$Treatment) < 2) { is_balanced = FALSE; break }
          
          m1 <- mean(sub_dt[Treatment == "A", get(target_col)])
          m2 <- mean(sub_dt[Treatment == "B", get(target_col)])
          p_sub <- t.test(get(target_col) ~ Treatment, data = sub_dt)$p.value
          
          stats_list[[paste0(v, "_", lvl)]] <- data.table::data.table(
            Variable = paste0(target_col, " in ", lvl),
            Mean_A = m1,
            Mean_B = m2,
            SE = NA,
            EffectSize_d = NA,
            p_value = p_sub
          )
          
          if (p_sub <= 0.2) { is_balanced <- FALSE; break }
        }
        if (!is_balanced) break
      }
    }
    
    if (is_balanced) {
      success <- TRUE
      break
    }
  }
  
  if (success) {
    message(sprintf("Balanced groups found after %d attempts.", attempt))
    return(list(
      data = working_dt,
      summary_stats = data.table::rbindlist(stats_list)
    ))
  } else {
    stop("Could not achieve balance (p > 0.2) within max_attempts. Try relaxing constraints.")
  }
}