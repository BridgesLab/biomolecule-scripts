# Cage-Based Pseudorandomization Tool

This repository contains an R function designed for the **pseudorandomization of mice** into experimental groups. It ensures that baseline metrics (like OLM scores, Body Weight, or Age) are balanced between Treatment A and Treatment B while strictly maintaining **cage integrity**.

## Why Use This?
In animal research, mice housed in the same cage are more similar to each other than to mice in other cages ("cage effects"). To maintain statistical independence:

1.  **Cage Integrity:** All mice in a single cage must be assigned to the same treatment group to avoid cross-contamination of treatment effects and to account for shared variance.
2.  **Balanced Baselines:** Groups should not have significant differences at the start of an experiment. This tool uses **constrained randomization**—iteratively testing random assignments until one is found where all group differences (including subgroups) have a $p > 0.2$.

---

##  Features
* **Automatic Cage-Grouping:** All animals with the same `Cage` ID are moved as a single unit.
* **Subgroup Balancing:** Ensures balance across categorical variables (e.g., Sex, Age).
* **Granular Analytics:** Reports Mean `target_col` (e.g., OLM) specifically for subgroups (e.g., "OLM in Males") to ensure no hidden biases.
* **Constraint Checking:** Automatically rejects iterations where group sizes differ by more than $\pm 2$ or p-values are $\leq 0.2$.
* **Metadata Preservation:** The final output contains all original columns (IDs, Genotypes, etc.).

---

## How to Use

Ensure you have the `data.table` package installed:
```r
install.packages("data.table")
library(data.table)
source("pseudorandomize_cages.R")

# Example: Balancing OLM scores across Cages, while checking for Sex and Age balance
results <- pseudorandomize_cages(
  dt = baseline_data, 
  target_col = "OLM", # the variable column you want balanced across groups
  extra_vars = c("Sex", "Age") #the extra groups you want balance over (beyond treatments A/B)
)
```

### Access the data with the new 'Treatment' column

```r
final_dt <- results$data
```

### View the statistical balance report

```r
print(results$summary_stats)
```

### Understanding the Output Table

| Column | Description |
| :--- | :--- |
| **Variable** | The specific group being tested (e.g., "Overall OLM" or "OLM in Females"). |
| **Mean_A / Mean_B** | The average score for that variable in each treatment group. |
| **SE** | Standard Error of the difference (derived from the linear model). |
| **EffectSize_d** | Cohen’s $d$. A value close to 0 indicates highly similar groups. |
| **p_value** | The probability that differences are due to chance. The tool requires $p > 0.2$. |

## Limitations & Tips

- Convergence: If you have few cages but many metadata variables (e.g., trying to balance OLM, Age, Sex, and Weight with only 8 cages), it may be mathematically impossible to find a $p > 0.2$ for every combination. If the function hits max_attempts, try removing one of the extra_vars.
- Reproducibility: To ensure you can recreate the exact same group assignments later, run set.seed(123) (or any number) before calling the function