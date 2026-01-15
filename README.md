# Code for Section 3 and Section 4 of “Conditional estimations for seamless phase II/III clinical trials involving multi-stage early stopping”

This repository contains simulation code and illustrative examples supporting the paper  
**“Conditional estimations for seamless phase II/III clinical trials involving multi-stage early stopping.”**

The repository is organized to correspond directly to the structure of the paper.

---

## Section 3: Simulation Studies

The following files contain the full simulation framework used in **Section 3 of the manuscript**:

- `simulation-normal.R`  
- `simulation-binary (logOR).R`  
- `simulation-survival (log HR).R`  

These scripts implement:
- Data generation under different outcome types (normal, binary, survival) under a seamless phase II/III design with multi-stage early stopping
- Estimation methods including **naive**, **stage2**, **CBAE-SI**, **CBAE-MI**, **CMUE-MLE**, **CMUE-ZERO**, and **RB**
- Parallel computation for large-scale Monte Carlo simulations

Each script is self-contained and corresponds to a specific outcome type discussed in Section 3.

---

## Section 4: Illustrative Examples

The following files are used for the illustrative example in **Section 4 of the manuscript**:

- `example.R`  
- `example-normal.csv`  

The file `example.R` provides step-by-step code demonstrating how to implement the proposed methods on a single simulated dataset.  
The file `example-normal.csv` provides example data used in these illustrative analyses.

---

## Simulation Script Inputs (Section 3)

### Common parameters (used across simulation scripts)
- `K` : number of experimental treatment arms  
- `J` : number of stages  
- `t` : information fractions at each stage  
- `unit_n` : per-arm sample size per stage (normal/binary outcomes)  
- `unit_D` : expected number of events per stage (survival outcomes)  
- `n` / `D_n` : cumulative sample size or event counts at each stage  
- `mc_num` / `mc_n` : Monte Carlo size for conditional estimation  
- `sim` : number of Monte Carlo replicates  
- `setting` : scenario indicator (`NULL`, `H1`, `peak`, `line`)  

### Distribution-specific parameters
- **Normal outcome**
  - `mu_control`, `mu_treatment`
  - `sigma_control`, `sigma_treatment`

- **Binary outcome (log odds ratio)**
  - `p_control`
  - `p_treatment`

- **Survival outcome (log hazard ratio)**
  - `lambda_C` : control hazard
  - `lambda_E` : treatment hazards
  - `shape` : Weibull shape parameter
  - `a_t` : accrual time
  - `C` : accrual rate

---

## Example Script Inputs (Section 4)

### Common inputs
- `select_index` : index of the experimental arm selected at the end of phase II  
- `K` : number of experimental treatment arms  
- `J` : number of stages  
- `up_bround` : upper stopping boundary (efficacy boundary)  
- `low_bround` : lower stopping boundary (futility boundary)  
- `n` : cumulative control-arm sample size at each stage  
- `sigma_control` / `sigma_treatment` : standard deviations for control and treatment arms  

### Method-specific inputs
- **CBAE-MI**
  - `initial` : initial value for the iterative bias correction (here taken as the trial-end MLE)
  - `max_iterations = 20` : maximum number of iterations
  - `tol = 0.001` : convergence tolerance

- **CBAE-SI**
  - `initial` : initial value for the single-step bias correction (here taken as the trial-end MLE)

- **CMUE-MLE / CMUE-ZERO**
  - `sub = "MLE"` : substitution option; choose `"MLE"` or `"ZERO"`
  - `end_stage` : stage at which the trial stops

- **RB_mc**
  - `s_all` : length `K + 1` vector  
    - first `K` elements: stage-1 score statistics for each experimental arm  
    - element `K + 1`: score statistic for the selected arm at the stopping stage  
  - `end_stage` : stage at which the trial stops  
  - `mc_num` : Monte Carlo sample size used in Step 1 of the RB procedure  

---
