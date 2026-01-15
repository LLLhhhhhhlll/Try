# Try

This repository contains simulation code and examples supporting the numerical studies in our manuscript on multi-arm multi-stage (MAMS) designs with treatment selection and early stopping.

The repository is organized to correspond directly to the structure of the paper.

---

## Repository Structure



---

## Relation to the Manuscript

### Section 3: Simulation Studies

The files  
- `simulation-normal.R`  
- `simulation-binary (logOR).R`  
- `simulation-survival (log HR).R`  

contain the full simulation framework used in **Section 3 of the manuscript**.

These scripts implement:
- Data generation under different outcome types (normal, binary, survival)
- Multi-arm multi-stage designs with early stopping
- Treatment selection at interim analyses
- Estimation methods including naive, stage-2, SI, MI, MUE, UMVCUE, and RB estimators
- Parallel computation for large-scale Monte Carlo simulations

Each script is self-contained and corresponds to a specific outcome type discussed in Section 3.

---

### Section 4: Illustrative Examples

The folder  
- `example/`

contains example code used in **Section 4 of the manuscript** to illustrate:
- Step-by-step implementation of the proposed methods
- Practical usage on a single simulated dataset
- Interpretation of estimators and confidence intervals

The file `example-normal.csv` provides example data used in these illustrative analyses.

---

## Input Parameters (Input List)

Across the simulation scripts, the main input parameters include:

### General Design Parameters
- `K` : number of experimental treatment arms  
- `J` : number of stages  
- `t` : information fractions at each stage  
- `alpha` : one-sided type I error level  
- `aftype` : alpha-spending function type (e.g., `sfLDOF`)  

### Sample Size / Information Parameters
- `unit_n` : per-arm sample size per stage (normal/binary outcomes)  
- `unit_D` : expected number of events per stage (survival outcomes)  
- `n` / `D_n` : cumulative sample size or event counts at each stage  

### Outcome-Specific Parameters
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

### Simulation Parameters
- `sim` : number of Monte Carlo replicates  
- `setting` : scenario indicator (`NULL`, `H1`, `peak`, `line`)  
- `mc_num` / `mc_n` : Monte Carlo size for conditional estimation  

---

## How to Run

Each simulation script can be run directly in R, for example:

```r
source("simulation-normal.R")
