# Monte Carlo Statistical Simulations (MATLAB)

A small collection of **Monte Carlo simulations** written in MATLAB, comparing
**theoretical probability distributions** with **empirical results** obtained from simulation.

These projects focus on using simulation to validate probability models and
understand convergence behavior.

## Projects

### Coin Toss Simulation
- Simulates repeated fair coin tosses
- Compares empirical PMF/CDF with the theoretical **Binomial distribution**
- Analyzes run-lengths of heads and compares them to a **Geometric distribution**
- Includes visual comparisons (histograms, CDFs)

### Lottery Simulation (6/49)
- Computes the exact distribution of matches using the **Hypergeometric distribution**
- Runs large-scale Monte Carlo simulations (millions of trials)
- Compares simulated probabilities to theoretical values
- Includes confidence intervals to illustrate sampling error and convergence

## Structure
- `coin-toss-binomial/` — coin toss simulation and report
- `lottery-hypergeometric/` — lottery simulation and report

Each folder contains the MATLAB code and a short report describing the approach,
results, and conclusions.

## Running the code
Open the relevant `.m` file in MATLAB and run it.  
Plots are generated automatically by the scripts.

---

### Note
This project was done as part of a university course. Because these are very common assignments, I had to prioritize writing an original solution that wouldn’t trigger plagiarism checks, rather than fully optimizing style or structure. Also, the lottery project is the more advanced one.
