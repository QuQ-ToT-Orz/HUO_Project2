# Mortality Risk Reflects Temporal Clustering of Physical Activity, Rather Than Frequency or Intensity: A Dispersion Index Approach Grounded in Hawkes Process Theory

[![R-CMD-check](https://github.com/QuQ-ToT-Orz/Project2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/QuQ-ToT-Orz/Project2/actions/workflows/R-CMD-check.yaml)

Analysis of physical activity temporal patterns using Hawkes processes and dispersion indices from NHANES accelerometer data, with application to all-cause mortality prediction.

## Overview

This project implements a novel framework for characterizing the temporal clustering ("burstiness") of physical activity using:

1. **Hawkes Process Models** - Self-exciting point processes that capture how physical activity events trigger subsequent activity
2. **Dispersion Index Analysis** - A computationally efficient alternative that quantifies deviation from Poisson (random) activity patterns

The key insight is that physical activity is not randomly distributed throughout the day but exhibits clustering behavior that may be independently associated with health outcomes.

## HUO (Hawkes-based Unfolding of Observed human behavior) package
```r
devtools::install_local("./HUO", force = TRUE, dependencies = TRUE)
```

## Key Metrics
### Hawkes Process Parameters
- **μ (mu)**: Baseline intensity - rate of spontaneous activity initiation
- **α (alpha)**: Excitation parameter - how much each event increases future event probability
- **β (beta)**: Decay rate - how quickly the excitation effect diminishes
- **n = α/β**: Branching ratio - average number of triggered events per spontaneous event

### Dispersion-Based Metrics
- **D**: Dispersion index (variance-to-mean ratio of event counts) are adjusted for circadian rhythm using 30-minute bin normalization to avoid confounding time-of-day effects with true clustering.
- **n = 1 - 1/√D**: Branching ratio derived from dispersion
- **μ\***: Immigration rate 

Interpretation:
- n > 0: Clustered/bursty activity (events trigger more events)
- n = 0: Poisson process (random timing)
- n < 0: Regular/periodic activity (more evenly spaced)

### Mortality (Survey Weighting)
All population-level estimates use NHANES survey weights and account for the complex survey design (stratified cluster sampling).

## Citation
[Citation to be added upon publication]