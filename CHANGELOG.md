# Changelog

All notable changes to GNM-ToolBox will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-05-11

### Added
- Initial public release
- Core GNM pipeline: `gnm` → `gnm_fit` → `gnm_predict`
- NUFFT-accelerated kernel regression engine (`external/fast_nufft_reg/`)
- Two-level hierarchical model with constant and functional batch configurations
- Declarative formula interface
- GCV bandwidth selection with 1-SE rule
- Link functions: identity, log, logit, sqrt, probit
- MCD-based outlier detection with t-SNE visualization
- Out-of-sample prediction via stored griddedInterpolant objects
- Demo scripts for ABIDE I and HarMNqEEG datasets
- Paper figure generation scripts (fig01-fig07)
- Complete documentation and README

### Fixed
- Spline oscillation in `fpp_yhat` for batches with sparse data (Germany batch overflow bug): added automatic clamp in `nureg_y.m` that bounds grid values to 3× the data range with warning
