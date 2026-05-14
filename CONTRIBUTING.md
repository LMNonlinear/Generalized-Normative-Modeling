# Contributing to GNM-ToolBox

Thanks for your interest in contributing to GNM-ToolBox.

## How to Contribute

### Reporting Bugs

1. Check [existing issues](https://github.com/LMNonlinear/Generalized-Normative-Modeling/issues) to see if the bug is already reported
2. If not, open a new issue with:
   - MATLAB version and OS
   - Minimal reproducible example
   - Expected vs. actual behavior
   - Full error message and stack trace

### Suggesting Enhancements

Open a GitHub issue with the `enhancement` label. Describe:
- The problem you're trying to solve
- The proposed solution
- Alternative solutions considered

### Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes following the code style below
4. Add/update tests if applicable
5. Update CHANGELOG.md with your changes
6. Commit with clear messages: `git commit -m "feat: add X feature"`
7. Push to your fork: `git push origin feature/my-feature`
8. Open a pull request against `main`

## Code Style

### MATLAB

- Function names: `lowercase_with_underscores`
- Variable names: `camelCase` or `lowercase`
- Constants: `UPPERCASE`
- Every public function must have an H1 comment line describing its purpose
- Document all parameters and return values using MATLAB's standard help format
- Prefer vectorized operations over loops
- Use `set_defaults` for optional parameters

Example function header:
```matlab
function T_new = gnm_predict(mnhs_trained, path_csv_new, varargin)
% gnm_predict  Predict batch-corrected z-scores for new data.
%
%   T_new = gnm_predict(mnhs_trained, path_csv_new)
%   T_new = gnm_predict(mnhs_trained, path_csv_new, 'Name', Value, ...)
%
%   Inputs:
%     mnhs_trained - trained GNM struct from gnm_fit
%     path_csv_new - path to CSV file with new observations
%
%   Name-Value parameters:
%     'tag' - project identifier tag (default: 'predict')
%
%   Output:
%     T_new - table with predicted z-scores and harmonized data
```

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

- `feat:` new feature
- `fix:` bug fix
- `docs:` documentation changes
- `refactor:` code restructuring without behavior change
- `perf:` performance improvement
- `test:` adding or fixing tests
- `chore:` maintenance (dependencies, build, etc.)

Example: `fix(predict): clamp fpp_yhat values to prevent spline oscillation`

## Testing

Before submitting a pull request, run the demo scripts to verify nothing is broken:

```matlab
run('setup.m')
run('demo/test_general_pipeline.m')
run('demo/test_formula_pipeline.m')
run('demo/test_link_pipeline.m')
```

## Questions?

Open a GitHub discussion or contact the maintainers.
