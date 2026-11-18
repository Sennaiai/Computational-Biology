# Development and verification

- Discuss consequential design decisions and acceptance criteria with the human.
- For each task, add or update appropriate tests and run relevant checks.
  Report commands, results, failures, skipped checks, and remaining risks.
- The shared Python check is a baseline. Add project-specific reference cases,
  pipeline integration tests, controlled dependencies, and reproducibility
  evidence before claiming scientific or production readiness.
- Add Python tests under `tests/`. The `tools/FeNNol-PMC` submodule is outside
  this repository's lint and test scope and requires its own validation.
- Use synthetic data for tests unless explicitly authorized otherwise. Protect
  genomic data and credentials. Use external AI review only on request and
  only with authorization for the information shared.
- Request human review before merging or releasing.
