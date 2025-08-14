---
title: Development
---

Development
========

## Testing strategy

Automated unit tests are found in the `test` directory. These tests check for
pipeline or scripts correctness, but do not run dependency programs such as
assemblers. They are run automatically as part of the CI pipeline. Once Aviary
is setup, tests can be run as below:

```bash
pixi run -e dev pytest -v
```

Integration tests are found in `test/test_integration.py`. These tests are
specific to the Centre for Microbiome Research server setup and check for
integration between Aviary and its dependency programs such as assemblers.
They are excluded from standard pytest runs and the CI pipeline, but can be run
as below:

```bash
pixi run -e dev pytest -v --run-expensive --run-qsub
```
