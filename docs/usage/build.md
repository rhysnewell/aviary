---
title: aviary build
---

# aviary build

Build the Pixi environments used by Aviary's workflow rules, then exit.

```bash
aviary build
```

Unlike the analysis commands, `build` does not accept the shared performance,
output or Snakemake options.

## Examples

Build the CPU environments:

```bash
aviary build
```

Build CPU and GPU environments:

```bash
aviary build --gpu
```

The analysis commands also expose `--build` and `--build-gpu`, which build
environments and exit before starting that command's workflow. Prefer the
dedicated `aviary build` command when no analysis invocation is otherwise
needed.

## Options

**`--gpu`** `[BOOL]`

  Also build GPU-enabled environments. A bare `--gpu` means `true`; an
  explicit boolean value is also accepted. A compatible GPU must be present.
  [default: false]

**`-w`**, **`--workflow`** TARGET [TARGET ...]

  Snakemake target to run. This is exposed for parser consistency; the default
  and intended target is `build`. [default: `build`]
