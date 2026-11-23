# Porous Shell

Page text follows.

```@raw html
<!-- VALIDATION_RECORD_START: porous_shell -->
```

**Validation record**: `benchmarks/records/porous_shell.toml` (Source commit `d7f75b66`, config SHA-256 `d3d57425a447`)

| Case | Metric | Grid `Nx` | Measured | Reference | Tolerance | Units | Status |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `porous_shell_t2.5_phi0.25` | Median shell depth of the rock markers against the 1-D shell depth (Nx=33) | 33 | `6.0000e+03` | `4.5900e+03` | `atol=1.50e+03, rtol=0.00e+00` | m | **PASS** |
| `porous_shell_t2.5_phi0.25` | Median shell depth of the rock markers against the 1-D shell depth (Nx=65) | 65 | `5.2500e+03` | `4.5900e+03` | `atol=7.50e+02, rtol=0.00e+00` | m | **PASS** |
| `porous_shell_t2.5_phi0.25` | Shell depth error ratio abs(d33 - d_ref) / abs(d65 - d_ref), first order in dx | 65 | `2.1364e+00` | `2.0000e+00` | `atol=5.00e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `porous_shell_t2.5_phi0.25` | Relative L2 error of the binned peak temperature rise against the 1-D reference (Nx=65) | 65 | `4.9815e-02` | `0.0000e+00` | `<= 1.00e-01` | dimensionless | **PASS** |
| `porous_shell_t2.5_phi0.25` | Scaled energy residual of the last thermal step (Nx=65) | 65 | `2.2975e-10` | `0.0000e+00` | `<= 1.00e-09` | dimensionless | **PASS** |

```@raw html
<!-- VALIDATION_RECORD_END: porous_shell -->
```
