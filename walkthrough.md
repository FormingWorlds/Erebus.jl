
## Fix Architectural Disjoints (Constants Consolidation)

- **Constants Consolidation**: Removed shadowed globals (`dsubgrids`, `dsubgridt`, `hr_al`, `hr_fe`, `rcrust`, `psurface`, `dt_longest`, `DTmax`, `etamin`) from `src/constants.jl`. Threaded these variables properly into `save_state`, `simulation_loop`, and `compute_marker_properties!`.
- **Plasticity Loop Fix**: Discovered and fixed a major array-bounds bug where `iplast` in `simulation_loop` was mistakenly iterating to the global timestep iteration limit (`cfg.solver.titermax`) instead of the viscoplastic iteration limit (`cfg.solver.nplast`).
- **Tests Re-architected**: Tests in `test_numerics.jl`, `test_physics.jl`, and `test_particles.jl` were updated to define necessary physical constants locally that were stripped from global space, isolating tests from internal global dependencies.
- **Review Verification**: The `/rip` skill successfully ran and identified invalid struct accesses (`cfg.solver.dxymax`, `cfg.thermal.hr_al`, and `cfg.solver.dphimax` as kwargs), which were patched before pushing. 
