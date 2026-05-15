# pseudospectral-band-unfolding

This repository contains MATLAB scripts used to generate data and figures for the paper:

> Band Unfolding via the Quadratic Pseudospectrum
> Christopher A. Bairnsfather, Ralph M. Kaufmann, Terry A. Loring, Alexander Cerjan
> https://arxiv.org/abs/2605.05423

The code is associated with the numerical examples in the paper. The main purpose of this repository is to make the scripts used for data generation and figure preparation publicly available.

## Repository contents

```
.
├── README.md
├── .gitignore
├── breathing_honeycomb/
│   ├── breathing_heatmap.m
│   ├── c6v_tci_hamil_pos.m
│   ├── data_over_triangle.m
│   ├── HTTTXQ_looping.m
│   ├── HTXQ_looping.m
│   ├── smoothed_HTXQ_looping.m
│   └── spatially_resolved_breathing_heatmap.m
├── fibonacci_and_ssh/
│   ├── fibonacci_two_values.m
│   ├── ssh_m_is_three_eigenvectors.m
│   ├── ssh_model_improved.m
│   └── terry_method_m_is_three.m
├── plotting_scripts/
│   ├── peanut_plot_data_from_eigenvectors.m
│   ├── plot_isosurface_on.m
│   ├── plot_state_on_lattice.m
│   ├── split_data_and_axes.m
│   └── standard_plot_settings.m
└── utilities/
    ├── alt_eigs.m
    ├── fibonacci_word.m
    ├── idos.m
    ├── maybe_create_parpool.m
    ├── maybe_create_progress_bar.m
    ├── ssh_model_exact_solution.m
    ├── update_progress_bar.m
    └── validate_input.m
```

## Requirements

The scripts were used with MATLAB R2026a and the parallel computing toolbox. 

Some scripts may also require generated data files produced by other scripts in this repository.

## Notes

PLACEHOLDER_ADD_ANY_RELEVANT_NOTES_HERE

Examples of useful notes:

- The maybe_create_parpool script needs minor modification for use on Windows
- The progress bar will not appear if using MATLAB via the command line
- Generated data and images are not included in this repository

## Citation

If you use this code, please cite the associated paper:

```bibtex
PLACEHOLDER_BIBTEX_ENTRY_FOR_AFTER_ACCEPTANCE
```

## Contact

For questions about the code, contact:

Christopher A. Bairnsfather
cbairnsf at purdue dot edu
