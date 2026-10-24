# Symmetry-resolved Hubbard entanglement

Data and analysis scripts for

> G. Bellomia, C. Mejuto-Zaera, M. Capone and A. Amaricci,
> *Symmetry-resolved quasilocal entanglement as a key to Mott and pseudogap physics in plaquette CDMFT*

The repository contains the zero-temperature 2×2 plaquette CDMFT/ED solutions of the two-dimensional Hubbard model discussed in the paper, the CRISB reference data of Appendix C, and the MATLAB scripts that compute all entanglement measures and draw the figures.

## Contents

```
half-filled/            CDMFT at half filling, scan in U          (Sec. 4.1, Fig. 5, Figs. 2-3, 9, 11)
half-filled/CRISBreference/   CRISB at half filling, scan in U    (Appendix C, Fig. 11)
doped/                  CDMFT at U = 2.3D, scan in μ              (Sec. 4.2, Figs. 6-8)
scripts/                MATLAB scripts and functions
lib/qcmp-lab/           frozen copy of the QcmP-lab toolbox (data I/O and colors)
```

## Quick start

Requirements: MATLAB (tested with R2025b) with no additional toolboxes.

```matlab
cd scripts
half_sym_entanglement      % e.g. draws Fig. 5
```

Every script must be run from inside `scripts/`: it adds `../lib/qcmp-lab` to the path by itself. Scripts that read data move into the relevant data folder and return to `scripts/` at the end; if one stops with an error, `cd` back to `scripts/` before running it again.

While reading the data, the scripts (re)write the summary `*.txt` tables at the top of the data folders (see [Summary tables](#summary-tables)); the versions shipped here are the same files.

## Figures and scripts

| Paper | Script | Data |
|---|---|---|
| Figs. 2, 3, 9 | `spy_symmetry_resolution.m` | analytic Hubbard dimer + `half-filled/` |
| Fig. 4 | `analytic_dimer.m` | none (exact isolated dimer) |
| Eq. 27 | `noninteracting.m` | none (noninteracting square lattice) |
| Fig. 5 | `half_sym_entanglement.m` | `half-filled/` |
| Sec. 4.1, next-nearest neighbors | `nnn_sym_entanglement.m` | `half-filled/`, `doped/` |
| Fig. 6 | `doped_nonfreeness.m` | `doped/` |
| Fig. 7 | `doped_sym_entanglement.m` | `doped/` |
| Fig. 8 | `doped_dominant_states.m` | `doped/` |
| Fig. 11 | `crisb_vs_cdmft.m` | `half-filled/`, `half-filled/CRISBreference/`, `scripts/reduced_density_matrix_2sites_noninteracting.dat` |

Notes:

- `spy_symmetry_resolution.m` asks for a value of U among the computed ones. The sparsity patterns do not depend on it. It also draws a few intermediate orderings of the density matrix that are not shown in the paper.
- `noninteracting.m` prints the values of Eq. 27 and writes the nearest-neighbor density matrix of the noninteracting lattice to `scripts/`: `reduced_density_matrix_2sites_noninteracting.dat` (infinite lattice, used for the stars in Fig. 11) and `reduced_density_matrix_2sites_noninteracting_size.dat` (100×100 lattice). At L = 100 the script needs well over 10 GB of RAM; both files are already included.
- The manuscript figures were sometimes restyled in Inkscape after export; the inset of Fig. 7 was added there.

### Next-nearest neighbors (Sec. 4.1)

`nnn_sym_entanglement.m` computes the symmetry-resolved negativities and relative entropies of the two next-nearest-neighbor density matrices, i.e. the plaquette diagonals (1,4) and (2,3), at every computed point. It reports the points where a superselected sector is entangled, i.e. where its negativity exceeds $10^{-12}$. The two diagonals agree to three digits. At half filling:

| Quantity | Result |
|---|---|
| $\mathcal{N}^{\uparrow\downarrow}$, $E^{\uparrow\downarrow}$ | exactly 0 for all U |
| $\mathcal{N}^{\mathrm{hd}}$, $E^{\mathrm{hd}}$ | 0 up to U/D = 1.9, finite from U/D = 2.0 on |
| largest $\mathcal{N}^{\mathrm{hd}}$ | 4.70·10⁻⁴ bit, at U/D = 2.4 |
| largest $E^{\mathrm{hd}}$ | 7.45·10⁻⁵ bit, at U/D = 2.9 |
| largest $\mathcal{N}^{\mathrm{F}}$ | 9.45·10⁻³ bit, at U/D = 1.6 |

The hopping negativity $\mathcal{N}^{\mathrm{F}t}$ exceeds $\mathcal{N}^{\mathrm{hd}}$ at every U, but the holon-doublon share grows with U: at U/D = 8, $\mathcal{N}^{\mathrm{F}t} = 3.19·10^{-5}$ bit and $\mathcal{N}^{\mathrm{hd}} = 1.49·10^{-5}$ bit. The holon-doublon sector is entangled when $|\rho(6,11)|^2 > \rho(1,1)\,\rho(16,16)$ (Eq. 19); the ratio of the two sides grows smoothly with U, from 0.11 at U/D = 1.5 to 1.06 at U/D = 2.0 and about 470 at U/D = 8.

In the doped dataset (U/D = 2.3) the antiferromagnetic sector vanishes everywhere as well, while the holon-doublon sector is entangled only at μ/D = −0.48, −0.47 and −0.46, with $\mathcal{N}^{\mathrm{hd}} \le 1.42·10^{-4}$ bit and $E^{\mathrm{hd}} \le 1.02·10^{-5}$ bit.

### Shared functions in `scripts/`

| Function | Computes |
|---|---|
| `partial_transpose.m` | standard and fermionic partial transpose of a two-site density matrix (Eqs. 6-7) |
| `negativity.m` | full logarithmic negativity, standard and fermionic (Eqs. 4-8) |
| `sym_negativity.m` | symmetry-resolved negativities: spin-singlet, hopping and holon-doublon (Eqs. 12-13) |
| `sym_ree.m` | symmetry-resolved relative entropy of entanglement: spin-singlet and holon-doublon (Eq. 20) |
| `build_SSR.m` | relative entropy of entanglement under N-SSR and P-SSR (cross-check of `sym_ree.m`) |
| `vonNeumann.m` | von Neumann entropy |
| `slater_condon.m` | many-body matrix elements of one-body operators, used to get ⟨c†c⟩ from a density matrix |
| `RDM0.m` | Gaussian density matrix from a one-body correlation matrix (Eqs. 25-26) |

Naming in the scripts follows the paper: `Nud`, `Nft`, `Nhd` are the symmetry-resolved negativities $\mathcal{N}^{\uparrow\downarrow}$, $\mathcal{N}^{\mathrm{F}t}$, $\mathcal{N}^{\mathrm{hd}}$ (plotted as $\log_2(2\mathcal{N}+1)$); `Eud`, `Ehd` are the symmetry-resolved relative entropies $E^{\uparrow\downarrow}$, $E^{\mathrm{hd}}$; `logN` is the full fermionic negativity $\mathcal{N}^{\mathrm{F}}$; `MI`/`Iij` is the mutual information $I_{ij}$. All entropies and entanglement measures are in bits.

## Data

### Model and conventions

Hubbard model on the square lattice (Eq. 1), solved with 2×2 plaquette CDMFT and an exact-diagonalization solver at zero temperature, in the paramagnetic phase, with two bath replicas of the cluster (Appendix B). Energies are in units of the half-bandwidth $D = 4t$, with $t = 0.25$, so $D = 1$.

| Dataset | Fixed | Scanned |
|---|---|---|
| `half-filled/` | μ = 0 (half filling) | U/D = 0.1 … 3.0 (step 0.1), 3.5 … 8.0 (step 0.5) |
| `doped/` | U/D = 2.3 | μ/D = −1.00 … −0.36 (step 0.01) |
| `half-filled/CRISBreference/` | μ = 0 (half filling), 100×100 momenta | U/D = 0.1 … 4.0 (step 0.1) |

Each parameter point has its own folder, named `U=<value>` or `xmu=<value>` (e.g. `U=2.000000`, `xmu=-0.500000`). `U_list.txt`/`xmu_list.txt` and `U_conv.txt`/`xmu_conv.txt` list the computed and the converged points.

The point μ/D = −0.36 is not converged: its self-consistency loop ended at an error of about 5·10⁻⁴, above the 10⁻⁵ threshold, so it is missing from `xmu_conv.txt`. It lies in the Mott plateau, its density matrices differ from those at μ/D = −0.37 by about 10⁻³, and it is kept for completeness. The message of the paper is not affected by either including or discarding it.

### Density matrices

The files used by the scripts are the reduced density matrices in each CDMFT run folder:

| File | Content |
|---|---|
| `reduced_density_matrix_1sites.dat` | single-site density matrix (4×4) |
| `reduced_density_matrix_i<X>l1.dat` | single-site density matrix of cluster site X (4×4, `doped/` only) |
| `reduced_density_matrix_2sites.dat` | two-site density matrix of sites 1 and 2 |
| `reduced_density_matrix_i<X>l1_i<Y>l1.dat` | two-site density matrix of cluster sites X and Y |

Plaquette sites are numbered

```
3 4
1 2
```

so (1,2), (3,4), (1,3), (2,4) are nearest neighbors and (1,4), (2,3) are next-nearest neighbors. The scripts use the pair (1,2); `l1` is the orbital index (single band).

Two-site files from CDMFT are mostly 32×16 text matrices: rows 1-16 hold the real part and rows 17-32 the imaginary part of the 16×16 density matrix (read them with `QcmP.post.get_Hloc`). In `doped/`, some of the `i<X>l1_i<Y>l1` files contain only the real 16×16 matrix (read them with `load`); the scripts try `get_Hloc` first and fall back to `load`. The basis is the one of Table 1 in the paper: state $s = 1 + n_{X\uparrow} + 2n_{Y\uparrow} + 4n_{X\downarrow} + 8n_{Y\downarrow}$. Single-site matrices use $s = 1 + n_\uparrow + 2n_\downarrow$, i.e. the order $|\bullet\rangle, |\uparrow\rangle, |\downarrow\rangle, |\uparrow\downarrow\rangle$.

The CRISB folders contain `reduced_density_matrix_i1l1.dat` (4×4) and `reduced_density_matrix_i1l1_i2l1.dat` (16×16, real), in the same bases, together with the CRISB input (`input.grisb`, `lattice.input`) and convergence files.

### Other files in each CDMFT folder

- `inputHM2D.conf`, `used.inputHM2D.conf`: input parameters of the run.
- `eentropy_1sites.dat`, `eentropy_2sites.dat`: von Neumann entropies (bits) of `reduced_density_matrix_1sites.dat` and `reduced_density_matrix_2sites.dat`.
- `LOG.subtraces`: log of the program that traced the cluster density matrix down to the `i<X>l1_i<Y>l1` files.
- Everything else (`sig_*`, `fit_weiss_*`, `observables_*`, `energy_*`, …) is standard output of the ED solver, kept for completeness.

### Summary tables

The `*.txt` files at the top of `half-filled/` and `doped/` (e.g. `dens__site001.txt`, `nup_1_site001.txt`, `eentropy_line_1sites.txt`) collect one value per run folder, sorted by increasing U or μ. They are produced by QcmP-lab when the scripts read the data.

### Helper scripts in the data folders

- `check_subtraces.m`: checks that the four nearest-neighbor and the two next-nearest-neighbor density matrices of the plaquette coincide (Appendix B). Run it from inside `half-filled/` or `doped/`, after `addpath ../lib/qcmp-lab`.
- `subtraces.m`: calls the external tracing program that produced the `i<X>l1_i<Y>l1` files. The program is not included; its output is.
- `half-filled/gather_eigvals.m`: plots the eigenvalues listed by the ED solver (`eigenvalues_list.ed`) as a function of U.

## Library

`lib/qcmp-lab` is a frozen copy of [QcmP-lab](https://github.com/beddalumia/QcmP-LAB), which reads the solver output, together with the `colorlab` color library it bundles (MIT license, see `lib/qcmp-lab/lib/colorlab/LICENSE`).
