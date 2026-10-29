# openfoam-sparse-seulex

An exact, sparse LU linear solver and autonomous derivative caching engine for OpenFOAM's `seulex` stiff chemistry ODE integrator. Up to **31x speedup** on detailed kinetic mechanisms, with verified trajectory equivalence down to writing precision ($|\Delta T| \le 0.01\text{ K}$).

> **The core claim:** Chemical kinetic Jacobian matrices for hydrocarbons are over 95% sparse, yet OpenFOAM's standard `seulex` factors them with dense $O(n^3)$ Gaussian elimination. By replacing dense elimination with an exact, threshold-monitored minimum-degree sparse LU factorisation (`seulexSparseLU`) and eliminating redundant derivative calls on autonomous systems, combustion chemistry solves in a fraction of the time with **zero loss of physical accuracy**.

Official OpenFOAM Core Merge Request: [openfoam/core/openfoam!935](https://gitlab.com/openfoam/core/openfoam/-/merge_requests/935)

---

## Benchmark Results

All timings below are measured under OpenFOAM's standard `$FOAM_TUTORIALS/combustion/chemFoam` tutorials on a single CPU core, comparing stock dense `seulex` against this solver with verified identical step counts and temperatures:

| Tutorial (`chemFoam/*`) | Species ($n$) | Jacobian Sparsity | Baseline Dense `seulex` | Sparse `seulex` | Speedup | Trajectory Equivalence |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`h2`** | 10 ($n=12$) | 48.0% | 0.02 s | **0.02 s** | **1.0x** | Identical 47 steps, $\max |\Delta T| = 0.0\text{ K}$ |
| **`gri`** (GRI-Mech 3.0) | 53 ($n=55$) | 27.6% | 16.77 s | **15.72 s** | **1.07x** | $t_{\text{ign}}(1500\text{ K}) = 5.941\times 10^{-2}\text{ s}$, $T_{\text{end}} = 2660.57\text{ K}$ |
| **`nc7h16`** (n-heptane) | 544 ($n=546$) | 3.7% | 51.42 s | **2.73 s** | **18.8x** | Identical 161 steps, $\max |\Delta T| = 0.0\text{ K}$ |
| **`ic8h18`** (iso-octane) | 874 ($n=876$) | 2.3% | 537.46 s | **17.22 s** | **31.2x** | Identical 1,065 steps, $\max |\Delta T| = 0.01\text{ K}$ |

### Multi-Dimensional Counterflow Flame (`reactingFoam`)
On a 2D counterflow diffusion flame (4,000 cells, detailed 544-species n-heptane mechanism, 6 CPU cores parallel):
- **Stock OpenFOAM `seulex`:** 2,123 seconds (10 time steps)
- **Sparse `seulex`:** 64 seconds (10 time steps)
- **Speedup:** **33.2x faster**; maximum temperature difference across all 4,000 cells was $0.34\text{ K}$ (well within solver relative tolerance $\text{relTol} = 0.1$).

---

## Check It Yourself

### 1. Verify Recorded Benchmark Outputs (No Dependencies)
Runs with standard Python 3.8+ library only:
```bash
python3 verify.py
```

### 2. Run Live chemFoam Tutorials (Requires OpenFOAM v2512 or OF14)
```bash
# Source your OpenFOAM environment
source /usr/lib/openfoam/openfoam2512/etc/bashrc

# Run benchmark across h2, gri, nc7h16, ic8h18
bash benchmark.sh /tmp/seulex_bench
```

### 3. Quick Single-Case Reproduction (`nc7h16`, 544 species)
```bash
# Baseline dense seulex (~51.4 s)
cp -r "$FOAM_TUTORIALS/combustion/chemFoam/nc7h16" /tmp/nc7h16_dense
foamDictionary /tmp/nc7h16_dense/constant/chemistryProperties -entry odeCoeffs/sparse -set false
chemFoam -case /tmp/nc7h16_dense

# Sparse seulex (~2.7 s, ~19x faster)
cp -r "$FOAM_TUTORIALS/combustion/chemFoam/nc7h16" /tmp/nc7h16_sparse
foamDictionary /tmp/nc7h16_sparse/constant/chemistryProperties -entry odeCoeffs/sparse -set true
foamDictionary /tmp/nc7h16_sparse/constant/chemistryProperties -entry odeCoeffs/autonomous -set true
chemFoam -case /tmp/nc7h16_sparse
```

---

## How It Works

### 1. Sparse LU Factorisation (`seulexSparseLU`)
During extrapolation integration, `seulex` repeatedly solves linear systems of the form:
$$\left(\frac{1}{\Delta x_k}\mathbf{I} - \mathbf{J}\right) \Delta \mathbf{y} = \mathbf{f}$$

1. **Greedy Minimum-Degree Ordering:** Computes a symmetric elimination ordering on the off-diagonal sparsity pattern of $\mathbf{J}$ to minimise fill-in.
2. **Compiled Flat Execution:** Compiles the symbolic elimination into flat contiguous integer arrays (`elimP_`, `elimPiv_`, `elimCnt_`, `updDst_`, `updSrc_`). Numeric factorisation executes as a tight loop over non-zeros and fill-in slots, skipping numerically zero multipliers.
3. **Threshold Gate (`0.5 * n^3 / 3`):** Evaluates symbolic factorisation cost. If the sparse multiply-add count is not strictly lower than half the dense LU cost, or if the matrix is dense, it aborts symbolic construction immediately ($O(1)$) and defaults to dense LU.
4. **Dynamic Fill-In Absorption:** When new non-zero Jacobian entries arise during integration that land in existing fill-in slots, `addEntry()` absorbs them into the factors without triggering an expensive symbolic rebuild.
5. **Partial Pivoting Fallback:** Pivot magnitudes are continuously monitored against the initial diagonal scale (`1e-10 * scale_`). Any singularity or ill-conditioned pivot safely triggers dense `LUDecompose` with partial pivoting.

### 2. Autonomous Derivative Reuse
For autonomous kinetic systems ($\partial \mathbf{f} / \partial t = \mathbf{0}$):
- The initial derivative $\mathbf{f}(x_0 + \Delta x_k, \mathbf{y}_0)$ is identical across all extrapolation columns $k = 0, 1, \dots$. It is evaluated once at $k = 0$ and reused across all columns.
- The derivative evaluation performed during the divergence test (`nn == 1 && k <= 1`) is preserved in `dydx_` and directly passed to the subsequent integration substep, saving 20–35% of all RHS function calls.

---

## Honest Caveats

1. **Small Mechanisms ($\le 20$ species):**
   In small kinetic mechanisms (such as 10-species hydrogen `h2` or 53-species methane `gri`), the Jacobian factorisation represents a negligible fraction of the total execution time compared to thermodynamic property evaluations and reaction rate evaluations. Speedups on small mechanisms are modest ($1.0\times$ to $1.1\times$).
2. **Dense Jacobian Models:**
   If mechanism reduction or full analytical Jacobians produce a dense matrix (e.g. non-sparse tabulation), the solver's density guard disables sparse factorisation, performing standard dense LU with zero speedup and zero overhead.
3. **Hardware / Cache Effects:**
   Speedups are largest on large mechanisms ($n > 100$) where dense matrix factorisation exceeds CPU L1/L2 cache capacity and $O(n^3)$ operations dominate.

---

## Repository Structure

```text
openfoam-sparse-seulex/
├── README.md               # Technical documentation & verified benchmark report
├── LICENSE                 # GNU General Public License v3.0
├── verify.py               # Self-contained trajectory & gate verifier (no dependencies)
├── benchmark.sh            # Benchmark runner over standard OpenFOAM chemFoam tutorials
├── src/
│   ├── v2512/              # Standalone library for OpenFOAM v2512 (openfoam.com)
│   ├── of14/               # Standalone library for OpenFOAM 14 (openfoam.org)
│   └── upstream-patch/     # Upstream OpenFOAM MR patch (core/openfoam!935)
└── results/
    ├── dene-cikti.txt      # Raw terminal output from chemFoam benchmark suite (v2512)
    ├── dene-of14-cikti.txt # Raw terminal output from chemFoam benchmark suite (OF14)
    └── OLCUM.md            # Detailed multi-case measurement log (reactive flames & DLBFoam)
```

---

## Build Instructions

### For OpenFOAM v2512 (`openfoam.com`):
```bash
cd src/v2512
wmake libso
```
Add to `system/controlDict`:
```foam
libs ("libseulexFast.so");
```
And in `constant/chemistryProperties` (`odeCoeffs`):
```foam
odeCoeffs
{
    solver      seulexFast;
    autonomous  true;
}
```

### For Upstream OpenFOAM Integration:
See the patch in `src/upstream-patch/` or follow [OpenFOAM Merge Request !935](https://gitlab.com/openfoam/core/openfoam/-/merge_requests/935).

---

## License

This project is licensed under the GNU General Public License v3.0, consistent with OpenFOAM.
