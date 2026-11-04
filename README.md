# perturbdyngen

perturbdyngen simulates knockdown (CRISPRi) and activation (CRISPRa) Perturb-seq screens on
[dyngen](https://github.com/dynverse/dyngen) gene regulatory network models.

- **Paired cells.** Every cell is simulated twice, wild type and perturbed. Both copies share
  one random number per reaction and time step (common random numbers). A gene upstream of the
  perturbation therefore has identical counts in both copies, and any difference is caused by
  the perturbation.
- **Perturb-seq observation.** Guide UMIs and RNA UMIs are drawn on top of the true molecule
  counts. The RNA layer captures spliced and unspliced counts with a negative binomial whose
  mean scales with a per-cell library size shared with the guide layer.
- **Benchmark networks.** The package includes the 11 networks and the 22 configuration files
  behind the simulated screens of the CausalODE benchmark, so each dataset can be regenerated
  exactly.

## Installation

```r
# install.packages("remotes")
remotes::install_github("welch-lab/perturbdyngen")
```

Requirements:

- R ≥ 4.1, [dyngen](https://cran.r-project.org/package=dyngen) and
  [GillespieSSA2](https://cran.r-project.org/package=GillespieSSA2).
- A C++ compiler. GillespieSSA2 compiles each network's reactions when a simulation starts.
- To write `.h5ad` files: the R package `reticulate` and a Python with `anndata`. Choose the
  Python with `options(perturbdyngen.python = "/path/to/python")` or the
  `PERTURBDYNGEN_PYTHON` environment variable.

The benchmark data were produced with R 4.4.0, dyngen 1.0.5 and GillespieSSA2 0.3.0.

## Regenerate a benchmark dataset

```r
library(perturbdyngen)

bundled_configs("datasets")                 # the 22 screens
simulate_screen("sp99_knockdown", out_dir = "sp99_knockdown", n_cores = 36)
export_h5ad("sp99_knockdown", "sp99_knockdown/sp99_full_KD0p0050.h5ad", tag = "sp99",
            experiments = grep("_KD0p0050$", rownames(read_run("sp99_knockdown")$design_onehot),
                               value = TRUE))
```

From the shell:

```bash
Rscript $(Rscript -e 'cat(system.file("scripts/perturbdyngen.R", package = "perturbdyngen"))') \
    sp99_knockdown sp99_knockdown --n_cores=36
```

**Speed and size.** A knockdown screen has 380 experiments × 300 cells. On 36 cores it takes
about 4 h for sp99 and about 10 h for sp70, needs about 40 GB of memory and writes about
20 GB. Experiments are saved one by one, and each stage (molecule counts, guide UMIs, RNA UMIs)
is skipped when its file exists, so an interrupted run resumes where it stopped.

**Subsets.** `experiments = c("EXP_1_KD0p0050", ...)` runs only those experiments, and gives the
same files as the full run. Every cell has its own seed and every experiment's seed comes from
its row in the full design, so lowering `n_sims` reproduces the first cells of each experiment.

## Bundled networks

Every network has the same 222 transcription factors. All networks except `legacy` reuse its
differentiation path, simulation settings and per-cell seeds.

| Network | In the paper | Edges | Self-loops | Out-degree CV | Repressive edges | Origin |
|---|---|---|---|---|---|---|
| `legacy` | Branching | 472 | 20 | 1.43 | 24.8% | dyngen branching backbone (the template) |
| `dyngen_trifurcating` | Trifurcating | 353 | 12 | 1.49 | 25.2% | dyngen backbone, hierarchical trifurcation |
| `dyngen_cycle` | Cycle | 383 | 22 | 1.16 | 21.1% | dyngen backbone, linked cycles |
| `collectri_core222` | CollecTRI core | 4,309 | 7 | 0.95 | 11.8% | TF → TF edges of 222 human TFs in CollecTRI |
| `sp99` | sp99 | 473 | 20 | 1.57 | 22.7% | sampled, sparsity 0.99 |
| `sp90_hubLo` / `sp90_hubHi` | sp90 hubLo / hubHi | 4,928 | 20 | 0.92 / 1.90 | 12.8% / 12.9% | sampled, sparsity 0.90 |
| `sp80_hubLo` / `sp80_hubHi` | sp80 hubLo / hubHi | 9,857 | 20 | 0.78 / 1.36 | 6.7% / 6.6% | sampled, sparsity 0.80 |
| `sp70_hubLo` / `sp70_hubHi` | sp70 hubLo / hubHi | 14,785 | 20 | 0.63 / 1.08 | 4.5% / 4.4% | sampled, sparsity 0.70 |

Each directory in `inst/networks/<name>/` holds:

| File | Content |
|---|---|
| `model.rds` | dyngen model, without its stored simulations and gold standard |
| `seeds.rds` | seed of each simulated cell |
| `design_random.rds` | 57 random experiments, each TF perturbed in two, up to 8 TFs per experiment |
| `design_hub.rds` | 19 hub experiments: each perturbs the activating regulators of one high in-degree gene |
| `design_hub_original.rds` | `legacy` only: the original hub design, which mixes regulators of several hubs |
| `grn_truth.rds` | true network, `truth[regulator, target] = effect × strength` |
| `module_network.tsv` | edge list with sign, strength and Hill coefficient |
| `gene_to_node.csv` | external topologies: which gene was placed on which node |

The CollecTRI edges come from Müller-Dott et al., *Nucleic Acids Research* 51, 10934 (2023).
`load_network("sp99")` returns these as a list. `network_dir()` gives the path.

## Bundled dataset configs

The datasets are `<network>_knockdown` and `<network>_activation`, for each of the 11 networks.

| Config | Levels | Experiments | Snapshots |
|---|---|---|---|
| `<network>_knockdown` | transcription × 0.005, 0.2, 0.4, 0.6, 0.8 | (57 + 19) × 5 = 380 | t = 80–400, perturbation at t = 100 |
| `legacy_knockdown` | × 0.005 | 57 + 19 original hub = 76 | t = 100–400 |
| `<network>_activation` | transcription × 2, 5, 10 | (57 + 19) × 3 = 228 | t = 80–400, perturbation at t = 100 |

`legacy_knockdown` reproduces the original Branching knockdown screen: one level, the original hub
design, no snapshots before the perturbation.

A config is an R file that defines `cfg <- list(...)`. Every field can be changed, and any field
left out takes its value from `default_dataset_config()`:

| Field | Meaning |
|---|---|
| `network` | bundled network name, or a directory written by `build_network()` |
| `hub_design` | `"pure"` (`design_hub.rds`) or `"original"` (`design_hub_original.rds`) |
| `perturbation`, `levels` | `"knockdown"` (multipliers in (0, 1)) or `"activation"` (> 1); one run of every experiment per level |
| `experiments` | `"all"`, `"random"` or `"hub"` |
| `n_sims`, `seed` | cells per experiment; run seed (experiment seed = `seed + 1000 × row`) |
| `pre_kd_window` | time recorded before the perturbation, as a fraction of the total time |
| `snap_time` | readout snapshot for the RNA UMIs: `"last"` or a time |
| `guide_sigma_s2` | variance of the log library-size factor |
| `guide_nb_theta`, `umi_capture_rate` | guide UMI negative binomial size and capture probability |
| `sigmoid_k`, `sigmoid_u0` | sigmoid linking guide UMIs to knockdown level |
| `act_guide_kd_equiv` | activation: guide UMIs are simulated as for this knockdown level |
| `p_rna`, `p_premrna`, `rna_nb_theta` | capture probabilities of spliced and unspliced RNA; RNA negative binomial size |

Change fields without editing the file:

```r
run_config("sp99_knockdown", "out", overrides = list(n_sims = 20, levels = 0.005))
```

On the command line, add `--n_sims=20`.

## Model

**Simulation.** Each cell runs dyngen's reaction system with GillespieSSA2 tau-leaping (τ = 300/3600) through the
burn-in and up to t = 100. That is one quarter of the total time of 400. From there two copies continue
with `crn_etl_paired()`:

```
u_r ~ U(0, 1)                            one draw per reaction r and step
n_r(WT)   = qpois(u_r, a_r(x_WT)   · τ)
n_r(pert) = qpois(u_r, a_r(x_pert) · τ)
```

The perturbed copy has the transcription rate of each perturbed gene multiplied by the level. The
gene keeps its own regulation. Snapshots are recorded every 5 time units.

**Guide UMIs** (per cell c and targeted gene g):

```
s_c        ~ LogNormal(0, guide_sigma_s2)
true_cg    ~ NB(mu = s_c · u / E[s], size = guide_nb_theta),   u = u0 + log(1/level − 1) / k
observed   ~ Binomial(true_cg, umi_capture_rate)
```

**RNA UMIs** at the readout snapshot, for the perturbed and the wild-type copy:

```
observed_cg = min(NB(mu = count_cg · p · s_c, size = rna_nb_theta), count_cg)
p = p_rna (spliced) or p_premrna (unspliced)
```

The same `s_c` is used for guides and RNA.

## Outputs of `simulate_screen()`

| Path | Content |
|---|---|
| `checkpoints/<exp>.rds` | true counts: `mat_wt`, `mat_kd` (molecules × snapshots × cells), `snap_times`, `mol_ids`, `perturbed_genes` |
| `perturb_seq/individual/<exp>.rds` | guide UMIs: `true_umis`, `obs_umis`, `kd_continuous`, size factors `s_cells` |
| `perturb_seq/rna_individual/<exp>.rds` | RNA UMIs at the readout: `obs_spliced_kd`, `obs_unspliced_kd`, `obs_spliced_wt`, `obs_unspliced_wt` (cells × genes) |
| `design/` | one-hot and multiplier designs, feature info, true network |
| `model_index.csv`, `run_manifest.csv`, `perturb_seq/*.csv` | experiment index and all settings |

Molecule ids are `mol_premrna_<gene>` (unspliced), `mol_mrna_<gene>` (spliced) and `mol_protein_<gene>`.

## Export to AnnData

`export_h5ad(run_dir, file)` writes one row per cell and snapshot:

| Slot | Content |
|---|---|
| `X`, `layers["spliced"]`, `layers["unspliced"]` | noisy UMIs (CSR, int32) |
| `layers["spliced_true"]`, `layers["unspliced_true"]` | true molecule counts |
| `obs` | experiment, condition (`KD` / `ACT` / `WT`), level, `sim_time`, `readout`, `size_factor`, `kd_continuous`, … |
| `obsm["perturb_design"]`, `obsm["guide_obs_umis"]`, `obsm["guide_true_umis"]` | perturbed genes and guide UMIs |
| `uns["grn_truth"]` | true network |

Away from the readout snapshot, the noisy counts are drawn with their own seeds, as if the same cell were sequenced at every time
point. Rows from a real readout have `obs["readout"] == True`. The wild type is written once by default
(`wt = "single"`). All experiments use the same cell seeds, so the true wild-type counts are the same
in every experiment.

## Build a network

The configs in `inst/configs/networks/` record how each network except `legacy` was built:

```r
build_network("sp99", out_dir = "my_sp99")          # rebuild from its config
cfg <- read_config("sp90_hubHi"); cfg$sdlog <- 1.5; cfg$seed <- 1L
build_network(cfg, "my_network")                    # a new network
simulate_screen(list(network = "my_network", perturbation = "knockdown"), "my_screen")
```

Steps:

1. Sample out-degrees from a log-normal (`sdlog` sets hub dominance), scaled to the edge budget and capped at `max_out`.
2. Draw targets, signs (75% activating), Hill coefficients (2 or 4) and self-loops on the highest out-degree nodes.
   Allow at most 3 repressive inputs per gene.
3. Add the burn-in scaffolding: a root, a 5-node activating spine and 10% of nodes as burn-in.
4. Restore the exact edge budget, and add edges until the root reaches every node.
5. Build the dyngen model with the template's path and settings. Protein half-lives are × 3.
   Rebuild the hub design for the new topology.

`topology = "external"` reads a `regulator`, `target`, `sign` table instead of sampling. The bundled
tables are in `inst/topologies/`. Rebuilding a bundled network from its config gives the same model
under R 4.4 with dyngen 1.0.5.

## Tests

```r
testthat::test_local()
```
