# Three-Member Community Flux Sampling

Python workflow for flux sampling of an already assembled three-member
microbial community model in SBML format. The script uses the
Artificially Centered Hit-and-Run (ACHR) sampler implemented in
COBRApy.

The workflow was adapted from the joint-community sampling logic in
[FinleyLabUSC/FluxSampling_MetabolicInteractions](https://github.com/FinleyLabUSC/FluxSampling_MetabolicInteractions).
It samples one pooled community model; it does not reconstruct the
community from separate organism models.

## Contents

-   `sample_three_member_community.py` --- sampling script.
-   `mcaps_oryzae_agri.xml` --- SBML model.

## Input model

The script expects a single, already assembled SBML community model. The model must contain:

-   one community biomass reaction;
-   one biomass reaction for each of the three community members;
-   finite lower and upper bounds for all reactions.

The current script is configured for the following biomass reactions:

  Community member   Expected SBML reaction ID
  ------------------ ------------------------------------
  A. oryzae             `R_oryzae_bio1`
  M .capsulatus         `R_capsulatus_BIOMASS_Mcapsulatus`
  B. agri               `R_agri_bio1`
  Community          `R_community_biomass`

COBRApy may remove the `R_` prefix when importing SBML. The script
therefore checks both the original ID and the ID without that prefix.

**Important:** these reaction IDs are model-specific. If you use another
community model, update `EXPECTED_MEMBER_BIOMASS` and
`COMMUNITY_BIOMASS` near the beginning of the script.

## Run

Place the script and SBML model in the same directory, or provide the
full path to the model.

Run with all main parameters specified:

``` bash
python -u sample_three_member_community.py \
    --model mcaps_oryzae_agri.xml \
    --outdir sampling_results \
    --n-samples 5000 \
    --thinning 100 \
    --batch-size 250
```

`-u` disables Python output buffering, which is useful when monitoring a
long run or redirecting the output to a log file.

### Command-line arguments

  --------------------------------------------------------------------------
  Argument                                     Default Description
  ----------------------- ---------------------------- ---------------------
  `--model`                    `mcaps_oryzae_agri.xml` Input SBML/XML model
                                                       path

  `--outdir`                        `sampling_results` Directory for output
                                                       files

  `--n-samples`                                 `5000` Total number of
                                                       samples to generate

  `--thinning`                                   `100` ACHR steps between
                                                       retained samples

  `--batch-size`                                 `250` Number of samples
                                                       written to each
                                                       compressed CSV file
  --------------------------------------------------------------------------

## What the script does

1.  Loads the SBML community model with COBRApy.
2.  Resolves the community and member biomass reaction IDs.
3.  Checks that all reaction bounds are finite.
4.  Optimizes the model with FBA and reports the community and member
    biomass fluxes.
5.  Initializes one COBRApy `ACHRSampler` and reuses it for all batches.
6.  Samples the joint feasible flux space subject to the model's
    stoichiometric constraints and reaction bounds.
7.  Writes each batch to a compressed CSV file and calculates
    reaction-level summary statistics across all generated samples.

The model is sampled as one coupled community model. The script does not
sample each organism independently.

## Output files

The output directory contains:

  -----------------------------------------------------------------------
  File                                Description
  ----------------------------------- -----------------------------------
  `flux_samples_0001.csv.gz`,         Sampled flux vectors. Each row is
  `flux_samples_0002.csv.gz`, ...     one sample and each column is a
                                      reaction.

  `flux_summary.csv`                  Reaction-level mean, population
                                      standard deviation, minimum and
                                      maximum across all samples.

  `reaction_order.csv`                Reaction IDs and their column
                                      positions in the sampled batch
                                      files.

  `run_report.json`                   Model information, sampling
                                      parameters, biomass reaction IDs
                                      and fluxes from the FBA validation,
                                      solver information, and elapsed
                                      time.
  -----------------------------------------------------------------------

The number of batch files is determined by `n-samples` and `batch-size`.
For 5,000 samples with a batch size of 250, the script writes 20 batch
files.

Example output structure:

``` text
sampling_results/
├── flux_samples_0001.csv.gz
├── flux_samples_0002.csv.gz
├── ...
├── flux_samples_0020.csv.gz
├── flux_summary.csv
├── reaction_order.csv
└── run_report.json
```

## Reference

[FluxSampling_MetabolicInteractions](https://github.com/FinleyLabUSC/FluxSampling_MetabolicInteractions).
