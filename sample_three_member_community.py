#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, List

import cobra
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from cobra.sampling import ACHRSampler


EXPECTED_MEMBER_BIOMASS = {
    "oryzae": "R_oryzae_bio1",
    "capsulatus": "R_capsulatus_BIOMASS_Mcapsulatus",
    "agri": "R_agri_bio1",
}
COMMUNITY_BIOMASS = "R_community_biomass"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ACHR sampling of a joint 3-member community model")
    p.add_argument("--model", default="mcaps_oryzae_agri.xml", help="Input SBML/XML model")
    p.add_argument("--outdir", default="sampling_results", help="Output directory")
    p.add_argument("--n-samples", type=int, default=5000, help="Total number of samples (repository pairwise workflow uses 5000)")
    p.add_argument("--thinning", type=int, default=100,
                   help="ACHR steps between retained samples (100 is a starting point, not equivalent to RHMC)")
    p.add_argument("--batch-size", type=int, default=250,
                   help="Samples per output file; limits RAM use")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--solver", default=None, help="Optional COBRApy solver, e.g. glpk or cplex")
    p.add_argument("--min-growth", type=float, default=0.001,
                   help="Minimum lower bound for each organism biomass reaction; use 0 to disable")
    p.add_argument("--skip-member-check", action="store_true",
                   help="Do not require all three member biomass fluxes to be positive in FBA")
    return p.parse_args()


def resolve_reaction_id(model: cobra.Model, sbml_reaction_id: str) -> str:

    candidates = [sbml_reaction_id]
    if sbml_reaction_id.startswith("R_"):
        candidates.append(sbml_reaction_id[2:])

    for candidate in candidates:
        if candidate in model.reactions:
            return candidate

    raise KeyError(
        f"Reaction from SBML '{sbml_reaction_id}' was not found after import. "
        f"Tried IDs: {candidates}"
    )


def finite_bounds_check(model: cobra.Model) -> None:
    bad = [r.id for r in model.reactions
           if not math.isfinite(float(r.lower_bound)) or not math.isfinite(float(r.upper_bound))]
    if bad:
        raise ValueError(
            f"Sampling requires finite reaction bounds. {len(bad)} reactions have infinite bounds; "
            f"examples: {bad[:20]}"
        )


def validate_model(model: cobra.Model, skip_member_check: bool, min_growth: float) -> Dict[str, float]:
    global COMMUNITY_BIOMASS, EXPECTED_MEMBER_BIOMASS

    COMMUNITY_BIOMASS = resolve_reaction_id(model, COMMUNITY_BIOMASS)
    EXPECTED_MEMBER_BIOMASS = {
        member: resolve_reaction_id(model, rid)
        for member, rid in EXPECTED_MEMBER_BIOMASS.items()
    }

    print(f"Model ID: {model.id}")
    print(f"Community biomass reaction: {COMMUNITY_BIOMASS}")
    print(f"Member biomass reactions: {EXPECTED_MEMBER_BIOMASS}")
    print(f"Reactions: {len(model.reactions):,}")
    print(f"Metabolites: {len(model.metabolites):,}")
    print(f"Genes: {len(model.genes):,}")
    print(f"Objective: {model.objective.expression}")

    community_rxn = model.reactions.get_by_id(COMMUNITY_BIOMASS)
    if not (math.isclose(community_rxn.lower_bound, 1.0, abs_tol=1e-9)
            and math.isclose(community_rxn.upper_bound, 1.0, abs_tol=1e-9)):
        print("WARNING: community biomass bounds are not fixed at 1:",
              community_rxn.bounds, file=sys.stderr)

    if min_growth < 0:
        raise ValueError("min-growth cannot be negative")
    for rid in EXPECTED_MEMBER_BIOMASS.values():
        rxn = model.reactions.get_by_id(rid)
        rxn.lower_bound = max(float(rxn.lower_bound), float(min_growth))

    finite_bounds_check(model)

    solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError(f"The original model is not feasible/optimal (status={solution.status})")

    print(f"\nFBA status: {solution.status}")
    print(f"Community biomass flux: {solution.fluxes[COMMUNITY_BIOMASS]:.8g}")
    member_growth = {}
    for member, rid in EXPECTED_MEMBER_BIOMASS.items():
        value = float(solution.fluxes[rid])
        member_growth[member] = value
        print(f"{member} biomass flux ({rid}): {value:.8g}")
        if not skip_member_check and value <= 1e-9:
            raise RuntimeError(
                f"Member '{member}' has non-positive biomass flux in the FBA solution. "
                "Check the community biomass reaction and model constraints before sampling."
            )

    return member_growth


def main() -> None:
    args = parse_args()
    if args.n_samples < 1 or args.batch_size < 1 or args.thinning < 1:
        raise ValueError("n-samples, batch-size, and thinning must all be positive integers")

    model_path = Path(args.model).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    if not model_path.is_file():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading SBML: {model_path}")
    model = read_sbml_model(str(model_path))
    if args.solver:
        model.solver = args.solver
    member_growth = validate_model(model, args.skip_member_check, args.min_growth)

    print("\nInitializing ACHR sampler...")
    sampler = ACHRSampler(model, thinning=args.thinning, seed=args.seed)

    n_rxns = len(model.reactions)
    reaction_ids = [r.id for r in model.reactions]
    sum_v = np.zeros(n_rxns, dtype=np.float64)
    sum_v2 = np.zeros(n_rxns, dtype=np.float64)
    min_v = np.full(n_rxns, np.inf, dtype=np.float64)
    max_v = np.full(n_rxns, -np.inf, dtype=np.float64)
    n_written = 0
    batch_index = 0
    start = time.time()

    while n_written < args.n_samples:
        n_now = min(args.batch_size, args.n_samples - n_written)
        print(f"Sampling batch {batch_index + 1}: {n_now} samples "
              f"({n_written}/{args.n_samples} completed)")
        samples = sampler.sample(n_now)
        samples = samples.reindex(columns=reaction_ids)
        values = samples.to_numpy(dtype=np.float64, copy=False)
        if values.shape != (n_now, n_rxns):
            raise RuntimeError(f"Unexpected sample matrix shape: {values.shape}")
        if not np.isfinite(values).all():
            raise RuntimeError("Non-finite values detected in sampled flux matrix")

        sum_v += values.sum(axis=0)
        sum_v2 += np.square(values).sum(axis=0)
        min_v = np.minimum(min_v, values.min(axis=0))
        max_v = np.maximum(max_v, values.max(axis=0))

        batch_index += 1
        batch_file = outdir / f"flux_samples_{batch_index:04d}.csv.gz"
        samples.to_csv(batch_file, index=False, compression="gzip", float_format="%.10g")
        n_written += n_now
        print(f"  saved: {batch_file.name}")

    mean_v = sum_v / n_written
    var_v = np.maximum(sum_v2 / n_written - np.square(mean_v), 0.0)
    std_v = np.sqrt(var_v)
    summary = pd.DataFrame({
        "reaction": reaction_ids,
        "mean": mean_v,
        "std": std_v,
        "min": min_v,
        "max": max_v,
    })
    summary.to_csv(outdir / "flux_summary.csv", index=False, float_format="%.10g")

    pd.DataFrame({"column_index": np.arange(n_rxns), "reaction": reaction_ids}).to_csv(
        outdir / "reaction_order.csv", index=False
    )

    report = {
        "input_model": str(model_path),
        "model_id": model.id,
        "n_reactions": n_rxns,
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
        "sampling_method": "ACHR (COBRApy)",
        "n_samples": n_written,
        "thinning": args.thinning,
        "batch_size": args.batch_size,
        "random_seed": args.seed,
        "solver": str(model.solver.interface),
        "objective_reaction": COMMUNITY_BIOMASS,
        "community_biomass_bounds": list(model.reactions.get_by_id(COMMUNITY_BIOMASS).bounds),
        "member_biomass_reactions": EXPECTED_MEMBER_BIOMASS,
        "minimum_member_growth_bound": args.min_growth,
        "member_biomass_fluxes_in_fba": member_growth,
        "elapsed_seconds": round(time.time() - start, 2),
        "outputs": [
            "flux_samples_*.csv.gz: rows are samples, columns are reactions",
            "flux_summary.csv: mean, population standard deviation, min, max",
            "reaction_order.csv: reaction column order",
            "run_report.json: settings and validation information",
        ],
    }
    with open(outdir / "run_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=2)

    print("\nSampling completed.")
    print(f"Samples: {n_written:,}")
    print(f"Output directory: {outdir}")
    print(f"Elapsed time: {report['elapsed_seconds']} s")
    print("NOTE: ACHR is COBRApy's sampler; it is not a bit-for-bit implementation of MATLAB RHMC.")


if __name__ == "__main__":
    main()
