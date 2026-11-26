"""Finite-copy intracellular reactions and a deliberately reduced cell lineage.

Time is hours, volume is fL, and production rates are nM/h in the named
compartment. Piecewise-constant expression and exponential cell growth are
handled by exact stochastic thinning, not a fixed simulation time step.
"""

from dataclasses import asdict, dataclass
import math
import numbers

import numpy as np

AVOGADRO = 6.02214076e23
MOLECULES_PER_NM_FL = AVOGADRO * 1e-24
POOLS = ("cytoplasmic_Cas9", "nuclear_Cas9", "gRNA1", "gRNA2", "RNP1", "RNP2")
DNA_STATES = (
    "intact",
    "PAM_bound",
    "R_loop",
    "cleaved_Cas9_bound",
    "DSB_released",
    "synthetic_disruptive_allele_A",
    "synthetic_disruptive_allele_B",
    "synthetic_other_resistant_allele",
)


def _number(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be real")
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"Invalid {name}")


def _integer(value, name, low, high):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    if not low <= value <= high:
        raise ValueError(f"{name} outside [{low}, {high}]")


@dataclass(frozen=True)
class Locus:
    name: str
    guide: int
    binding_per_nM_h: float = 6.0
    unbinding_per_h: float = 108.0
    loop_per_h: float = 216.0
    unwind_per_h: float = 54.0
    cleavage_per_h: float = 90.0

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("A locus needs a name")
        _integer(self.guide, "guide", 0, 1)
        for name in (
            "binding_per_nM_h",
            "unbinding_per_h",
            "loop_per_h",
            "unwind_per_h",
            "cleavage_per_h",
        ):
            _number(getattr(self, name), name)


@dataclass(frozen=True)
class CellParameters:
    cas9_production_nM_h: float = 0.12
    guide_production_nM_h: tuple = (1.0, 1.0)
    expression_on_h: float = 0.0
    expression_off_h: float = 12.0
    import_per_h: float = 0.5
    export_per_h: float = 0.05
    protein_decay_per_h: float = 0.1
    guide_decay_per_h: float = 0.2
    rnp_decay_per_h: float = 0.1
    assembly_per_nM_h: float = 0.5
    disassembly_per_h: float = 0.2
    product_release_per_h: float = 1.0
    precise_repair_per_h: float = 0.0791
    mutation_repair_per_h: float = 0.0335
    disruptive_given_mutation: float = 0.7
    nucleus_birth_fL: float = 30.0
    cytoplasm_birth_fL: float = 150.0
    cycle_h: float = 8.0
    horizon_h: float = 24.0
    max_generations: int = 3

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name not in ("guide_production_nM_h", "max_generations"):
                _number(
                    value,
                    name,
                    positive=name
                    in (
                        "nucleus_birth_fL",
                        "cytoplasm_birth_fL",
                        "cycle_h",
                        "horizon_h",
                    ),
                )
        if (
            not isinstance(self.guide_production_nM_h, tuple)
            or len(self.guide_production_nM_h) != 2
        ):
            raise ValueError("Two guide production rates are required")
        for value in self.guide_production_nM_h:
            _number(value, "guide production")
        _integer(self.max_generations, "max_generations", 0, 6)
        if (
            self.expression_off_h < self.expression_on_h
            or self.disruptive_given_mutation > 1
        ):
            raise ValueError("Invalid expression window or repair fraction")
        if self.horizon_h > 96 or self.cycle_h < 0.25:
            raise ValueError("Outside the bounded prototype's time range")


DEFAULT_LOCI = (
    Locus("ABCG18_case_study", 0),
    Locus("ABCG19_case_study", 1),
    Locus("illustrative_off_target", 0, loop_per_h=0.0216),
)


def volume_scale(time, born, parameters):
    """Grow to twice birth volume; growth pauses during the simplified checkpoint."""
    return math.exp(
        math.log(2) * min(max(time - born, 0), parameters.cycle_h) / parameters.cycle_h
    )


def molecule_totals(pools, dna, loci):
    """Cas9 and each guide, including target-bound complexes."""
    totals = np.array(
        [
            pools[0] + pools[1] + pools[4] + pools[5],
            pools[2] + pools[4],
            pools[3] + pools[5],
        ],
        dtype=np.int64,
    )
    for row, locus in zip(dna, loci):
        bound = int(np.count_nonzero((row >= 1) & (row <= 3)))
        totals[0] += bound
        totals[locus.guide + 1] += bound
    return totals


def reactions(pools, dna, loci, parameters, nuclear_fL, cytoplasmic_fL, expression):
    """(hazard, pool delta, DNA update, birth ledger, loss ledger).

    Expression supplies protein in cytoplasm and guide RNA in nucleus. Nuclear
    guide production lumps transcription/processing. RNP assembly is nuclear.
    """
    p = parameters
    n = MOLECULES_PER_NM_FL * nuclear_fL
    c = MOLECULES_PER_NM_FL * cytoplasmic_fL
    result = []

    def add(rate, changes=(), dna_update=None, birth=None, loss=None):
        if rate < 0 or not math.isfinite(rate):
            raise ArithmeticError("Invalid reaction hazard")
        if rate > 0:
            result.append((rate, changes, dna_update, birth, loss))

    if expression:
        add(p.cas9_production_nM_h * c, ((0, 1),), birth=0)
        for guide in range(2):
            add(p.guide_production_nM_h[guide] * n, ((2 + guide, 1),), birth=guide + 1)
    add(p.import_per_h * pools[0], ((0, -1), (1, 1)))
    add(p.export_per_h * pools[1], ((1, -1), (0, 1)))
    for index in (0, 1):
        add(p.protein_decay_per_h * pools[index], ((index, -1),), loss=(0,))
    for guide in range(2):
        g, r = 2 + guide, 4 + guide
        add(p.guide_decay_per_h * pools[g], ((g, -1),), loss=(guide + 1,))
        add(p.rnp_decay_per_h * pools[r], ((r, -1),), loss=(0, guide + 1))
        add(
            p.assembly_per_nM_h * float(pools[1]) * pools[g] / n,
            ((1, -1), (g, -1), (r, 1)),
        )
        add(p.disassembly_per_h * pools[r], ((1, 1), (g, 1), (r, -1)))
    for li, locus in enumerate(loci):
        r = 4 + locus.guide
        for ai, state in enumerate(dna[li]):

            def update(new):
                return li, ai, new

            if state == 0:
                add(locus.binding_per_nM_h * pools[r] / n, ((r, -1),), update(1))
            elif state == 1:
                add(locus.unbinding_per_h, ((r, 1),), update(0))
                add(locus.loop_per_h, dna_update=update(2))
            elif state == 2:
                add(locus.unwind_per_h, dna_update=update(1))
                add(locus.cleavage_per_h, dna_update=update(3))
            elif state == 3:
                add(p.product_release_per_h, ((r, 1),), update(4))
            elif state == 4:
                add(p.precise_repair_per_h, dna_update=update(0))
                add(
                    p.mutation_repair_per_h * p.disruptive_given_mutation / 2,
                    dna_update=update(5),
                )
                add(
                    p.mutation_repair_per_h * p.disruptive_given_mutation / 2,
                    dna_update=update(6),
                )
                add(
                    p.mutation_repair_per_h * (1 - p.disruptive_given_mutation),
                    dna_update=update(7),
                )
            if state in (1, 2, 3):
                add(
                    p.rnp_decay_per_h,
                    dna_update=update(4 if state == 3 else 0),
                    loss=(0, locus.guide + 1),
                )
    return result


def partition(pools, dna, loci, rng):
    """Resolve reversible binding, replicate DNA states, partition free molecules.

    Division is blocked while a DSB or cleaved, Cas9-bound target remains.
    Reversible pre-cleavage complexes dissociate before replication. These are
    explicit prototype assumptions, not calibrated cell-cycle biology.
    """
    if np.any((dna == 3) | (dna == 4)):
        raise ValueError("Unrepaired breaks block division in this prototype")
    free, copied = pools.copy(), dna.copy()
    for li, locus in enumerate(loci):
        bound = (copied[li] == 1) | (copied[li] == 2)
        free[4 + locus.guide] += int(bound.sum())
        copied[li, bound] = 0
    left = rng.binomial(free, 0.5)
    return (left, copied.copy()), (free - left, copied.copy())


def _genotype(dna):
    # Intact binding states share genotype 0; -1 denotes an unresolved break.
    codes = np.where(dna < 3, 0, np.where(dna < 5, -1, dna))
    return tuple(tuple(sorted(map(int, row))) for row in codes)


def simulate_founder(
    parameters=CellParameters(),
    *,
    loci=DEFAULT_LOCI,
    times_h=None,
    seed=7,
    event_budget=200000,
    initial_pools=None,
    initial_dna=None,
):
    """One independently sampled founder and all modeled descendants."""
    if not isinstance(parameters, CellParameters):
        raise ValueError("Expected validated cell parameters")
    loci = tuple(loci)
    if any(not isinstance(x, Locus) for x in loci) or len(
        {x.name for x in loci}
    ) != len(loci):
        raise ValueError("Loci must be validated and uniquely named")
    _integer(seed, "seed", 0, 2**32 - 1)
    _integer(event_budget, "event_budget", 1, 1000000)
    times = np.asarray([parameters.horizon_h] if times_h is None else times_h)
    if times.dtype.kind not in "iuf" or times.ndim != 1 or not 1 <= len(times) <= 1000:
        raise ValueError("Invalid observation grid")
    if (
        not np.isfinite(times).all()
        or np.any(times < 0)
        or np.any(np.diff(times) <= 0)
        or times[-1] > parameters.horizon_h
    ):
        raise ValueError("Observation times must increase within the experiment")
    pools = (
        np.zeros(6, dtype=np.int64)
        if initial_pools is None
        else np.asarray(initial_pools)
    )
    dna = (
        np.zeros((len(loci), 2), dtype=np.int64)
        if initial_dna is None
        else np.asarray(initial_dna)
    )
    if (
        pools.shape != (6,)
        or pools.dtype.kind not in "iu"
        or np.any(pools < 0)
        or np.any(pools > 50000)
    ):
        raise ValueError("Invalid initial molecule counts")
    if (
        dna.shape != (len(loci), 2)
        or dna.dtype.kind not in "iu"
        or np.any(dna > 7)
        or np.any(dna < 0)
    ):
        raise ValueError("Invalid initial diploid DNA states")
    initial_total = molecule_totals(pools, dna, loci)
    rng = np.random.default_rng(seed)
    births, losses = np.zeros(3, dtype=np.int64), np.zeros(3, dtype=np.int64)
    # Columns: cells, double-biallelic proxy cells, off-target-mutant cells,
    # cells with DSBs, RNP1 count, RNP2 count, total nuclear volume in fL.
    observed = np.zeros((len(times), 9))
    terminals, counters = (
        [],
        {
            "accepted_reactions": 0,
            "thinning_proposals": 0,
            "divisions": 0,
            "cut_events": 0,
            "precise_repair_events": 0,
            "mutation_repair_events": 0,
        },
    )
    pending = [(0.0, 0, pools.copy(), dna.copy())]
    while pending:
        born, generation, pools, dna = pending.pop()
        now, index = born, int(np.searchsorted(times, born, side="left"))
        due = born + parameters.cycle_h

        def observe(stop, inclusive=False):
            nonlocal index
            while index < len(times) and (
                times[index] <= stop if inclusive else times[index] < stop
            ):
                scale = volume_scale(times[index], born, parameters)
                double = len(loci) >= 2 and bool(
                    np.all((dna[:2] == 5) | (dna[:2] == 6))
                )
                off = len(loci) >= 3 and bool(np.any(dna[2:] >= 5))
                mutated = [
                    int(np.count_nonzero(dna[li] >= 5)) if len(loci) > li else 0
                    for li in range(2)
                ]
                observed[index] += (
                    1,
                    double,
                    off,
                    np.any((dna == 3) | (dna == 4)),
                    pools[4],
                    pools[5],
                    parameters.nucleus_birth_fL * scale,
                    mutated[0],
                    mutated[1],
                )
                index += 1

        while True:
            broken = bool(np.any((dna == 3) | (dna == 4)))
            if now >= due and generation < parameters.max_generations and not broken:
                observe(now)
                before = molecule_totals(pools, dna, loci)
                children = partition(pools, dna, loci, rng)
                if not np.array_equal(
                    sum((molecule_totals(a, b, loci) for a, b in children)), before
                ):
                    raise ArithmeticError("Division created or lost molecular material")
                for child_pools, child_dna in children:
                    pending.append((now, generation + 1, child_pools, child_dna))
                counters["divisions"] += 1
                break
            if now >= parameters.horizon_h:
                observe(parameters.horizon_h, inclusive=True)
                terminals.append(
                    {
                        "pools": pools.tolist(),
                        "dna": dna.tolist(),
                        "genotype": [list(row) for row in _genotype(dna)],
                        "generation": generation,
                    }
                )
                break
            boundaries = [parameters.horizon_h]
            if now < due:
                boundaries.append(due)
            boundaries.extend(
                t
                for t in (parameters.expression_on_h, parameters.expression_off_h)
                if t > now
            )
            boundary = min(boundaries)
            expression = parameters.expression_on_h <= now < parameters.expression_off_h
            v0 = volume_scale(now, born, parameters)
            v1 = volume_scale(boundary, born, parameters)
            # Birth hazards increase with volume; bimolecular hazards decrease.
            # Their endpoint maxima bound every propensity on this segment.
            current = reactions(
                pools,
                dna,
                loci,
                parameters,
                parameters.nucleus_birth_fL * v0,
                parameters.cytoplasm_birth_fL * v0,
                expression,
            )
            upper = sum(x[0] for x in current)
            if expression:
                production = (
                    parameters.cas9_production_nM_h * parameters.cytoplasm_birth_fL
                    + sum(parameters.guide_production_nM_h)
                    * parameters.nucleus_birth_fL
                ) * MOLECULES_PER_NM_FL
                upper += production * (v1 - v0)
            if upper == 0:
                observe(boundary)
                now = boundary
                continue
            counters["thinning_proposals"] += 1
            if counters["thinning_proposals"] > event_budget:
                raise RuntimeError(
                    "Cell reaction budget exhausted; no partial experiment returned"
                )
            event_time = now + float(rng.exponential(1 / upper))
            if event_time >= boundary:
                observe(boundary)
                now = boundary
                continue
            observe(event_time)
            now = event_time
            scale = volume_scale(now, born, parameters)
            available = reactions(
                pools,
                dna,
                loci,
                parameters,
                parameters.nucleus_birth_fL * scale,
                parameters.cytoplasm_birth_fL * scale,
                expression,
            )
            total = sum(x[0] for x in available)
            if total > upper * (1 + 1e-12):
                raise ArithmeticError("Stochastic thinning envelope was violated")
            draw = rng.random() * upper
            if draw >= total:
                continue
            cumulative = 0.0
            for hazard, changes, update, birth, loss in available:
                cumulative += hazard
                if draw < cumulative:
                    for pool, delta in changes:
                        pools[pool] += delta
                    if update is not None:
                        previous = dna[update[0], update[1]]
                        counters["cut_events"] += int(previous == 2 and update[2] == 3)
                        counters["precise_repair_events"] += int(
                            previous == 4 and update[2] == 0
                        )
                        counters["mutation_repair_events"] += int(
                            previous == 4 and update[2] >= 5
                        )
                        dna[update[0], update[1]] = update[2]
                    if birth is not None:
                        births[birth] += 1
                    if loss is not None:
                        for component in loss:
                            losses[component] += 1
                    counters["accepted_reactions"] += 1
                    break
            if np.any(pools < 0) or np.any(pools > 50000):
                raise ArithmeticError(
                    "Molecule counts left the prototype's permitted range"
                )
    final_total = sum(
        (
            molecule_totals(np.asarray(t["pools"]), np.asarray(t["dna"]), loci)
            for t in terminals
        )
    )
    ledger_error = final_total - (initial_total + births - losses)
    if np.any(ledger_error) or np.any(observed[:, 0] < 1):
        raise ArithmeticError("Molecular conservation or observation coverage failed")
    genotypes = {_genotype(np.asarray(t["dna"])) for t in terminals}
    resolved = [
        t
        for t in terminals
        if not np.any((np.asarray(t["dna"]) == 3) | (np.asarray(t["dna"]) == 4))
    ]
    resolved_genotypes = {_genotype(np.asarray(t["dna"])) for t in resolved}
    return {
        "times_h": times.tolist(),
        "observations": observed.tolist(),
        "terminal_cells": terminals,
        "terminal_cells_count": len(terminals),
        "mixed_terminal_states": len(genotypes) > 1,
        "resolved_terminal_cells_count": len(resolved),
        "mosaicism_evaluable": len(resolved) >= 2,
        "mosaicism_proxy": len(resolved_genotypes) > 1,
        "distinct_resolved_synthetic_genotypes": len(resolved_genotypes),
        "births_Cas9_g1_g2": births.tolist(),
        "losses_Cas9_g1_g2": losses.tolist(),
        "terminal_Cas9_g1_g2": final_total.tolist(),
        "molecular_ledger_error": ledger_error.tolist(),
        **counters,
    }
