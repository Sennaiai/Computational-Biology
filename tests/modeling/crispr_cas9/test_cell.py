from dataclasses import replace
import math
import unittest

import numpy as np

from crispr_kinetics.cell import (
    CellParameters,
    DEFAULT_LOCI,
    Locus,
    MOLECULES_PER_NM_FL,
    molecule_totals,
    partition,
    reactions,
    simulate_founder,
    volume_scale,
)


def silent(**kwargs):
    return replace(
        CellParameters(), cas9_production_nM_h=0, guide_production_nM_h=(0, 0), **kwargs
    )


class CellTests(unittest.TestCase):
    def test_concentration_conversion_and_division_volumes(self):
        self.assertAlmostEqual(MOLECULES_PER_NM_FL * 30, 18.06642228)
        p = CellParameters()
        self.assertEqual(volume_scale(0, 0, p), 1)
        self.assertEqual(volume_scale(p.cycle_h, 0, p), 2)
        self.assertEqual(volume_scale(3 * p.cycle_h, 0, p), 2)

    def test_each_reaction_reconciles_free_and_bound_material(self):
        p = CellParameters()
        pools = np.array([3, 4, 5, 6, 7, 8])
        for state in range(8):
            dna = np.full((3, 2), state, dtype=int)
            before = molecule_totals(pools, dna, DEFAULT_LOCI)
            for hazard, changes, update, birth, loss in reactions(
                pools, dna, DEFAULT_LOCI, p, 30, 150, True
            ):
                self.assertGreater(hazard, 0)
                x, y = pools.copy(), dna.copy()
                for i, delta in changes:
                    x[i] += delta
                if update:
                    y[update[0], update[1]] = update[2]
                expected = before.copy()
                if birth is not None:
                    expected[birth] += 1
                for component in loss or ():
                    expected[component] -= 1
                self.assertTrue(np.all(x >= 0))
                np.testing.assert_array_equal(
                    molecule_totals(x, y, DEFAULT_LOCI), expected
                )

    def test_daughter_partition_conserves_molecules_and_inherits_diploid_edits(self):
        pools = np.array([11, 7, 9, 5, 3, 1])
        dna = np.array([[1, 5], [2, 6], [7, 0]])
        before = molecule_totals(pools, dna, DEFAULT_LOCI)
        children = partition(pools, dna, DEFAULT_LOCI, np.random.default_rng(19))
        np.testing.assert_array_equal(
            sum(molecule_totals(a, b, DEFAULT_LOCI) for a, b in children), before
        )
        np.testing.assert_array_equal(children[0][1], [[0, 5], [0, 6], [7, 0]])
        np.testing.assert_array_equal(children[0][1], children[1][1])
        dna[0, 0] = 4
        with self.assertRaises(ValueError):
            partition(pools, dna, DEFAULT_LOCI, np.random.default_rng(1))

    def test_no_Cas9_or_missing_guide_prevents_corresponding_edits(self):
        zero = replace(
            CellParameters(), cas9_production_nM_h=0, horizon_h=4, max_generations=0
        )
        result = simulate_founder(zero)
        self.assertEqual(result["cut_events"], 0)
        no_g1 = replace(
            CellParameters(),
            guide_production_nM_h=(0, 1),
            horizon_h=4,
            max_generations=0,
        )
        result = simulate_founder(no_g1)
        for terminal in result["terminal_cells"]:
            np.testing.assert_array_equal(np.asarray(terminal["dna"])[[0, 2]], 0)
        no_import = replace(
            CellParameters(), import_per_h=0, horizon_h=4, max_generations=0
        )
        self.assertEqual(simulate_founder(no_import)["cut_events"], 0)

    def test_clonal_founder_edits_are_preserved_without_mosaicism(self):
        p = silent(
            protein_decay_per_h=0,
            guide_decay_per_h=0,
            rnp_decay_per_h=0,
            import_per_h=0,
            export_per_h=0,
            assembly_per_nM_h=0,
            disassembly_per_h=0,
        )
        result = simulate_founder(
            p,
            initial_dna=np.full((3, 2), 5),
            initial_pools=np.array([10, 0, 0, 0, 0, 0]),
            times_h=[0, 8, 16, 24],
        )
        np.testing.assert_array_equal(
            np.asarray(result["observations"])[:, 0], [1, 2, 4, 8]
        )
        np.testing.assert_array_equal(
            np.asarray(result["observations"])[:, 1], [1, 2, 4, 8]
        )
        self.assertEqual(result["terminal_cells_count"], 8)
        self.assertFalse(result["mosaicism_proxy"])
        self.assertEqual(result["terminal_Cas9_g1_g2"], [10, 0, 0])

    def test_growth_dependent_births_match_nonhomogeneous_Poisson_reference(self):
        p = silent(
            horizon_h=2,
            cycle_h=8,
            max_generations=0,
            nucleus_birth_fL=2,
            cytoplasm_birth_fL=2,
            import_per_h=0,
            export_per_h=0,
            protein_decay_per_h=0,
            guide_decay_per_h=0,
            rnp_decay_per_h=0,
            assembly_per_nM_h=0,
            disassembly_per_h=0,
        )
        p = replace(p, cas9_production_nM_h=0.4)
        values = [
            simulate_founder(p, loci=(), seed=i)["terminal_Cas9_g1_g2"][0]
            for i in range(1500)
        ]
        growth = math.log(2) / p.cycle_h
        expected = 0.4 * 2 * MOLECULES_PER_NM_FL * math.expm1(growth * 2) / growth
        self.assertLess(
            abs(np.mean(values) - expected), 6 * math.sqrt(expected / len(values))
        )
        self.assertLess(abs(np.var(values) - expected), 0.2)

    def test_repair_only_matches_competing_exponential_branches(self):
        p = silent(
            horizon_h=3,
            max_generations=0,
            precise_repair_per_h=0.3,
            mutation_repair_per_h=0.2,
        )
        counts = np.zeros(8)
        for seed in range(2000):
            r = simulate_founder(
                p, loci=(Locus("test", 0),), initial_dna=np.array([[4, 7]]), seed=seed
            )
            counts[r["terminal_cells"][0]["dna"][0][0]] += 1
        repaired = 1 - math.exp(-0.5 * 3)
        expected = np.zeros(8)
        expected[[0, 4, 5, 6, 7]] = [
            repaired * 0.6,
            1 - repaired,
            repaired * 0.4 * 0.35,
            repaired * 0.4 * 0.35,
            repaired * 0.4 * 0.3,
        ]
        epsilon = math.sqrt(math.log(16 / 1e-6) / 4000)
        self.assertLess(float(abs(counts / 2000 - expected).max()), epsilon)

    def test_precise_repair_enables_recutting_and_unrepaired_breaks_block_division(
        self,
    ):
        p = silent(
            horizon_h=2,
            max_generations=0,
            precise_repair_per_h=100,
            mutation_repair_per_h=0,
            product_release_per_h=100,
            rnp_decay_per_h=0,
            disassembly_per_h=0,
        )
        locus = Locus(
            "test",
            0,
            unbinding_per_h=0,
            loop_per_h=100,
            unwind_per_h=0,
            cleavage_per_h=100,
        )
        r = simulate_founder(
            p, loci=(locus,), initial_pools=np.array([0, 0, 0, 0, 100, 0]), seed=4
        )
        self.assertGreater(r["cut_events"], 2)
        self.assertGreater(r["precise_repair_events"], 2)
        self.assertEqual(r["mutation_repair_events"], 0)
        blocked = silent(horizon_h=16, precise_repair_per_h=0, mutation_repair_per_h=0)
        r = simulate_founder(blocked, initial_dna=np.full((3, 2), 4))
        self.assertEqual(r["terminal_cells_count"], 1)
        self.assertEqual(r["divisions"], 0)

    def test_reproducibility_budget_and_bad_inputs(self):
        p = silent(horizon_h=1, max_generations=0)
        active = replace(CellParameters(), horizon_h=2, max_generations=0)
        self.assertEqual(
            simulate_founder(active, seed=17), simulate_founder(active, seed=17)
        )
        np.random.seed(41)
        before = np.random.get_state()
        simulate_founder(active, seed=17)
        after = np.random.get_state()
        np.testing.assert_array_equal(before[1], after[1])
        self.assertEqual(before[2:], after[2:])
        with self.assertRaises(RuntimeError):
            simulate_founder(replace(p, cas9_production_nM_h=100), event_budget=1)
        for changes in (
            {"import_per_h": -1},
            {"cycle_h": 0},
            {"horizon_h": float("nan")},
            {"disruptive_given_mutation": 2},
            {"max_generations": True},
        ):
            with self.assertRaises(ValueError):
                CellParameters(**changes)
        for times in ([], [1, 0], [float("nan")], [2], [True]):
            with self.assertRaises(ValueError):
                simulate_founder(p, times_h=times)

    def test_late_stochastic_edits_can_produce_distinct_daughter_genotypes(self):
        parameters = replace(
            CellParameters(),
            expression_on_h=8,
            expression_off_h=20,
            precise_repair_per_h=0,
            mutation_repair_per_h=2,
            product_release_per_h=4,
            max_generations=2,
        )
        result = simulate_founder(
            parameters, loci=(Locus("synthetic_test_locus", 0),), seed=42
        )
        self.assertTrue(result["mosaicism_evaluable"])
        self.assertTrue(result["mosaicism_proxy"])
        self.assertGreater(result["distinct_resolved_synthetic_genotypes"], 1)


if __name__ == "__main__":
    unittest.main()
