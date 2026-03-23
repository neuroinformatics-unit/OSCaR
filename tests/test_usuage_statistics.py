"""
Tests for animal usage statistics (issue #3).

Covers:
  - SurplusType.from_sacrifice_reason mapping
  - _surplus_stats_from_series counts and proportions
  - _calculate_usage_statistics per-genotype and overall breakdown
  - Integration via calculate_historical_stats_for_line
"""

import pytest
import pandas as pd

from oscar.historical_stats import (
    SurplusType,
    LineUsageStatistics,
    _surplus_stats_from_series,
    _calculate_usage_statistics,
    calculate_historical_stats_for_line,
)
from oscar.breeding_scheme import Genotype


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def minimal_line_df():
    """Minimal standardised dataframe for a single line with two genotypes."""
    return pd.DataFrame(
        {
            "ID_offspring": range(10),
            "line_name": ["LineA"] * 10,
            "date_of_birth": ["2024-01-01"] * 5 + ["2024-01-02"] * 5,
            "ID_father": ["F1"] * 10,
            "ID_mother": ["M1"] * 10,
            "genotype_father": ["het"] * 10,
            "genotype_mother": ["het"] * 10,
            "genotype_offspring": (
                ["wt"] * 3 + ["het"] * 5 + ["hom"] * 2
            ),
            "sacrifice_reason": (
                # wt × 3
                ["surplus", "surplus", "found dead"]
                # het × 5
                + ["experimental", "experimental", "surplus", "illness", None]
                # hom × 2
                + ["used for experiment", "unknown_reason"]
            ),
            "n_mutations": [1] * 10,
            "mutations": ["GeneX"] * 10,
        }
    )


# ---------------------------------------------------------------------------
# SurplusType.from_sacrifice_reason
# ---------------------------------------------------------------------------

class TestSurplusTypeMapping:

    def test_experimental_reasons(self):
        for reason in ("experimental", "Experimental", "terminal",
                       "used for experiment", "scientific procedure"):
            assert SurplusType.from_sacrifice_reason(reason) == SurplusType.EXPERIMENTAL

    def test_unavoidable_reasons(self):
        for reason in ("found dead", "died", "illness",
                       "humane endpoint", "welfare"):
            assert SurplusType.from_sacrifice_reason(reason) == SurplusType.UNAVOIDABLE

    def test_avoidable_reasons(self):
        for reason in ("surplus", "cull", "culled",
                       "excess", "wrong genotype", "not required"):
            assert SurplusType.from_sacrifice_reason(reason) == SurplusType.AVOIDABLE

    def test_unknown_reason_maps_to_other(self):
        assert SurplusType.from_sacrifice_reason("unknown_reason") == SurplusType.OTHER

    def test_none_maps_to_other(self):
        assert SurplusType.from_sacrifice_reason(None) == SurplusType.OTHER

    def test_nan_maps_to_other(self):
        import math
        assert SurplusType.from_sacrifice_reason(float("nan")) == SurplusType.OTHER

    def test_case_insensitive(self):
        assert SurplusType.from_sacrifice_reason("SURPLUS") == SurplusType.AVOIDABLE
        assert SurplusType.from_sacrifice_reason("  Experimental  ") == SurplusType.EXPERIMENTAL


# ---------------------------------------------------------------------------
# _surplus_stats_from_series
# ---------------------------------------------------------------------------

class TestSurplusStatsFromSeries:

    def test_all_four_types_always_present(self):
        series = pd.Series([SurplusType.EXPERIMENTAL])
        stats = _surplus_stats_from_series(series)
        assert set(stats.n_per_surplus_type.keys()) == set(SurplusType)

    def test_counts_are_correct(self):
        series = pd.Series(
            [SurplusType.EXPERIMENTAL] * 4 +
            [SurplusType.AVOIDABLE] * 3 +
            [SurplusType.UNAVOIDABLE] * 2 +
            [SurplusType.OTHER] * 1
        )
        stats = _surplus_stats_from_series(series)
        assert stats.n_per_surplus_type[SurplusType.EXPERIMENTAL] == 4
        assert stats.n_per_surplus_type[SurplusType.AVOIDABLE] == 3
        assert stats.n_per_surplus_type[SurplusType.UNAVOIDABLE] == 2
        assert stats.n_per_surplus_type[SurplusType.OTHER] == 1
        assert stats.total_n == 10

    def test_proportions_sum_to_one(self):
        series = pd.Series(
            [SurplusType.EXPERIMENTAL] * 6 +
            [SurplusType.AVOIDABLE] * 4
        )
        stats = _surplus_stats_from_series(series)
        total_pct = sum(stats.proportion_per_surplus_type.values())
        assert abs(total_pct - 1.0) < 1e-9

    def test_proportion_values_correct(self):
        series = pd.Series(
            [SurplusType.EXPERIMENTAL] * 3 +
            [SurplusType.AVOIDABLE] * 1
        )
        stats = _surplus_stats_from_series(series)
        assert abs(stats.proportion_per_surplus_type[SurplusType.EXPERIMENTAL] - 0.75) < 1e-9
        assert abs(stats.proportion_per_surplus_type[SurplusType.AVOIDABLE] - 0.25) < 1e-9


# ---------------------------------------------------------------------------
# _calculate_usage_statistics
# ---------------------------------------------------------------------------

class TestCalculateUsageStatistics:

    def test_returns_line_usage_statistics(self, minimal_line_df):
        result = _calculate_usage_statistics(minimal_line_df)
        assert isinstance(result, LineUsageStatistics)

    def test_per_genotype_keys_match_offspring_genotypes(self, minimal_line_df):
        result = _calculate_usage_statistics(minimal_line_df)
        expected_genotypes = {
            (Genotype.WT,), (Genotype.HET,), (Genotype.HOM,)
        }
        assert set(result.surplus_per_genotype.keys()) == expected_genotypes

    def test_overall_total_equals_all_offspring(self, minimal_line_df):
        result = _calculate_usage_statistics(minimal_line_df)
        assert result.overall_surplus.total_n == len(minimal_line_df)

    def test_per_genotype_totals_sum_to_overall(self, minimal_line_df):
        result = _calculate_usage_statistics(minimal_line_df)
        genotype_total = sum(
            s.total_n for s in result.surplus_per_genotype.values()
        )
        assert genotype_total == result.overall_surplus.total_n

    def test_wt_surplus_counts(self, minimal_line_df):
        # wt rows: surplus, surplus, found dead
        result = _calculate_usage_statistics(minimal_line_df)
        wt_stats = result.surplus_per_genotype[(Genotype.WT,)]
        assert wt_stats.n_per_surplus_type[SurplusType.AVOIDABLE] == 2
        assert wt_stats.n_per_surplus_type[SurplusType.UNAVOIDABLE] == 1
        assert wt_stats.total_n == 3

    def test_het_surplus_counts(self, minimal_line_df):
        # het rows: experimental, experimental, surplus, illness, None→OTHER
        result = _calculate_usage_statistics(minimal_line_df)
        het_stats = result.surplus_per_genotype[(Genotype.HET,)]
        assert het_stats.n_per_surplus_type[SurplusType.EXPERIMENTAL] == 2
        assert het_stats.n_per_surplus_type[SurplusType.AVOIDABLE] == 1
        assert het_stats.n_per_surplus_type[SurplusType.UNAVOIDABLE] == 1
        assert het_stats.n_per_surplus_type[SurplusType.OTHER] == 1
        assert het_stats.total_n == 5

    def test_overall_proportions_sum_to_one(self, minimal_line_df):
        result = _calculate_usage_statistics(minimal_line_df)
        total = sum(result.overall_surplus.proportion_per_surplus_type.values())
        assert abs(total - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# Integration: calculate_historical_stats_for_line
# ---------------------------------------------------------------------------

class TestIntegration:

    def test_usage_statistics_present_in_line_stats(self, minimal_line_df):
        stats = calculate_historical_stats_for_line(minimal_line_df, "LineA")
        assert isinstance(stats.usage_statistics, LineUsageStatistics)

    def test_overall_surplus_total_matches_offspring_count(self, minimal_line_df):
        stats = calculate_historical_stats_for_line(minimal_line_df, "LineA")
        assert stats.usage_statistics.overall_surplus.total_n == stats.total_n_offspring

    def test_all_surplus_types_present_in_overall(self, minimal_line_df):
        stats = calculate_historical_stats_for_line(minimal_line_df, "LineA")
        assert set(stats.usage_statistics.overall_surplus.n_per_surplus_type.keys()) == set(SurplusType)
