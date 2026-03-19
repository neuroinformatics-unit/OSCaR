"""Tests for the OSCaR frontend Flask app (issue #6)."""

import json
import math

import pytest

from oscar.frontend.app import app as flask_app


@pytest.fixture
def client():
    """Create a test client for the Flask app."""
    flask_app.config["TESTING"] = True
    with flask_app.test_client() as client:
        yield client


class TestIndexRoute:
    def test_index_returns_200(self, client):
        response = client.get("/")
        assert response.status_code == 200

    def test_index_returns_html(self, client):
        response = client.get("/")
        assert b"OSCaR" in response.data


class TestMendelianRatioEndpoint:
    def _post(self, client, parent_1, parent_2):
        return client.post(
            "/api/mendelian_ratio",
            data=json.dumps({"parent_1": parent_1, "parent_2": parent_2}),
            content_type="application/json",
        )

    def test_het_x_het_ratios(self, client):
        """HET x HET should give 25% HOM, 50% HET, 25% WT."""
        response = self._post(client, "het", "het")
        assert response.status_code == 200
        data = response.get_json()
        assert data["ratios"]["hom"] == pytest.approx(0.25)
        assert data["ratios"]["het"] == pytest.approx(0.5)
        assert data["ratios"]["wt"] == pytest.approx(0.25)

    def test_hom_x_hom_gives_only_hom(self, client):
        """HOM x HOM should give 100% HOM offspring."""
        response = self._post(client, "hom", "hom")
        data = response.get_json()
        assert data["ratios"] == {"hom": 1.0}
        assert data["target_ratio"] == 1.0

    def test_wt_x_wt_gives_only_wt(self, client):
        response = self._post(client, "wt", "wt")
        data = response.get_json()
        assert data["ratios"] == {"wt": 1.0}
        assert data["target_ratio"] == 0.0

    def test_ratios_sum_to_one(self, client):
        for scheme in [("het", "het"), ("hom", "het"), ("wt", "het")]:
            response = self._post(client, *scheme)
            data = response.get_json()
            assert sum(data["ratios"].values()) == pytest.approx(1.0)

    def test_target_ratio_and_animals_per_target_consistent(self, client):
        """animals_per_target should equal 1 / target_ratio."""
        response = self._post(client, "het", "het")
        data = response.get_json()
        assert data["animals_per_target"] == pytest.approx(
            1 / data["target_ratio"]
        )

    def test_surplus_per_target_consistent(self, client):
        """surplus_per_target == animals_per_target - 1 (issue #14 fix)."""
        response = self._post(client, "het", "het")
        data = response.get_json()
        assert data["surplus_per_target"] == pytest.approx(
            data["animals_per_target"] - 1
        )

    def test_two_mutations(self, client):
        response = self._post(client, "het_het", "het_het")
        data = response.get_json()
        assert data["n_mutations"] == 2
        assert data["ratios"]["hom_hom"] == pytest.approx(0.0625)

    def test_invalid_genotype_returns_400(self, client):
        response = self._post(client, "invalid", "het")
        assert response.status_code == 400

    def test_mismatched_lengths_returns_400(self, client):
        response = self._post(client, "het_het", "het")
        assert response.status_code == 400


class TestSurplusEndpoint:
    def _post(self, client, **kwargs):
        defaults = {
            "parent_1": "het",
            "parent_2": "het",
            "litter_size": 5.58,
            "n_required": 20,
        }
        defaults.update(kwargs)
        return client.post(
            "/api/surplus",
            data=json.dumps(defaults),
            content_type="application/json",
        )

    def test_returns_200(self, client):
        assert self._post(client).status_code == 200

    def test_matings_needed_is_positive_integer(self, client):
        data = self._post(client).get_json()
        assert isinstance(data["matings_needed"], int)
        assert data["matings_needed"] > 0

    def test_surplus_int_is_ceil_of_surplus_float(self, client):
        """Issue #14: surplus_int must equal math.ceil(surplus_float)."""
        data = self._post(client).get_json()
        assert data["surplus_int"] == math.ceil(data["surplus_float"])

    def test_surplus_is_not_per_genotype_ceil(self, client):
        """Issue #14: surplus_int must NOT be sum of per-genotype ceils."""
        data = self._post(client).get_json()
        per_genotype_ceil_sum = sum(
            info["expected_ceil"]
            for info in data["breakdown"].values()
        )
        # The two approaches should differ (this is the whole point of #14)
        required = 20
        per_geno_surplus = per_genotype_ceil_sum - required
        single_ceil_surplus = data["surplus_int"]
        # Single-ceil surplus <= per-genotype-ceil surplus
        assert single_ceil_surplus <= per_geno_surplus

    def test_impossible_cross_returns_400(self, client):
        response = self._post(client, parent_1="wt", parent_2="wt")
        assert response.status_code == 400

    def test_breakdown_ratios_sum_to_one(self, client):
        data = self._post(client).get_json()
        total = sum(info["ratio"] for info in data["breakdown"].values())
        assert total == pytest.approx(1.0)


class TestSchemesEndpoint:
    def test_one_mutation_returns_schemes(self, client):
        response = client.get("/api/schemes?n_mutations=1")
        data = response.get_json()
        assert response.status_code == 200
        assert data["total"] > 0
        assert len(data["schemes"]) == data["total"]

    def test_schemes_sorted_by_efficiency(self, client):
        data = client.get("/api/schemes?n_mutations=1").get_json()
        scores = [
            s["animals_per_target"]
            for s in data["schemes"]
            if s["animals_per_target"] is not None
        ]
        assert scores == sorted(scores)

    def test_best_scheme_for_one_mutation_is_hom_x_hom(self, client):
        """HOM x HOM should always be the most efficient for HOM target."""
        data = client.get("/api/schemes?n_mutations=1&target=hom").get_json()
        best = data["schemes"][0]
        assert best["animals_per_target"] == pytest.approx(1.0)

    def test_invalid_n_mutations_returns_400(self, client):
        response = client.get("/api/schemes?n_mutations=10")
        assert response.status_code == 400

    def test_two_mutations(self, client):
        data = client.get("/api/schemes?n_mutations=2").get_json()
        assert data["total"] > 0
        assert all(s["parent_1"].count("_") == 1 for s in data["schemes"])