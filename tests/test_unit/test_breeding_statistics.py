import pytest

from oscar.breeding_statistics import Genotype, generate_breeding_schemes


@pytest.mark.parametrize(
    "n_mutations, expected_schemes",
    [
        pytest.param(
            1,
            [
                [(Genotype.WT,), (Genotype.HOM,)],
                [(Genotype.WT,), (Genotype.HET,)],
                [(Genotype.HOM,), (Genotype.HOM,)],
                [(Genotype.HOM,), (Genotype.HET,)],
                [(Genotype.HET,), (Genotype.HET,)],
            ],
            id="one mutation",
        ),
        pytest.param(
            2,
            [
                [(Genotype.WT, Genotype.WT), (Genotype.HOM, Genotype.HOM)],
                [
                    (Genotype.WT, Genotype.WT),
                    (
                        Genotype.HOM,
                        Genotype.HET,
                    ),
                ],
                [
                    (Genotype.WT, Genotype.WT),
                    (
                        Genotype.HET,
                        Genotype.HET,
                    ),
                ],
                [
                    (
                        Genotype.WT,
                        Genotype.HOM,
                    ),
                    (Genotype.HOM, Genotype.HOM),
                ],
                [
                    (
                        Genotype.WT,
                        Genotype.HOM,
                    ),
                    (Genotype.HOM, Genotype.HET),
                ],
                [
                    (
                        Genotype.WT,
                        Genotype.HOM,
                    ),
                    (Genotype.HET, Genotype.HET),
                ],
                [
                    (
                        Genotype.WT,
                        Genotype.HET,
                    ),
                    (Genotype.HET, Genotype.HET),
                ],
                [
                    (
                        Genotype.HOM,
                        Genotype.HOM,
                    ),
                    (Genotype.HOM, Genotype.HOM),
                ],
                [
                    (
                        Genotype.HOM,
                        Genotype.HOM,
                    ),
                    (Genotype.HOM, Genotype.HET),
                ],
                [
                    (
                        Genotype.HOM,
                        Genotype.HOM,
                    ),
                    (Genotype.HET, Genotype.HET),
                ],
                [
                    (
                        Genotype.HOM,
                        Genotype.HET,
                    ),
                    (Genotype.HET, Genotype.HET),
                ],
                [
                    (
                        Genotype.HET,
                        Genotype.HET,
                    ),
                    (Genotype.HET, Genotype.HET),
                ],
            ],
            id="two mutations",
        ),
    ],
)
def test_generate_breeding_schemes(
    n_mutations: int, expected_schemes: list[list[tuple]]
):
    schemes = generate_breeding_schemes(2)
    print(schemes)
