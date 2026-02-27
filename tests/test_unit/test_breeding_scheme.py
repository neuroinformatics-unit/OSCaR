import pytest

from oscar.breeding_scheme import (
    BreedingScheme,
    Genotype,
    generate_breeding_schemes,
)


@pytest.mark.parametrize(
    "n_mutations, expected_schemes",
    [
        pytest.param(
            1,
            [
                BreedingScheme((Genotype.WT,), (Genotype.HOM,)),
                BreedingScheme((Genotype.WT,), (Genotype.HET,)),
                BreedingScheme((Genotype.HOM,), (Genotype.HOM,)),
                BreedingScheme((Genotype.HOM,), (Genotype.HET,)),
                BreedingScheme((Genotype.HET,), (Genotype.HET,)),
            ],
            id="one mutation",
        ),
        pytest.param(
            2,
            [
                BreedingScheme(
                    (Genotype.WT, Genotype.WT), (Genotype.HOM, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.WT), (Genotype.HOM, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.WT), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.WT), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HOM), (Genotype.HOM, Genotype.WT)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HOM), (Genotype.HOM, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HOM), (Genotype.HOM, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HOM), (Genotype.HET, Genotype.WT)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HOM), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HOM), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HET), (Genotype.HOM, Genotype.WT)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HET), (Genotype.HOM, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HET), (Genotype.HOM, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HET), (Genotype.HET, Genotype.WT)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HET), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.WT, Genotype.HET), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.WT), (Genotype.HOM, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.WT), (Genotype.HOM, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.WT), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.WT), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HOM), (Genotype.HOM, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HOM), (Genotype.HOM, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HOM), (Genotype.HET, Genotype.WT)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HOM), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HOM), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HET), (Genotype.HOM, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HET), (Genotype.HET, Genotype.WT)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HET), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HOM, Genotype.HET), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HET, Genotype.WT), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HET, Genotype.WT), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HET, Genotype.HOM), (Genotype.HET, Genotype.HOM)
                ),
                BreedingScheme(
                    (Genotype.HET, Genotype.HOM), (Genotype.HET, Genotype.HET)
                ),
                BreedingScheme(
                    (Genotype.HET, Genotype.HET), (Genotype.HET, Genotype.HET)
                ),
            ],
            id="two mutations",
        ),
    ],
)
def test_generate_breeding_schemes(n_mutations, expected_schemes):
    schemes = generate_breeding_schemes(n_mutations)
    assert schemes == expected_schemes


@pytest.mark.parametrize(
    "parent_1_genotype, parent_2_genotype, expected_ratios",
    [
        pytest.param(
            (Genotype.WT,),
            (Genotype.HET,),
            {(Genotype.WT,): 0.5, (Genotype.HET,): 0.5},
            id="WT x HET",
        ),
        pytest.param(
            (Genotype.WT,),
            (Genotype.HOM,),
            {(Genotype.HET,): 1},
            id="WT x HOM",
        ),
        pytest.param(
            (Genotype.HET,),
            (Genotype.HOM,),
            {(Genotype.HOM,): 0.5, (Genotype.HET,): 0.5},
            id="HET x HOM",
        ),
        pytest.param(
            (Genotype.HOM,),
            (Genotype.HOM,),
            {
                (Genotype.HOM,): 1,
            },
            id="HOM x HOM",
        ),
        pytest.param(
            (Genotype.HET,),
            (Genotype.HET,),
            {
                (Genotype.WT,): 0.25,
                (Genotype.HET,): 0.5,
                (Genotype.HOM,): 0.25,
            },
            id="HET x HET",
        ),
    ],
)
def test_mendelian_ratio_single_mutation(
    parent_1_genotype, parent_2_genotype, expected_ratios
):
    scheme = BreedingScheme(
        parent_1_genotype=parent_1_genotype,
        parent_2_genotype=parent_2_genotype,
    )
    assert scheme.mendelian_ratio() == expected_ratios


@pytest.mark.parametrize(
    "parent_1_genotype, parent_2_genotype, expected_ratios",
    [
        pytest.param(
            (Genotype.WT, Genotype.HET),
            (Genotype.HET, Genotype.HOM),
            {
                (Genotype.WT, Genotype.HET): 0.25,
                (Genotype.WT, Genotype.HOM): 0.25,
                (Genotype.HET, Genotype.HET): 0.25,
                (Genotype.HET, Genotype.HOM): 0.25,
            },
            id="WT,HET x HET,HOM",
        ),
        pytest.param(
            (Genotype.HET, Genotype.HET),
            (Genotype.HET, Genotype.HET),
            {
                (Genotype.WT, Genotype.WT): 0.0625,
                (Genotype.WT, Genotype.HET): 0.125,
                (Genotype.WT, Genotype.HOM): 0.0625,
                (Genotype.HET, Genotype.WT): 0.125,
                (Genotype.HET, Genotype.HET): 0.25,
                (Genotype.HET, Genotype.HOM): 0.125,
                (Genotype.HOM, Genotype.WT): 0.0625,
                (Genotype.HOM, Genotype.HET): 0.125,
                (Genotype.HOM, Genotype.HOM): 0.0625,
            },
            id="HET,HET x HET,HET",
        ),
        pytest.param(
            (Genotype.HET, Genotype.HOM),
            (Genotype.HOM, Genotype.HET),
            {
                (Genotype.HET, Genotype.HET): 0.25,
                (Genotype.HET, Genotype.HOM): 0.25,
                (Genotype.HOM, Genotype.HET): 0.25,
                (Genotype.HOM, Genotype.HOM): 0.25,
            },
            id="HET,HOM x HOM,HET",
        ),
        pytest.param(
            (Genotype.HET, Genotype.HET),
            (Genotype.HET, Genotype.WT),
            {
                (Genotype.WT, Genotype.WT): 0.125,
                (Genotype.WT, Genotype.HET): 0.125,
                (Genotype.HET, Genotype.WT): 0.25,
                (Genotype.HET, Genotype.HET): 0.25,
                (Genotype.HOM, Genotype.WT): 0.125,
                (Genotype.HOM, Genotype.HET): 0.125,
            },
            id="HET,HET x HET,WT",
        ),
    ],
)
def test_mendelian_ratio_two_mutations(
    parent_1_genotype, parent_2_genotype, expected_ratios
):
    scheme = BreedingScheme(
        parent_1_genotype=parent_1_genotype,
        parent_2_genotype=parent_2_genotype,
    )
    assert scheme.mendelian_ratio() == expected_ratios


@pytest.mark.parametrize(
    "parent_1_genotype, parent_2_genotype, expected_ratios",
    [
        pytest.param(
            (Genotype.WT, Genotype.HET, Genotype.WT),
            (Genotype.HET, Genotype.HET, Genotype.HET),
            {
                (Genotype.WT, Genotype.WT, Genotype.WT): 0.0625,
                (Genotype.WT, Genotype.WT, Genotype.HET): 0.0625,
                (Genotype.WT, Genotype.HET, Genotype.WT): 0.125,
                (Genotype.WT, Genotype.HET, Genotype.HET): 0.125,
                (Genotype.WT, Genotype.HOM, Genotype.WT): 0.0625,
                (Genotype.WT, Genotype.HOM, Genotype.HET): 0.0625,
                (Genotype.HET, Genotype.WT, Genotype.WT): 0.0625,
                (Genotype.HET, Genotype.WT, Genotype.HET): 0.0625,
                (Genotype.HET, Genotype.HET, Genotype.WT): 0.125,
                (Genotype.HET, Genotype.HET, Genotype.HET): 0.125,
                (Genotype.HET, Genotype.HOM, Genotype.WT): 0.0625,
                (Genotype.HET, Genotype.HOM, Genotype.HET): 0.0625,
            },
            id="WT,HET,WT x HET,HET,HET",
        ),
        pytest.param(
            (Genotype.HET, Genotype.HET, Genotype.HET),
            (Genotype.HET, Genotype.HET, Genotype.HET),
            {
                (Genotype.WT, Genotype.WT, Genotype.WT): 0.015625,
                (Genotype.WT, Genotype.WT, Genotype.HET): 0.03125,
                (Genotype.WT, Genotype.WT, Genotype.HOM): 0.015625,
                (Genotype.WT, Genotype.HET, Genotype.WT): 0.03125,
                (Genotype.WT, Genotype.HET, Genotype.HET): 0.0625,
                (Genotype.WT, Genotype.HET, Genotype.HOM): 0.03125,
                (Genotype.WT, Genotype.HOM, Genotype.WT): 0.015625,
                (Genotype.WT, Genotype.HOM, Genotype.HET): 0.03125,
                (Genotype.WT, Genotype.HOM, Genotype.HOM): 0.015625,
                (Genotype.HET, Genotype.WT, Genotype.WT): 0.03125,
                (Genotype.HET, Genotype.WT, Genotype.HET): 0.0625,
                (Genotype.HET, Genotype.WT, Genotype.HOM): 0.03125,
                (Genotype.HET, Genotype.HET, Genotype.WT): 0.0625,
                (Genotype.HET, Genotype.HET, Genotype.HET): 0.125,
                (Genotype.HET, Genotype.HET, Genotype.HOM): 0.0625,
                (Genotype.HET, Genotype.HOM, Genotype.WT): 0.03125,
                (Genotype.HET, Genotype.HOM, Genotype.HET): 0.0625,
                (Genotype.HET, Genotype.HOM, Genotype.HOM): 0.03125,
                (Genotype.HOM, Genotype.WT, Genotype.WT): 0.015625,
                (Genotype.HOM, Genotype.WT, Genotype.HET): 0.03125,
                (Genotype.HOM, Genotype.WT, Genotype.HOM): 0.015625,
                (Genotype.HOM, Genotype.HET, Genotype.WT): 0.03125,
                (Genotype.HOM, Genotype.HET, Genotype.HET): 0.0625,
                (Genotype.HOM, Genotype.HET, Genotype.HOM): 0.03125,
                (Genotype.HOM, Genotype.HOM, Genotype.WT): 0.015625,
                (Genotype.HOM, Genotype.HOM, Genotype.HET): 0.03125,
                (Genotype.HOM, Genotype.HOM, Genotype.HOM): 0.015625,
            },
            id="HET,HET,HET x HET,HET,HET",
        ),
    ],
)
def test_mendelian_ratio_three_mutations(
    parent_1_genotype, parent_2_genotype, expected_ratios
):
    scheme = BreedingScheme(
        parent_1_genotype=parent_1_genotype,
        parent_2_genotype=parent_2_genotype,
    )
    assert scheme.mendelian_ratio() == expected_ratios


def test_breeding_scheme_equality():
    """
    Check breeding schemes are marked as equal independent of parent 1 vs
    2 order.
    """
    parent_1_genotype = (Genotype.WT, Genotype.HET)
    parent_2_genotype = (Genotype.HOM, Genotype.HOM)

    assert BreedingScheme(
        parent_1_genotype, parent_2_genotype
    ) == BreedingScheme(parent_1_genotype, parent_2_genotype)
    assert BreedingScheme(
        parent_1_genotype, parent_2_genotype
    ) == BreedingScheme(parent_2_genotype, parent_1_genotype)
    assert BreedingScheme(
        parent_1_genotype, parent_2_genotype
    ) != BreedingScheme(parent_2_genotype, parent_2_genotype)
