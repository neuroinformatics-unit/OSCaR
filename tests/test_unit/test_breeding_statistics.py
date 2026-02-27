import pytest

from oscar.breeding_statistics import (
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
                [(Genotype.WT, Genotype.WT), (Genotype.HOM, Genotype.HET)],
                [(Genotype.WT, Genotype.WT), (Genotype.HET, Genotype.HOM)],
                [(Genotype.WT, Genotype.WT), (Genotype.HET, Genotype.HET)],
                [(Genotype.WT, Genotype.HOM), (Genotype.HOM, Genotype.WT)],
                [(Genotype.WT, Genotype.HOM), (Genotype.HOM, Genotype.HOM)],
                [(Genotype.WT, Genotype.HOM), (Genotype.HOM, Genotype.HET)],
                [(Genotype.WT, Genotype.HOM), (Genotype.HET, Genotype.WT)],
                [(Genotype.WT, Genotype.HOM), (Genotype.HET, Genotype.HOM)],
                [(Genotype.WT, Genotype.HOM), (Genotype.HET, Genotype.HET)],
                [(Genotype.WT, Genotype.HET), (Genotype.HOM, Genotype.WT)],
                [(Genotype.WT, Genotype.HET), (Genotype.HOM, Genotype.HOM)],
                [(Genotype.WT, Genotype.HET), (Genotype.HOM, Genotype.HET)],
                [(Genotype.WT, Genotype.HET), (Genotype.HET, Genotype.WT)],
                [(Genotype.WT, Genotype.HET), (Genotype.HET, Genotype.HOM)],
                [(Genotype.WT, Genotype.HET), (Genotype.HET, Genotype.HET)],
                [(Genotype.HOM, Genotype.WT), (Genotype.HOM, Genotype.HOM)],
                [(Genotype.HOM, Genotype.WT), (Genotype.HOM, Genotype.HET)],
                [(Genotype.HOM, Genotype.WT), (Genotype.HET, Genotype.HOM)],
                [(Genotype.HOM, Genotype.WT), (Genotype.HET, Genotype.HET)],
                [(Genotype.HOM, Genotype.HOM), (Genotype.HOM, Genotype.HOM)],
                [(Genotype.HOM, Genotype.HOM), (Genotype.HOM, Genotype.HET)],
                [(Genotype.HOM, Genotype.HOM), (Genotype.HET, Genotype.WT)],
                [(Genotype.HOM, Genotype.HOM), (Genotype.HET, Genotype.HOM)],
                [(Genotype.HOM, Genotype.HOM), (Genotype.HET, Genotype.HET)],
                [(Genotype.HOM, Genotype.HET), (Genotype.HOM, Genotype.HET)],
                [(Genotype.HOM, Genotype.HET), (Genotype.HET, Genotype.WT)],
                [(Genotype.HOM, Genotype.HET), (Genotype.HET, Genotype.HOM)],
                [(Genotype.HOM, Genotype.HET), (Genotype.HET, Genotype.HET)],
                [(Genotype.HET, Genotype.WT), (Genotype.HET, Genotype.HOM)],
                [(Genotype.HET, Genotype.WT), (Genotype.HET, Genotype.HET)],
                [(Genotype.HET, Genotype.HOM), (Genotype.HET, Genotype.HOM)],
                [(Genotype.HET, Genotype.HOM), (Genotype.HET, Genotype.HET)],
                [(Genotype.HET, Genotype.HET), (Genotype.HET, Genotype.HET)],
            ],
            id="two mutations",
        ),
    ],
)
def test_generate_breeding_schemes(
    n_mutations: int, expected_schemes: list[list[tuple]]
):
    schemes = generate_breeding_schemes(n_mutations)
    assert schemes == expected_schemes


@pytest.mark.parametrize(
    "parent_1_genotype, parent_2_genotype, expected_ratios",
    [
        pytest.param(
            (Genotype.WT,),
            (Genotype.HET,),
            {Genotype.WT: 0.5, Genotype.HET: 0.5},
            id="WT x HET",
        ),
        pytest.param(
            (Genotype.WT,),
            (Genotype.HOM,),
            {Genotype.HET: 1},
            id="WT x HOM",
        ),
        pytest.param(
            (Genotype.HET,),
            (Genotype.HOM,),
            {Genotype.HOM: 0.5, Genotype.HET: 0.5},
            id="HET x HOM",
        ),
        pytest.param(
            (Genotype.HOM,),
            (Genotype.HOM,),
            {
                Genotype.HOM: 1,
            },
            id="HOM x HOM",
        ),
        pytest.param(
            (Genotype.HET,),
            (Genotype.HET,),
            {Genotype.WT: 0.25, Genotype.HET: 0.5, Genotype.HOM: 0.25},
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
