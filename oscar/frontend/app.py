"""Simple Flask frontend for interacting with the OSCaR package."""

from flask import Flask, jsonify, render_template, request

from oscar.breeding_scheme import BreedingScheme, Genotype, generate_breeding_schemes

app = Flask(__name__)


@app.route("/")
def index():
    """Render the main frontend page."""
    return render_template("index.html")


@app.route("/api/mendelian_ratio", methods=["POST"])
def mendelian_ratio():
    """Calculate the Mendelian ratio for a given breeding scheme.

    Expects JSON body:
        {
            "parent_1": "het_het",
            "parent_2": "het_het"
        }

    Returns
    -------
    JSON
        A dict with genotype strings as keys and ratios as values,
        plus summary statistics (target_ratio, surplus_per_target).
    """
    data = request.get_json()
    parent_1_str = data.get("parent_1", "het")
    parent_2_str = data.get("parent_2", "het")

    try:
        scheme = BreedingScheme(parent_1_str, parent_2_str)
    except (ValueError, KeyError) as e:
        return jsonify({"error": str(e)}), 400

    ratios = scheme.mendelian_ratio()

    # Convert tuple keys to readable strings
    result = {}
    for genotype_tuple, ratio in ratios.items():
        key = "_".join(g.name.lower() for g in genotype_tuple)
        result[key] = round(ratio, 6)

    # Target is all-HOM (most common experimental target)
    target_tuple = tuple(Genotype.HOM for _ in scheme.parent_1_genotype)
    target_ratio = ratios.get(target_tuple, 0.0)

    return jsonify(
        {
            "ratios": result,
            "scheme": repr(scheme),
            "n_mutations": scheme.n_mutations,
            "target_ratio": round(target_ratio, 6),
            "animals_per_target": round(1 / target_ratio, 4) if target_ratio > 0 else None,
            "surplus_per_target": round((1 - target_ratio) / target_ratio, 4)
            if target_ratio > 0
            else None,
        }
    )


@app.route("/api/surplus", methods=["POST"])
def surplus():
    """Calculate surplus animals for a breeding scheme.

    Expects JSON body:
        {
            "parent_1": "het",
            "parent_2": "het",
            "litter_size": 5.58,
            "n_required": 20
        }

    Returns
    -------
    JSON
        Matings needed, total offspring, surplus (using single ceil on total
        per Issue #14 fix), and per-genotype breakdown.
    """
    import math

    data = request.get_json()
    parent_1_str = data.get("parent_1", "het")
    parent_2_str = data.get("parent_2", "het")
    litter_size = float(data.get("litter_size", 5.58))
    n_required = int(data.get("n_required", 20))

    try:
        scheme = BreedingScheme(parent_1_str, parent_2_str)
    except (ValueError, KeyError) as e:
        return jsonify({"error": str(e)}), 400

    ratios = scheme.mendelian_ratio()
    target_tuple = tuple(Genotype.HOM for _ in scheme.parent_1_genotype)
    target_ratio = ratios.get(target_tuple, 0.0)

    if target_ratio == 0:
        return jsonify({"error": "This cross cannot produce HOM offspring."}), 400

    matings_needed = math.ceil(n_required / (litter_size * target_ratio))
    total_offspring_float = litter_size * matings_needed

    # Issue #14 fix: single ceil on total, not per-genotype
    surplus_float = total_offspring_float - n_required
    surplus_int = math.ceil(surplus_float)

    breakdown = {}
    for genotype_tuple, ratio in ratios.items():
        key = "_".join(g.name.lower() for g in genotype_tuple)
        expected = total_offspring_float * ratio
        breakdown[key] = {
            "ratio": round(ratio, 6),
            "expected_float": round(expected, 2),
            "expected_ceil": math.ceil(expected),
        }

    return jsonify(
        {
            "matings_needed": matings_needed,
            "total_offspring_float": round(total_offspring_float, 2),
            "surplus_float": round(surplus_float, 4),
            "surplus_int": surplus_int,
            "breakdown": breakdown,
            "scheme": repr(scheme),
        }
    )


@app.route("/api/schemes", methods=["GET"])
def schemes():
    """Return all valid breeding schemes for a given number of mutations.

    Query params:
        n_mutations (int): number of mutations (default 1)
        target (str): target genotype string e.g. hom or het (default hom)

    Returns
    -------
    JSON
        List of schemes with their Mendelian ratios, ranked by efficiency.
    """
    n_mutations = int(request.args.get("n_mutations", 1))
    target_str = request.args.get("target", "hom")

    if n_mutations < 1 or n_mutations > 4:
        return jsonify({"error": "n_mutations must be between 1 and 4"}), 400

    all_schemes = generate_breeding_schemes(n_mutations)

    target_tuple = tuple(
        Genotype[target_str.upper()] for _ in range(n_mutations)
    )

    results = []
    for scheme in all_schemes:
        ratios = scheme.mendelian_ratio()
        target_ratio = ratios.get(target_tuple, 0.0)
        animals_per_target = (1 / target_ratio) if target_ratio > 0 else None
        surplus_per_target = (
            (1 - target_ratio) / target_ratio if target_ratio > 0 else None
        )

        ratios_serialisable = {
            "_".join(g.name.lower() for g in k): round(v, 6)
            for k, v in ratios.items()
        }

        results.append(
            {
                "scheme": repr(scheme),
                "parent_1": "_".join(
                    g.name.lower() for g in scheme.parent_1_genotype
                ),
                "parent_2": "_".join(
                    g.name.lower() for g in scheme.parent_2_genotype
                ),
                "ratios": ratios_serialisable,
                "target_ratio": round(target_ratio, 6),
                "animals_per_target": round(animals_per_target, 4)
                if animals_per_target
                else None,
                "surplus_per_target": round(surplus_per_target, 4)
                if surplus_per_target
                else None,
            }
        )

    # Sort by efficiency (ascending animals_per_target; impossible last)
    results.sort(
        key=lambda x: x["animals_per_target"]
        if x["animals_per_target"] is not None
        else float("inf")
    )

    return jsonify({"schemes": results, "total": len(results)})


def run():
    """Launch the OSCaR frontend in a browser."""
    import webbrowser

    webbrowser.open("http://localhost:5000")
    app.run(debug=False, port=5000)
