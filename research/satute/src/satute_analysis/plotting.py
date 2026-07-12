"""Small SVG helpers and shared styles for SatuTe diagnostic figures."""


SATURATION_FORMULAS = (
    "dominant",
    "eigenvalue_weighted",
    "mixture_likelihood_weighted",
)

SATURATION_FORMULA_STYLE = {
    "dominant": ("Dominant", "#0A3FBA", "", "circle"),
    "eigenvalue_weighted": (
        "Relative eigenvalue-weighted",
        "#E69504",
        "9,5",
        "square",
    ),
    "mixture_likelihood_weighted": (
        "Soft likelihood mixture",
        "#B6497D",
        "2,4",
        "triangle",
    ),
}


def scale_linear(value, lower, upper, start, length, invert=False):
    if upper <= lower:
        return start + length / 2.0
    fraction = (value - lower) / (upper - lower)
    if invert:
        fraction = 1.0 - fraction
    return start + fraction * length


def svg_path(points):
    if not points:
        return ""
    return "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in points)


def polygon(points):
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
