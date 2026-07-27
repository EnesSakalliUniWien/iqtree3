"""Small SVG helpers and shared styles for SatuTe diagnostic figures."""


SATURATION_FORMULAS = (
    "dominant",
    "eigenvalue_weighted",
)

SATURATION_FORMULA_STYLE = {
    "dominant": ("Published SatuTe", "#0A3FBA", "", "circle"),
    "eigenvalue_weighted": (
        "Eigenvalue-weighted SatuTe",
        "#E69504",
        "9,5",
        "square",
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
