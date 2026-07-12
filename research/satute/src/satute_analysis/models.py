"""Reversible substitution-model helpers shared by exact experiments."""

import numpy as np


def parse_model(model):
    if model == "JC":
        return np.ones(6, dtype=float), np.repeat(0.25, 4)
    if model.startswith("F81+F{") and model.endswith("}"):
        pi = np.array([float(value) for value in model[len("F81+F{") : -1].split(",")])
        if len(pi) != 4:
            raise ValueError(f"Unsupported F81 model dimensions: {model}")
        pi /= pi.sum()
        return np.ones(6, dtype=float), pi
    if model.startswith("K2P{") and model.endswith("}+FQ"):
        kappa = float(model[len("K2P{") : -len("}+FQ")])
        return np.array([1.0, kappa, 1.0, 1.0, kappa, 1.0]), np.repeat(0.25, 4)
    if not (model.startswith("GTR{") and "}+F{" in model and model.endswith("}")):
        raise ValueError(f"Unsupported model syntax: {model}")
    rate_text, frequency_text = model[4:-1].split("}+F{", 1)
    rates = np.array([float(value) for value in rate_text.split(",")], dtype=float)
    pi = np.array([float(value) for value in frequency_text.split(",")], dtype=float)
    expected = len(pi) * (len(pi) - 1) // 2
    if len(pi) not in {4, 20} or len(rates) != expected:
        raise ValueError(f"Unsupported model dimensions: {model}")
    pi /= pi.sum()
    return rates, pi


def build_q(model):
    rates, pi = parse_model(model)
    q = np.zeros((len(pi), len(pi)), dtype=float)
    rate_index = 0
    for i in range(len(pi) - 1):
        for j in range(i + 1, len(pi)):
            rate = rates[rate_index]
            rate_index += 1
            q[i, j] = rate * pi[j]
            q[j, i] = rate * pi[i]
    for i in range(len(pi)):
        q[i, i] = -float(np.sum(q[i, :]))
    q /= -float(np.dot(pi, np.diag(q)))
    return q, pi


def reversible_eigendecomposition(q, pi):
    sqrt_pi = np.sqrt(pi)
    inverse_sqrt_pi = 1.0 / sqrt_pi
    symmetric_q = np.diag(sqrt_pi) @ q @ np.diag(inverse_sqrt_pi)
    eigenvalues, basis = np.linalg.eigh(symmetric_q)
    eigenvectors = np.diag(inverse_sqrt_pi) @ basis
    inverse = basis.T @ np.diag(sqrt_pi)
    return eigenvalues, eigenvectors, inverse


def transition_matrix(eigenvalues, eigenvectors, inverse, length):
    transition = eigenvectors @ np.diag(np.exp(eigenvalues * length)) @ inverse
    transition[np.logical_and(transition < 0.0, transition > -1e-14)] = 0.0
    return transition
