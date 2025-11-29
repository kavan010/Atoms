# ------- import libs -------
import numpy as np
import json
from scipy.special import sph_harm, factorial
from math import sqrt, exp
import scipy.special as sp

a0 = 1.0

# Radial part R_{n,l}(r) for hydrogen (normalized)
def R(n, l, r):
    # using associated Laguerre polynomials
    rho = 2.0 * r / (n * a0)

    # Normalization constant
    norm = sqrt((2.0 / (n * a0))**3 * factorial(n - l - 1) /
                (2 * n * factorial(n + l)))

    L = sp.assoc_laguerre(rho, n - l - 1, 2 * l + 1)

    return norm * np.exp(-rho / 2.0) * rho**l * L


# Probability density = |ψ|²
def probability(n, l, m, r, theta, phi):
    Y = sph_harm(m, l, phi, theta)  # note: phi, theta
    return (abs(R(n, l, r) * Y))**2


def sample_points(n, l, m, N=50000):
    points = []

    # Rejection sampling
    max_prob = None

    for _ in range(N):
        # Sample position in a sphere ~ few Bohr radii
        r = np.random.exponential(scale=n**2)  # hydrogen radial scale
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2 * np.pi)

        p = probability(n, l, m, r, theta, phi)

        # Establish approximate max probability
        if max_prob is None or p > max_prob:
            max_prob = p

        # Rejection step
        if np.random.rand() < p / max_prob:
            # Convert spherical → cartesian
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            points.append([x, y, z])

    return points


def main():
    print("\n=== Hydrogen Orbital Generator ===\n")
    n = int(input("Enter n (1..7): "))
    l = int(input(f"Enter l (0..{n-1}): "))
    m = int(input(f"Enter m (-{l}..{l}): "))
    N = int(input("How many particle samples? (ex: 50000): "))

    print("\nGenerating samples... please wait...")
    pts = sample_points(n, l, m, N)

    import os
    os.makedirs("orbitals", exist_ok=True)
    filename = f"orbitals/orbital_n{n}_l{l}_m{m}.json"
    with open(filename, "w") as f:
        json.dump({
            "n": n,
            "l": l,
            "m": m,
            "points": pts
        }, f)

    print(f"\nSaved {len(pts)} samples to: {filename}")
    print("\nLoad this JSON in your C++ viewer to visualize!\n")


if __name__ == "__main__":
    main()
