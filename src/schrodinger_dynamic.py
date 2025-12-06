import numpy as np
import json
from scipy.special import sph_harm, factorial
from math import sqrt, exp
import scipy.special as sp

ħ = 1.054571817e-34    # reduced Planck constant
me = 9.10938356e-31    # electron mass
a0 = 1.0               # scale = 1 Bohr radius

# Radial part R_nl(r)
def R(n, l, r):
    rho = 2.0 * r / (n * a0)
    norm = sqrt((2.0 / (n * a0))**3 *
                factorial(n - l - 1) /
                (2 * n * factorial(n + l)))
    L = sp.assoc_laguerre(rho, n - l - 1, 2 * l + 1)
    return norm * np.exp(-rho / 2.0) * rho**l * L

# |ψ|^2
def prob(n, l, m, r, theta, phi):
    return abs(R(n, l, r) * sph_harm(m, l, phi, theta))**2

# Unit vector in +phi direction
def phi_hat(theta, phi):
    return np.array([-np.sin(phi), np.cos(phi), 0])

def sample_points(n, l, m, N=20000):
    points = []
    max_prob = None

    for _ in range(N):
        r = np.random.exponential(scale=n**2)
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2*np.pi)

        P = prob(n,l,m,r,theta,phi)

        if max_prob is None or P > max_prob:
            max_prob = P

        if np.random.rand() < P / max_prob:
            x = r*np.sin(theta)*np.cos(phi)
            y = r*np.sin(theta)*np.sin(phi)
            z = r*np.cos(theta)

            if abs(np.sin(theta)) < 1e-6:
                continue

            # Quantum current velocity (simplified)
            vel_mag = (ħ/me) * (m / (r*np.sin(theta)))
            v = vel_mag * phi_hat(theta, phi)

            points.append({
                "pos": [float(x), float(y), float(z)],
                "vel": [float(v[0]), float(v[1]), float(v[2])],
                "prob": float(P)
            })

    return points


def main():
    n = int(input("n: "))
    l = int(input(f"l (0..{n-1}): "))
    m = int(input(f"m (-{l}..{l}): "))
    N = int(input("Samples (ex: 30000): "))

    pts = sample_points(n,l,m,N)

    import os
    os.makedirs("orbitals_dynamic", exist_ok=True)
    filename = f"orbitals_dynamic/orbital_n{n}_l{l}_m{m}.json"

    with open(filename, "w") as f:
        json.dump({"n":n, "l":l, "m":m, "points":pts}, f)

    print(f"Saved {len(pts)} to {filename}")

if __name__ == "__main__":
    main()
