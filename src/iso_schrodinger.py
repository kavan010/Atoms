import numpy as np
from scipy.special import sph_harm, genlaguerre
from skimage import measure
import math;
import json

# ---------------- Physics ---------------- #
def radial_wavefunction(n, l, r):
    """Hydrogen radial wavefunction (in atomic units a0=1)"""
    rho = 2 * r / n
    # normalization factor
    if n - l - 1 < 0:
        raise ValueError(f"Invalid quantum numbers: n={n}, l={l}. Ensure n > l.")
    prefac = np.sqrt((2/n)**3 * math.factorial(n-l-1) / (2*n*math.factorial(n+l)))
    # Laguerre polynomial
    L = genlaguerre(n-l-1, 2*l+1)(rho)
    return prefac * np.exp(-r/n) * rho**l * L

def real_spherical_harm(l, m, theta, phi):
    """Compute real spherical harmonics for visualization"""
    if m == 0:
        return sph_harm(0, l, phi, theta).real
    elif m > 0:
        return np.sqrt(2) * (-1)**m * sph_harm(m, l, phi, theta).real
    else:  # m < 0
        return np.sqrt(2) * (-1)**m * sph_harm(-m, l, phi, theta).imag

def psi_squared(n, l, m, x, y, z):
    """Compute |psi|^2 at given Cartesian points"""
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.divide(z, r, out=np.zeros_like(z), where=r!=0))
    phi = np.arctan2(y, x)
    R = radial_wavefunction(n, l, r)
    Y = real_spherical_harm(l, m, theta, phi)
    psi = R * Y
    return psi**2

# ---------------- Grid ---------------- #
def generate_grid(grid_size=100, extent=20.0):
    """Create 3D Cartesian grid from -extent to +extent"""
    lin = np.linspace(-extent, extent, grid_size)
    X, Y, Z = np.meshgrid(lin, lin, lin, indexing='ij')
    return X, Y, Z

# ---------------- Isosurface ---------------- #
def extract_isosurface(P, level=0.01):
    """Extract mesh from 3D probability array using marching cubes"""
    verts, faces, _, _ = measure.marching_cubes(P, level=level)
    return verts, faces

# ---------------- Export ---------------- #
def export_mesh(vertices, faces, filename="mesh.json"):
    """Save mesh as JSON"""
    data = {"vertices": vertices.tolist(), "faces": faces.tolist()}
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Mesh exported to {filename}")

# ---------------- Main ---------------- #
if __name__ == "__main__":
    print("\n=== Hydrogen Orbital Generator ===\n")
    n = int(input("Enter n (1..7): "))
    l = int(input(f"Enter l (0..{n-1}): "))
    m = int(input(f"Enter m (-{l}..{l}): "))

    print("Generating grid...")
    X, Y, Z = generate_grid(grid_size=120, extent=20.0)

    print("Computing probability density...")
    P = psi_squared(n, l, m, X, Y, Z)

    # normalize probability to max=1 for easy thresholding
    P /= P.max()

    iso_level = 0.1
    print(f"Extracting isosurface at level {iso_level}...")
    verts, faces = extract_isosurface(P, level=iso_level)

    print(f"Vertices: {len(verts)}, Faces: {len(faces)}")
    export_mesh(verts, faces)
