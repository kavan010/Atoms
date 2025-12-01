import numpy as np
import json
import os
from math import sqrt
from pyscf import gto, scf

def build_h2o(basis="sto-3g"):
    # water geometry (Angstrom) — simple approx
    # Place O at origin, H atoms at some approximate positions
    mol = gto.M(
        atom = [
            ("O", (0.0, 0.0, 0.0)),
            ("H", ( 0.757,  0.586, 0.0)),
            ("H", (-0.757,  0.586, 0.0)),
        ],
        basis = basis,
        unit = "Angstrom",
        verbose = 0,
    )
    mol.build()
    return mol

def run_hf(mol):
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    return mf

def make_mo_wavefunc(mol, mo_coeff, mo_idx):
    """
    Returns a function psi(r) that evaluates the molecular orbital `mo_idx` at 3D point r = (x,y,z).
    We'll do this by evaluating all atomic basis functions at r, then sum c_i * phi_i(r).
    """
    ao_locations = mol.atom_coords()       # in Angstrom
    ao_basis = mol._basis                   # internal basis info (primitives, exponents, etc.)
    # but simplest: we use mol.eval_gto for AO value evaluation

    def psi(r):
        # r: 3-element array (x,y,z) in Angstrom
        # pyscf expects shape (Npts,3)
        val_ao = mol.eval_gto('GTOval_sph', np.array([r]), blas=True)[0]  # AO values at r
        # mo_coeff: shape (nao, nmo)
        coeffs = mo_coeff[:, mo_idx]
        return np.dot(coeffs, val_ao)

    return psi

def sample_from_density(psi, n_samples=50000, box_radius=5.0, max_tries=500000):
    """
    Simple rejection sampling inside a sphere of radius box_radius (Å).
    Returns list of points (x,y,z).
    """
    samples = []
    max_density = 0.0
    tries = 0

    while len(samples) < n_samples and tries < max_tries:
        # uniform inside sphere
        u = np.random.rand()
        cos_theta = np.random.uniform(-1,1)
        phi = np.random.uniform(0, 2*np.pi)
        r = box_radius * (u ** (1/3))
        sin_theta = np.sqrt(1 - cos_theta**2)

        x = r * sin_theta * np.cos(phi)
        y = r * sin_theta * np.sin(phi)
        z = r * cos_theta
        R = np.array([x, y, z])

        dens = abs(psi(R))**2
        if dens > max_density:
            max_density = dens
        # accept with probability density / max_density
        if np.random.rand() < dens / (max_density + 1e-20):
            samples.append([float(x), float(y), float(z)])
        tries += 1

    return samples

def main():
    print("Building H2O / running HF …")
    mol = build_h2o(basis="sto-3g")
    mf = run_hf(mol)

    nmo = mf.mo_coeff.shape[1]
    print(f"Number of MOs: {nmo}, energies: {mf.mo_energy}")

    # Choose which MO to sample — e.g. HOMO (highest occupied), or user picks index
    mo_idx = int(input(f"Enter MO index to sample (0 .. {nmo-1}): "))
    N = int(input("How many sample points? (e.g. 50000): "))
    box_radius = float(input("Sampling sphere radius (Å, e.g. 5.0): "))

    psi = make_mo_wavefunc(mol, mf.mo_coeff, mo_idx)

    print("Sampling … (this may take a while)")
    pts = sample_from_density(psi, n_samples=N, box_radius=box_radius)

    os.makedirs("molecular_orbitals", exist_ok=True)
    filename = f"molecular_orbitals/h2o_mo{mo_idx}_N{len(pts)}.json"
    with open(filename, "w") as f:
        json.dump({
            "molecule": "H2O",
            "basis": mol.basis,
            "mo_idx": mo_idx,
            "mo_energy": float(mf.mo_energy[mo_idx]),
            "points": pts,
        }, f)

    print(f"Saved {len(pts)} samples → {filename}")

if __name__ == "__main__":
    main()
