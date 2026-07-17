import os, time
from pathlib import Path

import numpy as np
from torch import cuda
from ase import units
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.plumed import Plumed
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from mace.calculators import MACECalculator

# ---- config ----
MODEL_PATH   = "../../mace01.model_r2scan0"
INPUT_XYZ    = "gly0ZA-init2.xyz"
TEMPERATURE  = 300.0                # K
NSTEPS       = 500000
DT_FS        = 0.1                  # fs
PRINT_STRIDE = 100
TRAJ_STRIDE  = 100
XYZ_STRIDE   = 100
TDAMP        = 100.0 * units.fs
PLUMED_KERNEL= "/home/k0107/k010734/0-master/plumed/plumed-2.9.1/lib/libplumedKernel.so"

# ---- setup ----
os.environ.setdefault("PLUMED_KERNEL", PLUMED_KERNEL)

xyz_path = Path("trajectory.xyz")
log_path = Path("thermo.log")

device    = "cuda" if cuda.is_available() else "cpu"

# ---- system ----
atoms = read(INPUT_XYZ, index=0).copy()
atoms.set_pbc((True, True, True))

np.random.seed(12345)
MaxwellBoltzmannDistribution(atoms, temperature_K=TEMPERATURE)

calculator = MACECalculator(model_paths=MODEL_PATH, device=device, enable_cueq=True)

# ---- read plumed input ----
PLUMED_FILE = "plumed.dat"
if not os.path.exists(PLUMED_FILE):
    raise FileNotFoundError(f"{PLUMED_FILE} not found!")

plumed_input = [
    line.strip() for line in open(PLUMED_FILE, "r")
    if line.strip() and not line.strip().startswith("#")
]

atoms.calc = Plumed(
    calc=calculator,
    input=plumed_input,
    timestep=DT_FS * units.fs,
    atoms=atoms,
    kT=TEMPERATURE * units.kB,
)

# ---- dynamics ----
dyn = NoseHooverChainNVT(
    atoms=atoms,
    timestep=DT_FS * units.fs,
    temperature_K=TEMPERATURE,
    tdamp=TDAMP,        # thermostat damping (adjust to match desired coupling)
    tchain=3,           # number of thermostats in chain
    tloop=1             # integration loops per step
)

# ---- outputs ----
dyn.attach(
    MDLogger(dyn, atoms, str(log_path), header=True, stress=False, peratom=True, mode="w"),
    interval=PRINT_STRIDE,
)
dyn.attach(lambda: write(str(xyz_path), atoms, format="xyz", append=True), interval=XYZ_STRIDE)

last = time.time()
def status():
    n = len(atoms)
    epot = float(np.asarray(atoms.get_potential_energy())) / n
    ekin = float(np.asarray(atoms.get_kinetic_energy())) / n
    temp = float(np.asarray(atoms.get_temperature()))
    etot = epot + ekin
    now = time.time()
    sim_ps = (PRINT_STRIDE * DT_FS * units.fs) / (1000.0 * units.fs)
    ns_per_day = (sim_ps * 1e-3) / max(now - last, 1e-12) * 86400.0
    print(f"E/atom: Epot={epot:.4f} eV  Ekin={ekin:.4f} eV  T={temp:.0f} K  Etot={etot:.4f} eV | {ns_per_day:.3f} ns/day")
    globals()["last"] = now
dyn.attach(status, interval=PRINT_STRIDE)

# ---- run ----
print(f"Start NVT (Nose-Hoover Chain) @ {int(TEMPERATURE)} K on {device.upper()} | dt={DT_FS:.3f} fs | steps={NSTEPS}")
dyn.run(NSTEPS)
print("Done")
