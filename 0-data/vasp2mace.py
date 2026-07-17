import os, re
from ase.io import read, write

OUT_FILE = "asp.xyz"
EXCLUDE = {}

if os.path.exists(OUT_FILE):
    os.remove(OUT_FILE)

n_written = 0
for i in range(1, 224):
    if i in EXCLUDE:
        print(f"Skipping spe-{i}: excluded")
        continue
    spe_dir = f"spe-{i}"
    path = f"{spe_dir}/OUTCAR"
    print(f"Processing {spe_dir}...")
    if not os.path.isfile(path):
        print(f"  [Skipping] {spe_dir}: OUTCAR not found")
        continue
    atoms = read(path, format="vasp-out")
    atoms.set_pbc([True, True, True])
    write(OUT_FILE, atoms, format="extxyz", write_results=True, append=True)
    n_written += 1
    print(f"  Appended frame from {spe_dir}")

if n_written:
    with open(OUT_FILE, "r", encoding="utf-8") as f:
        txt = f.read()
    txt = re.sub(r'stress="[^"]*"\s*', "", txt)
    txt = re.sub(r"free_energy=[^\s]*\s*", "", txt)
    with open(OUT_FILE, "w", encoding="utf-8") as f:
        f.write(txt)
    print(f"Created cleaned file: {OUT_FILE}")
else:
    print("No frames written; nothing to clean.")

