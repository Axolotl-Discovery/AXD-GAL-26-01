import numpy as np

# Centro del sitio ortostérico (desde el ligando co-cristalizado)
with open("01_Receptor/8EF5 (2).pdb") as f:
    hetatm = [l for l in f if l[:6].strip()=="HETATM" and l[17:20].strip() not in ("HOH","WAT")]

if hetatm:
    xs = [float(l[30:38]) for l in hetatm]
    ys = [float(l[38:46]) for l in hetatm]
    zs = [float(l[46:54]) for l in hetatm]
    cx, cy, cz = np.mean(xs), np.mean(ys), np.mean(zs)
    print(f"Ligando co-cristalizado: {hetatm[0][17:20].strip()} ({len(hetatm)} átomos)")
    print(f"\nSITIO  → cx={cx:.2f} cy={cy:.2f} cz={cz:.2f}")
else:
    print("No se encontraron HETATM no-agua")

# Centro del receptor completo (para barrido ciego)
with open("01_Receptor/receptor_clean.pdbqt") as f:
    atoms = [l for l in f if l[:4]=="ATOM"]
xs2 = [float(l[30:38]) for l in atoms]
ys2 = [float(l[38:46]) for l in atoms]
zs2 = [float(l[46:54]) for l in atoms]
bcx, bcy, bcz = np.mean(xs2), np.mean(ys2), np.mean(zs2)
sx = max(xs2)-min(xs2)+15; sy = max(ys2)-min(ys2)+15; sz = max(zs2)-min(zs2)+15
print(f"CIEGO  → cx={bcx:.2f} cy={bcy:.2f} cz={bcz:.2f}  size={sx:.0f}x{sy:.0f}x{sz:.0f}")
