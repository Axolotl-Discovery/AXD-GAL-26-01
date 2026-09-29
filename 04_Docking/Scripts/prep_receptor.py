CHAINS = {"R", "M"}
with open("01_Receptor/8EF5 (2).pdb") as f:
    lineas = f.readlines()
out = [l for l in lineas if (l[:4]=="ATOM" and l[21] in CHAINS) or l[:3] in ("TER","END")]
with open("01_Receptor/receptor_clean.pdb", "w") as f:
    f.writelines(out)
print(f"Guardado: {sum(1 for l in out if l[:4]=='ATOM')} átomos ATOM")
