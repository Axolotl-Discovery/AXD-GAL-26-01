#!/usr/bin/env python3
"""
Analisis de contactos proteina-ligando — AXD-GAL-26-01
Criterios (documentacion del entorno):
  H-Bond:       <= 3.5 A  entre N/O del receptor y N/O del ligando
  Hidrofobico:  <= 4.0 A  entre C del receptor y C del ligando
  Electrostatico: <= 4.0 A entre atomos cargados del receptor y heteroatomos del ligando
"""
import os, numpy as np
from collections import defaultdict

DIST_HB, DIST_HF, DIST_EL = 3.5, 4.0, 4.0

# Residuos clave del sitio ortostertico MOR (numeracion UniProt P35372)
KEY_RES = {147, 148, 151, 155, 219, 229, 241, 293, 297, 300, 322, 326, 329, 332}

def parse_pdbqt(path, model=1):
    atoms, cur, active = [], 0, False
    with open(path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "MODEL":
                cur += 1; active = (cur == model)
            elif rec == "ENDMDL":
                if active: break
                active = False
            elif rec in ("ATOM","HETATM") and (active or cur == 0):
                try:
                    atoms.append({
                        "name":    line[12:16].strip(),
                        "resname": line[17:20].strip(),
                        "chain":   line[21],
                        "resid":   int(line[22:26]),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "elem":    line[12:16].strip()[0],
                    })
                except: pass
    return atoms

def d(a, b):
    return np.sqrt((a["x"]-b["x"])**2+(a["y"]-b["y"])**2+(a["z"]-b["z"])**2)

def is_charged(a):
    r, n = a["resname"], a["name"]
    return ((r=="ASP" and n in ("OD1","OD2")) or
            (r=="GLU" and n in ("OE1","OE2")) or
            (r=="LYS" and n=="NZ") or
            (r=="ARG" and n in ("NH1","NH2","NE")) or
            (r=="HIS" and n in ("ND1","NE2")))

def contacts(rec, lig):
    result = []
    lig_NO = [a for a in lig if a["elem"] in ("N","O")]
    lig_C  = [a for a in lig if a["elem"]=="C"]
    lig_het= [a for a in lig if a["elem"] in ("N","O","S")]
    seen = set()
    for ra in rec:
        key_base = (ra["chain"], ra["resid"], ra["resname"])
        # H-Bond
        if ra["elem"] in ("N","O"):
            for la in lig_NO:
                dd = d(ra, la)
                if dd <= DIST_HB:
                    k = key_base+("HB",)
                    if k not in seen:
                        result.append({**{kk:ra[kk] for kk in ("chain","resid","resname")},
                                       "tipo":"H-Bond","dist":round(dd,2),
                                       "ar":ra["name"],"al":la["name"]})
                        seen.add(k)
        # Hidrofobico
        if ra["elem"]=="C":
            for la in lig_C:
                dd = d(ra, la)
                if dd <= DIST_HF:
                    k = key_base+("HF",)
                    if k not in seen:
                        result.append({**{kk:ra[kk] for kk in ("chain","resid","resname")},
                                       "tipo":"Hidrofobico","dist":round(dd,2),
                                       "ar":ra["name"],"al":la["name"]})
                        seen.add(k)
        # Electrostatico
        if is_charged(ra):
            for la in lig_het:
                dd = d(ra, la)
                if dd <= DIST_EL:
                    k = key_base+("EL",)
                    if k not in seen:
                        result.append({**{kk:ra[kk] for kk in ("chain","resid","resname")},
                                       "tipo":"Electrostatico","dist":round(dd,2),
                                       "ar":ra["name"],"al":la["name"]})
                        seen.add(k)
    return result

def mostrar(ctts, label):
    print(f"\n{'='*62}")
    print(f"  {label}")
    print(f"{'='*62}")
    if not ctts:
        print("  (sin contactos)"); return
    by_tipo = defaultdict(list)
    for c in ctts: by_tipo[c["tipo"]].append(c)
    for tipo in ["H-Bond","Electrostatico","Hidrofobico"]:
        grupo = by_tipo.get(tipo,[])
        if not grupo: continue
        print(f"\n  -- {tipo} ({len(grupo)}) --")
        for c in sorted(grupo, key=lambda x: x["resid"]):
            flag = " [CLAVE]" if c["resid"] in KEY_RES else ""
            print(f"    {c['resname']:3s}{c['resid']:<5d}[{c['chain']}] "
                  f"{c['ar']:4s}<->{c['al']:4s}  {c['dist']:.2f} A{flag}")
    clave = [c for c in ctts if c["resid"] in KEY_RES]
    print(f"\n  Total: {len(ctts)} contactos | Residuos clave: {len(clave)}")

# ── Main ──
rec = parse_pdbqt("01_Receptor/receptor_clean.pdbqt")
print(f"Receptor cargado: {len(rec)} atomos")

resumen = []
for lig in ["fenta","NF1","NF2"]:
    for corrida, suf in [("BLIND","_blind.pdbqt"),("SITE","_site.pdbqt")]:
        path = f"04_Docking/Resultados/{lig}{suf}"
        if not os.path.exists(path):
            print(f"[!] No encontrado: {path}"); continue
        latoms = parse_pdbqt(path, model=1)
        ctts = contacts(rec, latoms)
        mostrar(ctts, f"{lig.upper()} — {corrida}")
        hb = sum(1 for c in ctts if c["tipo"]=="H-Bond")
        hf = sum(1 for c in ctts if c["tipo"]=="Hidrofobico")
        el = sum(1 for c in ctts if c["tipo"]=="Electrostatico")
        ck = sum(1 for c in ctts if c["resid"] in KEY_RES)
        resumen.append((lig,corrida,len(ctts),hb,hf,el,ck))

print(f"\n{'='*70}")
print("  TABLA RESUMEN")
print(f"{'='*70}")
print(f"  {'Ligando':<8}{'Corrida':<8}{'Total':>6}{'H-Bond':>8}{'Hfob':>7}{'Elec':>6}{'Clave':>8}")
print(f"  {'-'*55}")
for r in resumen:
    print(f"  {r[0]:<8}{r[1]:<8}{r[2]:>6}{r[3]:>8}{r[4]:>7}{r[5]:>6}{r[6]:>8}")
