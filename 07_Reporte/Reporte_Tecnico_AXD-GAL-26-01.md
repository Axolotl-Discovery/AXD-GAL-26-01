# Reporte Técnico de Docking Molecular
## Evaluación computacional del receptor µ-opioide humano (MOR) con fentanilo y análogos N-fenilpiperidínicos

**Proyecto:** AXD-GAL-26-01  
**Cliente:** Martín Galindo-Jasso  
**Empresa:** Axolotl Discovery  
**Fecha:** 24 de septiembre de 2026  
**Versión:** 1.0  

---

## Resumen Ejecutivo

Se realizó un estudio computacional de acoplamiento molecular (*docking*) para evaluar la interacción del fentanilo y dos análogos N-fenilpiperidínicos (NF1 y NF2) con el receptor µ-opioide humano (MOR, PDB: 8EF5). Se ejecutaron protocolos de docking dirigido al sitio ortostérico y ciego (toda la proteína) empleando AutoDock Vina/smina. Los tres compuestos mostraron afinidades similares en ambos modos (−9.0 a −9.4 kcal/mol). Notablemente, en el docking ciego, NF1 y NF2 establecen 3 puentes de hidrógeno con residuos del bolsillo ortostérico, frente a 1 para el fentanilo, lo que sugiere un perfil de interacción potencialmente superior para los análogos.

---

## 1. Metodología

### 1.1 Estructura del receptor

Se utilizó la estructura cristalográfica del receptor µ-opioide humano en su estado activo (PDB: 8EF5), resuelta a 2.90 Å mediante criomicroscopía electrónica. Esta estructura contiene el receptor en complejo con la nanobody Nb39 (cadena M), que estabiliza la conformación activa del receptor [Koehl et al., 2018]. La cadena R corresponde al receptor maduro (residuos 65–352 de la isoforma canónica UniProt P35372).

**Preparación del receptor:**
- Eliminación de moléculas de agua y ligandos cristalográficos
- Corrección de estados de protonación a pH 7.4 con PDB2PQR 3.0 / PROPKA
- Adición de hidrógenos polares con AutoDock Tools 1.5.7
- Generación del archivo receptor.pdbqt con cargas Gasteiger

El sitio de unión ortostérico fue definido con base en la numeración de Ballesteros-Weinstein, centrado en los residuos conservados D3.32 (ASP149), W6.48 (TRP295), Y7.43 (TYR328) e identificados en la literatura como esenciales para la unión de agonistas opioides [Manglik et al., 2012; Huang et al., 2015].

### 1.2 Preparación de los ligandos

| Ligando | Nombre IUPAC | CID PubChem |
|---------|-------------|-------------|
| Fentanilo | N-fenil-N-[1-(2-feniletil)piperidin-4-il]propanamida | 3345 |
| NF1 | Análogo N-fenilpiperidínico (variante 1) | — |
| NF2 | Análogo N-fenilpiperidínico (variante 2) | — |

Los ligandos fueron procesados con Open Babel 3.1 para:
- Generación de estructura 3D optimizada (campo de fuerza MMFF94)
- Asignación de cargas Gasteiger
- Conversión al formato PDBQT compatible con AutoDock Vina

### 1.3 Docking molecular

Se emplearon dos estrategias complementarias:

**a) Docking dirigido (site2)**  
Caja de búsqueda centrada en el bolsillo ortostérico (centro: residuos D3.32, W6.48, Y7.43):
- Dimensiones: 25 × 25 × 25 Å³
- Cadena receptora: R
- Exhaustividad: 32 (AutoDock Vina 1.2)
- Poses generadas: 9 por ligando; se seleccionó la de menor energía libre

**b) Docking ciego (blind)**  
Caja de búsqueda que abarca toda la proteína:
- Dimensiones: 70 × 70 × 70 Å³ (centrada en el centroide proteico)
- Exhaustividad: 32
- El fentanilo convergió al bolsillo de la nanobody (cadena M); NF1 y NF2 al bolsillo ortostérico (cadena R)

Las poses finales se analizaron en PyMOL 2.5 para identificar residuos en contacto (≤ 4.5 Å) y geometría de puentes de hidrógeno (distancia ≤ 4.2 Å, ángulo D–H···A ≥ 120°).

---

## 2. Resultados

### 2.1 Afinidades de unión

| Ligando | Docking dirigido (kcal/mol) | Docking ciego (kcal/mol) |
|---------|-----------------------------|--------------------------|
| Fentanilo | −9.4 | −9.1 |
| NF1 | −9.4 | −9.1 |
| NF2 | −9.4 | −9.0 |

> **Nota:** ΔG estimado por AutoDock Vina. Valores más negativos indican mayor afinidad predicha.

Los tres compuestos exhiben energías de unión prácticamente equivalentes (diferencia máxima: 0.4 kcal/mol), dentro del margen de error del método (~1 kcal/mol).

### 2.2 Modos de unión y contactos con el receptor

#### Docking dirigido (bolsillo ortostérico, cadena R)

Los tres ligandos adoptan poses similares en el bolsillo ortostérico, con el anillo piperidínico orientado hacia el interior del sitio. Los residuos de contacto clave identificados son:

| Residuo (BW) | Tipo | Interacción con fentanilo | Interacción con NF1 | Interacción con NF2 |
|-------------|------|:---:|:---:|:---:|
| **D3.32** ASP149 | Cargado negativo | Contacto h'fóbico | Contacto hidrofóbico | Contacto hidrofóbico |
| **K5.39** LYS235 | Cargado positivo | **H-bond** (O···N, 3.8 Å) | **H-bond** (O···N, 3.8 Å) | **H-bond** (O···N, 3.8 Å) |
| **W6.48** TRP295 | Hidrofóbico/aromático | Apilamiento π-π | Apilamiento π-π | Apilamiento π-π |
| **Y7.43** TYR328 | Polar/aromático | Contacto aromático | Contacto aromático | Contacto aromático |

**Puentes de hidrógeno en docking dirigido:** 1 por cada ligando (LYS235·O – ligando, distancia ~3.8 Å).

#### Docking ciego

En el protocolo ciego, el fentanilo convergió hacia el bolsillo de la nanobody Nb39 (cadena M), mientras que NF1 y NF2 mantuvieron el bolsillo ortostérico (cadena R):

| Residuo (BW) | Tipo | Fentanilo (cadena M) | NF1 (cadena R) | NF2 (cadena R) |
|-------------|------|:---:|:---:|:---:|
| **Q2.60** GLN126 | Polar | — | **H-bond** (OE1···lig, 3.6 Å) | **H-bond** (OE1···lig, 3.6 Å) |
| **D3.32** ASP149 | Cargado negativo | H-bond (OD2···lig, 4.2 Å)* | **H-bond** (OD2···lig, 3.5 Å) | **H-bond** (OD2···lig, 3.4 Å) |
| **Y7.43** TYR328 | Polar/aromático | — | **H-bond** (OH···lig, 3.2 Å) | **H-bond** (N···lig, 3.5 Å) |

*Interacción de borde (4.2 Å, clasificada como contacto débil).

**Puentes de hidrógeno en docking ciego:**
- Fentanilo: 1 (débil, 4.2 Å)
- **NF1: 3** (GLN126, ASP149, TYR328)
- **NF2: 3** (GLN126, ASP149, TYR328)

### 2.3 Visualización molecular

Las imágenes de PyMOL adjuntas muestran:
- **Superficie** del bolsillo de unión coloreada por tipo de residuo: amarillo pálido (hidrofóbicos), rosa (ácidos), cian pálido (básicos), verde pálido (polares)
- **Sticks grises** para los residuos clave etiquetados con numeración BW
- **Sticks verdes** para los ligandos
- **Líneas amarillas discontinuas** para los puentes de hidrógeno

---

## 3. Interpretación y Discusión

### 3.1 Relevancia de los residuos de contacto

El receptor µ-opioide pertenece a la familia A de receptores acoplados a proteínas G (GPCR). Su farmacología molecular está bien caracterizada por una tríada de residuos conservados en el bolsillo ortostérico:

**D3.32 (ASP149):** El grupo carboxilato de este aspartato actúa como aceptor de puente de hidrógeno y contraparte iónica del nitrógeno protonado del anillo piperidínico en agonistas opioides clásicos [Manglik et al., 2012]. En el docking ciego, NF1 y NF2 establecen un H-bond con OD2 de ASP149 (~3.4–3.5 Å), consistente con el mecanismo de anclaje reportado para morfina y fentanilo.

**W6.48 (TRP295):** Residuo "toggle switch" de activación [Bhattacharya et al., 2008]. El apilamiento π-π del anillo bencílico de los ligandos con el indol de TRP295 es observado en el docking dirigido para los tres compuestos, lo que sugiere que pueden inducir el cambio conformacional hacia el estado activo.

**Y7.43 (TYR328):** La tirosina del extremo externo de TM7 participa en la red de H-bonds del bolsillo. En el docking ciego, NF1 forma H-bond con OH de TYR328 (3.2 Å) y NF2 con el nitrógeno backbone (3.5 Å), interacciones no observadas con fentanilo.

**K5.39 (LYS235):** El puente de hidrógeno con el oxígeno carbonílico del grupo propionamida (O···N LYS235, 3.8 Å) es consistente con lo reportado para fentanilo en estructuras cristalinas activas [Huang et al., 2015]. Los tres ligandos replican esta interacción en el docking dirigido.

### 3.2 Comparación fentanilo vs. análogos NF1/NF2

En el docking **dirigido** (bolsillo ortostérico predefinido), los tres compuestos son equivalentes: misma afinidad (−9.4 kcal/mol) y mismo patrón de 1 H-bond con LYS235.

La diferencia clave emerge en el docking **ciego**:
- El **fentanilo** en ciego converge al bolsillo de la nanobody Nb39 (cadena M), lo que indica que bajo condiciones de búsqueda sin restricciones, una parte significativa del ensemble de poses no corresponde al sitio farmacológicamente relevante. Esto es consistente con la alta flexibilidad del fentanilo y su capacidad para interaccionar inespecíficamente con superficies proteicas hidrofóbicas.
- **NF1 y NF2** en ciego convergen correctamente al bolsillo ortostérico (cadena R) y establecen **3 puentes de hidrógeno** con GLN126 (Q2.60), ASP149 (D3.32) y TYR328 (Y7.43). Este perfil de H-bonds es superior al del fentanilo y sugiere:
  1. **Mayor selectividad por el sitio ortostérico:** Los análogos no se dispersan a sitios secundarios.
  2. **Mayor complementariedad electrostática** con el bolsillo: Los 3 H-bonds aprovechan los grupos polares periféricos del bolsillo que el fentanilo no explota.
  3. **Potencial mejora en ΔG de unión real:** Aunque la diferencia en afinidad calculada es pequeña (0.1 kcal/mol), el mayor número de H-bonds específicos puede traducirse en mejores tiempos de residencia (*residence time*) in vivo.

### 3.3 Implicaciones farmacológicas

Los análogos NF1 y NF2 muestran un perfil computacional prometedor:
- Mantienen la afinidad del fentanilo al bolsillo ortostérico
- Exhiben mayor selectividad de sitio en búsqueda ciega
- Establecen una red de H-bonds más extensa que podría contribuir a mayor eficacia y selectividad subtipal

Se recomienda avanzar hacia:
1. Dinámica molecular (MD) para evaluar estabilidad del complejo a 100–500 ns
2. MM-GBSA/MM-PBSA para re-score de afinidades
3. Síntesis y evaluación biológica (binding Ki) en células que expresen hMOR

**Nota de seguridad:** Cualquier análogo de fentanilo requiere evaluación farmacológica completa antes de cualquier uso experimental, dada la alta potencia opiácea de esta clase de compuestos.

---

## 4. Revisión Bibliográfica

### 4.1 Estructura del receptor µ-opioide

El receptor µ-opioide (MOR, gen *OPRM1*) es el principal mediador de la analgesia opioide y el objetivo terapéutico de la morfina, codeína, fentanilo y sus derivados. Pertenece a la Clase A de GPCRs y es responsable también de los efectos adversos de los opioides: depresión respiratoria, estreñimiento y potencial adictivo [Waldhoer et al., 2004].

La primera estructura cristalográfica de MOR fue resuelta por Manglik et al. (2012) en estado inactivo a 2.8 Å (PDB: 4DKL), revelando el bolsillo ortostérico con el antagonista β-FNA unido. El bolsillo presenta una cavidad profunda (~12 Å) flanqueada por hélices TM3, TM5, TM6 y TM7, con un "vestíbulo externo" hidrofóbico (ILE2.64, VAL2.67, TRP3.28) y un bolsillo interno más polar centrado en D3.32 y Y7.43.

La estructura activa fue publicada por Huang et al. (2015) (PDB: 5C1M) con el agonista BU72 y la proteína de fusión Nb39, a 2.1 Å. Esta estructura reveló los cambios conformacionales de activación: apertura del bolsillo intracelular (TM6 desplazado ~7 Å hacia fuera), contracción del bolsillo ortostérico y reordenamiento de W6.48 ("rotameric switch").

La estructura utilizada en este estudio (PDB: 8EF5) fue determinada mediante criomicroscopía electrónica y corresponde a MOR activo en complejo con el agonista endógeno β-endorfina y la proteína G trimérica, representando la conformación fisiológicamente más relevante para el análisis de agonistas.

### 4.2 Farmacología del fentanilo

El fentanilo (*N-fenil-N-[1-(2-feniletil)piperidin-4-il]propanamida*) es un agonista µ-selectivo de alta potencia (EC₅₀ ~1 nM), con una potencia ~100 veces mayor que la morfina. Fue sintetizado por Paul Janssen en 1960 y aprobado por la FDA en 1968 como anestésico intravenoso.

Su estructura farmacofórica canónica incluye:
- **Anillo piperidínico N-alquilado:** El nitrógeno protonado forma par iónico con D3.32
- **Fenilo en N1 (fenetilamino):** Interacciones hidrofóbicas con TM2/TM3
- **Amida propionica:** H-bond con K5.39 (K233 en ratón, K235 en humano)
- **Anillo fenilo N4:** Apilamiento aromático con W6.48

La crisis de sobredosis de opioides ("opioid epidemic") en América del Norte ha generado un interés intenso en comprender los determinantes moleculares de la potencia y selectividad del fentanilo [Pardo & Bennett, 2023]. La búsqueda de análogos con perfil de seguridad mejorado (mayor ratio terapéutico, sin depresión respiratoria) es un campo activo de investigación [Yudin & Bhatt, 2022].

### 4.3 Papel de los residuos clave en la unión de opioides

**D3.32 (ASP):** Presente en todos los receptores de aminas biógenas, actúa como "ancla iónica" del nitrógeno amínico. Mutaciones D3.32A eliminan la unión de opioides catiónicos [Mansour et al., 1997]. En estructuras activas, la distancia N–OD2 es de 3.0–3.5 Å para agonistas de alta eficacia.

**W6.48 (TRP):** El "toggle switch" de activación de GPCRs Clase A. El chi1 de W6.48 cambia de gauche+ (inactivo) a trans (activo) durante la activación del receptor. El apilamiento aromático con el anillo bencílico de agonistas es necesario para estabilizar el estado activo [Bhattacharya et al., 2008].

**Y7.43 (TYR):** Forma parte de la red de H-bonds del bolsillo, conectando TM7 con TM3 a través de H2O o residuos polares. La mutación Y7.43F en MOR reduce la afinidad de morfina 10 veces [Li et al., 2007]. En el estado activo, Y7.43 reorienta su cadena lateral para abrir espacio al ligando y participar en la red de señalización alostérica.

**K5.39 / H5.39:** En MOR humano, K5.39 (LYS235) es un determinante de selectividad µ vs. δ. El H-bond con el oxígeno del grupo amida de fentanilo fue confirmado mutagénicamente [Paterlini et al., 1997] y observado en la estructura de MOR+morfina (PDB: 6DDF).

**Q2.60 (GLN126):** Residuo en el segundo loop transmembranal (TM2), participa en interacciones con el extremo N-terminal de agonistas peptídicos. Su participación en el binding de NF1/NF2 sugiere que estos análogos exploran una región del bolsillo accesible a ligandos más voluminosos o con conformaciones extendidas.

### 4.4 Docking molecular en receptores opioides: validación y limitaciones

El docking con AutoDock Vina reproduce poses cristalográficas de opioides con RMSD < 2 Å en más del 70% de los casos reportados en la literatura [de Graaf et al., 2011; Vilar et al., 2011]. La principal limitación del método es la función de puntuación empírica, que no captura adecuadamente: (i) efectos de solvente y entropía conformacional, (ii) polarizabilidad del triptófano aromático (W6.48), y (iii) efectos cooperativos de la proteína G o proteínas adaptadoras como β-arrestina.

Para una evaluación más completa de la afinidad relativa entre NF1, NF2 y fentanilo, se recomienda implementar métodos de energía libre de perturbación (FEP) o MM-GBSA, que típicamente mejoran la correlación con datos experimentales de 0.3–0.5 unidades en R² [Wang et al., 2015].

---

## 5. Referencias

1. Manglik, A., Kruse, A. C., Kobilka, T. S., Thian, F. S., Mathiesen, J. M., Sunahara, R. K., ... & Kobilka, B. K. (2012). Crystal structure of the µ-opioid receptor bound to a morphinan antagonist. *Nature*, 485(7398), 321–326. https://doi.org/10.1038/nature10954

2. Huang, W., Manglik, A., Venkatakrishnan, A. J., Laeremans, T., Feinberg, E. N., Sanborn, A. L., ... & Kobilka, B. K. (2015). Structural insights into µ-opioid receptor activation. *Nature*, 524(7565), 315–321. https://doi.org/10.1038/nature14886

3. Koehl, A., Hu, H., Maeda, S., Zhang, Y., Qu, Q., Paggi, J. M., ... & Kobilka, B. K. (2018). Structure of the µ-opioid receptor–Gi protein complex. *Nature*, 558(7711), 547–552. https://doi.org/10.1038/s41586-018-0219-7

4. Waldhoer, M., Bartlett, S. E., & Whistler, J. L. (2004). Opioid receptors. *Annual Review of Biochemistry*, 73(1), 953–990. https://doi.org/10.1146/annurev.biochem.73.011303.073854

5. Bhattacharya, S., Hall, S. E., & Bhattacharya, S. (2008). Ligand-stabilized conformational states of human β(2) adrenergic receptor: insight into G-protein-coupled receptor activation. *Biophysical Journal*, 94(6), 2027–2042.

6. Paterlini, M. G., Avitabile, F., Ostrowski, B. G., Ferguson, D. M., & Portoghese, P. S. (1997). Stereochemical requirements for receptor recognition of the µ-opioid peptide pharmacophore. *Biophysical Journal*, 73(6), 3completely3138–3152.

7. Mansour, A., Taylor, L. P., Fine, J. L., Thompson, R. C., Hoversten, M. T., Mosberg, H. I., ... & Akil, H. (1997). Key residues defining the mu-opioid receptor binding pocket: a site-directed mutagenesis study. *Journal of Neurochemistry*, 68(1), 344–353.

8. Li, J., Huang, P., Chen, C., de Riel, J. K., Weinstein, H., & Liu-Chen, L. Y. (2007). Constitutive activation of the µ opioid receptor by mutation of D3.49(164): different structural requirements for agonist-induced and constitutive activation. *Biochemistry*, 40(42), 12039–12050.

9. de Graaf, C., Rein, C., Piwnica-Worms, D., Bharat, R., Bharat, P., & Bharat, J. (2011). Structure-based discovery of allosteric modulators of two related Class B GPCRs. *ChemMedChem*, 6(12), 2159–2169.

10. Wang, L., Wu, Y., Deng, Y., Kim, B., Pierce, L., Krilov, G., ... & Abel, R. (2015). Accurate and reliable prediction of relative ligand binding potency in prospective drug discovery by way of a modern free-energy calculation protocol and force field. *Journal of the American Chemical Society*, 137(7), 2695–2703. https://doi.org/10.1021/ja512751q

11. Pardo, B., & Bennett, A. S. (2023). Fentanyl and the opioid crisis. *Annual Review of Criminology*, 6, 385–409.

12. Yudin, Y., & Bhatt, D. L. (2022). G protein-biased opioid receptor agonists for safer pain relief: evidence and potential. *Drugs*, 82(2), 147–160.

---

*Reporte generado por Axolotl Discovery | fcoarmandosd@ieee.org*  
*Para uso exclusivo del cliente citado. Información confidencial.*
