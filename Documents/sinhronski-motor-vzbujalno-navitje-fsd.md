# Sinhronski motor z vzbujalnim navitjem — načrtovanje in optimizacija (FSD)

## 1. Pregled sistema

### 1.1 Namen
Razvoj programskega orodja (Python) za analitično dimenzioniranje, večkriterijsko optimizacijo (z genetskimi algoritmi) in numerično verifikacijo (FEMM) 6-polnega 3-faznega sinhronskega elektromotorja z vzbujalnim navitjem na rotorju. Cilj je seminarska naloga pri predmetu *Konstruiranje električnih strojev* (FE UL).

### 1.2 Problemska izjava
Klasično ročno dimenzioniranje sinhronskega stroja z vzbujalnim navitjem je dolgotrajno in vrne le eno rešitev. Naloga zahteva sistematično preiskovanje konstrukcijskega prostora ter izbiro **Pareto optimalnih** rešitev glede na izkoristek in volumen aktivnega dela. Vsako od izbranih rešitev je treba dodatno verificirati v 2D FEM simulaciji (FEMM).

### 1.3 Uporabniki / deležniki
- **Avtor naloge** (študent): poganja program, interpretira rezultate, piše poročilo.
- **Mentor / nosilec predmeta**: ocenjuje delovanje, ustreznost analitičnega izračuna in FEMM modela.
- **Zagovorno komisija**: pričakuje, da študent zna pojasniti vsako enačbo, vsako optimizacijsko spremenljivko in vsako sliko v poročilu.

### 1.4 Cilji
- **C1**: Analitično izračunati popolno geometrijo 6-polnega stroja iz vhodnih podatkov (200 V faz., 50 kW @ 7000 vrt/min, n_max = 14000 vrt/min, material M330-35A).
- **C2**: Izvesti dvokriterijsko optimizacijo (max η, min V) z genetskim algoritmom in zgraditi Pareto fronto.
- **C3**: Izbrati 5 Pareto optimalnih rešitev in vsako parametrično izrisati v FEMM s pravilno geometrijo statorskih in rotorskih zob (rotorski zob z **zaokroženim sinusnim koncem**, kot na sliki `stroj.JPG`; statorski zob razširjen ob zračni reži, kot na sliki `utor.JPG`).
- **C4**: Iz FEMM simulacij izvesti: (i) prosti tek – izračun inducirane napetosti pri 7000 vrt/min, (ii) statorsko napajanje – navorna karakteristika in Fourierova analiza.
- **C5**: Generirati podrobno pisno poročilo z vsemi enačbami, tabelo spremenljivk, razlago GA, FEMM postopka in komentarji vseh 5 izbranih rešitev.

### 1.5 Necilji
- Termični in mehanski izračun (npr. centrifugalne napetosti rotorja) ni del naloge.
- 3D simulacija (npr. čelno polje navitij) ni del naloge.
- Generiranje poročila s strani programa **ni** zahtevano — poročilo piše študent ročno.
- Avtomatska izdelava izvedbenih risb (CAD) ni del naloge.

### 1.6 Tok visokonivojskega delovanja
```
[Vhodni podatki] → [Analitični izračun (modul GEOM)] → [GA optimizacija (modul OPT)]
   → [Pareto fronta] → [Izbor 5 rešitev] → [Parametrični izris v FEMM (modul DRAW)]
   → [FEMM analiza: prosti tek + navor (modul SIM)] → [Post-procesiranje (FFT, U_ind)]
   → [Izvoz rezultatov + slike za poročilo]
```

---

## 2. Arhitektura sistema

### 2.1 Logična arhitektura

Sistem je sestavljen iz **šestih ločenih modulov** (paketov), ki se izvajajo zaporedno v eni cevni liniji (pipeline). Vmesni rezultati se hranijo v slovarjih/`numpy` poljih in po želji serializirajo v `.npz`/`.json` za ponovljivost.

| Modul | Odgovornost | Vhod | Izhod |
|-------|-------------|------|-------|
| `inputs` | Privzete in CLI vhodne spremenljivke | argparse / `inputs.yaml` | `dict` parametrov |
| `analytical` | Analitični preračun geometrije, izgub in η | parametri stroja + GA gen | `MotorDesign` objekt |
| `optimization` | Genetski algoritem (NSGA-II), generiranje populacije, fitness | meje optimizacijskih spr. | populacija + Pareto fronta |
| `pareto` | Izbira 5 reprezentativnih rešitev s Pareto fronte | Pareto fronta | seznam 5 `MotorDesign` |
| `femm_model` | Parametrični izris geometrije v FEMM (`pyFEMM`) | `MotorDesign` | `.fem` datoteka |
| `femm_sim` | Vodenje simulacij prostega teka in navora | `.fem` + scenarij | navor(θ), Ψ(θ), B_max |
| `post` | FFT navora, U_ind, izgube Fe, izračun η iz FEMM | rezultati simulacij | grafi `.png`, CSV |

### 2.2 Programska / izvajalna platforma

- **Programski jezik**: Python ≥ 3.10 (priporočeno 3.11; primer1 dokazano deluje).
- **Operacijski sistem**: Windows (FEMM 4.2 je primarno Windows aplikacija). Pod macOS/Linux preko Wine (kot v `primer1/python/KES.py`).
- **Zunanji programi**:
  - **FEMM 4.2** — Finite Element Method Magnetics (David Meeker).
- **Ključne knjižnice**:
  - `pyFEMM` (`pip install pyfemm`) — Python ovojnica za FEMM ALC ukaze.
  - `numpy`, `scipy` — vektorizacija, FFT, polifit, optimizacija.
  - `pymoo` ali `DEAP` — NSGA-II genetski algoritem *(assumed: `pymoo`, ker je trenutno najbolj zrelo orodje za večkriterijsko GA)*.
  - `matplotlib` — grafi (Pareto fronta, U_ind, navor, FFT).
  - `pandas` — tabele rezultatov, izvoz v CSV.
  - `pyyaml` *(assumed)* — branje konfiguracije.
- **B-H krivulja**: shranjena v `BH.txt` (kot v primer1, format `B H` po vrsticah), naložena v `analytical` in `femm_model` modula.
- **Železove izgube**: polinom 2D (2,3) `f(B, f)` zgrajen iz tabele v `izgube_fem.txt`. Implementacija prevzeta iz MATLAB skripte in pretvorjena v Python (`numpy.polyfit` ali `scipy.optimize.curve_fit`).

### 2.3 Programska arhitektura — strukturni pregled

```
KES_v7/
├── src/
│   ├── inputs.py              # privzeti parametri (200 V, 50 kW, 7000 rpm…)
│   ├── analytical.py          # MotorDesign + analitični izračun
│   ├── losses.py              # polinom železovih izgub, izračun P_Fe, P_Cu
│   ├── optimization.py        # NSGA-II setup, fitness, omejitve
│   ├── pareto.py              # izbira 5 reprezentantov
│   ├── femm_model.py          # narisi_geometrijo() + koncaj_geometrijo()
│   ├── femm_sim.py            # prosti_tek(), navor() — vrtenje + integrali
│   ├── post.py                # FFT, U_ind, grafi
│   ├── plot.py                # matplotlib pomožne funkcije
│   └── main.py                # orkestracija (CLI)
├── BH.txt                     # B-H krivulja M330-35A
├── data/izgube_M330-35A.csv   # tabela P_Fe(B,f)
├── outputs/
│   ├── pareto.png
│   ├── fem/<design_id>.fem
│   ├── fem/<design_id>.ans
│   ├── results.csv
│   └── figures/<design_id>/{uind.png, navor.png, fft.png}
└── Documents/
    ├── sinhronski-motor-vzbujalno-navitje-fsd.md   # ta dokument
    └── porocilo.md / porocilo.tex                  # končno poročilo (ročno)
```

### 2.4 Zagonski tok
1. `main.py` prebere vhodne podatke.
2. Inicializira GA populacijo (`optimization.py`).
3. Za vsak osebek populacije pokliče `analytical.evaluate()` (vrne η, V, P_Fe, P_Cu, izvedljivost).
4. NSGA-II konvergira; vrne Pareto fronto.
5. `pareto.select_five()` izbere 5 reprezentantov (npr. ekstremi + 3 vmesni).
6. Za vsak reprezentant: `femm_model.build(design)` → `.fem` datoteka.
7. Za vsak `.fem`: `femm_sim.no_load()` in `femm_sim.torque()`.
8. `post.analyze()` izračuna grafe, FFT, primerjavo analitično vs. FEMM.
9. Izvozi `results.csv` in figure.

---

## 3. Implementacijske faze

### 3.1 Faza 1 — Analitični izračun in temelji
**Obseg**: Realizacija analitičnega modela enega stroja.
**Rezultati**:
- `inputs.py` z vsemi privzetimi parametri.
- `analytical.py`: razred `MotorDesign` z metodami `compute_geometry()`, `compute_emf()`, `compute_losses()`, `efficiency()`.
- `losses.py`: polinomska aproksimacija izgub (B, f) iz `izgube_fem.txt`.
- B-H krivulja naložena iz `BH.txt`.
- Unit testi za analitične enačbe (preverjanje konsistence z `PRIMER1/porocilo.pdf`).

**Izhodna merila**:
- Za vhodne podatke 200 V / 50 kW / 7000 rpm program v < 1 s vrne popolno geometrijo (D_r, l, b_ds, h_ds, h_ys, Q_s, N_s, q, b_dr, h_yr, vzbujalni tok I_m, vzbujalni ampere-ovoji N_r·I_m).
- Izračunani izkoristek je v razumnem območju (η > 0.9 za ciljno točko).
- Geometrija je izvedljiva: vsa navitja se prilegajo utorom (K_Cu ≤ 0.4 za stator, K_Cu_rot ≤ 0.7 za rotor *(assumed)*).

**Odvisnosti**: izhodišče (brez prejšnjih faz).

### 3.2 Faza 2 — Optimizacija in Pareto fronta
**Obseg**: Večkriterijska optimizacija stroja.
**Rezultati**:
- `optimization.py`: NSGA-II konfiguracija, dekoder gen → `MotorDesign`, fitness funkcija (η in V), trde omejitve (J_Cu_s ≤ 10 A/mm², J_Cu_r ≤ 5 A/mm², U_ind ≤ 0.95·U_grid).
- `pareto.py`: izbor 5 točk z Pareto fronte.
- Graf Pareto fronte (`outputs/pareto.png`).

**Izhodna merila**:
- GA konvergira v < 10 min na ciljnem PC (populacija ≥ 100, generacije ≥ 50) *(assumed: parametri primerni za nalogo te velikosti — prilagodljivo)*.
- Pareto fronta vsebuje vsaj 20 nedominiranih rešitev.
- 5 izbranih rešitev pokriva: 1× max η, 1× min V, 3× vmesne (enakomerno razporejene po krivulji).
- Vse trde omejitve so izpolnjene v vseh 5 rešitvah.

**Odvisnosti**: Faza 1.

### 3.3 Faza 3 — FEMM model in simulacije
**Obseg**: Parametrični izris in numerična verifikacija.
**Rezultati**:
- `femm_model.py`: funkciji `narisi_geometrijo(design)` in `koncaj_geometrijo(design)` po vzoru `primer1/python/KES_geometrija.py`.
  - Rotorski zob: točke generirane iz `zr_a = zr_min / cos(α·κ)` (sinusna porazdelitev zračne reže), zaokrožen rob.
  - Statorski zob: razširjen ob zračni reži (cevelj), kot v `utor.JPG`; klic `mi_createradius`.
  - Vzbujalno navitje: 2 bloka `Copper` na vsakem rotorskem zobu (en pozitiven, en negativen), z lastnostmi NVS DC tokokroga.
  - Statorsko navitje: razporeditev `["A","c","c","B","B","a","a","C","C","b","b","A", …]` za q = 2 (Q_s = 36); za q ≠ 2 generirano dinamično iz EMETOR korelacij.
  - Materiali: `M330-35A` z B-H krivuljo iz `BH.txt`, `Copper`, `Air`.
- `femm_sim.py`:
  - **Prosti tek** (točka L): I_stator = 0, vzbujanje = I_m, rotor vrteti za 360°/p (= 120° za p = 3) v inkrementih (npr. 5°), izračun Ψ_A(θ) → U_ind(θ) = dΨ/dt pri n = 7000 rpm.
  - **Navor** (točka M): I_A = √2·I_n, I_B = -I_A/2, I_C = -I_A/2, vzbujanje = I_m, rotor vrteti za 120°, integral `mo_blockintegral(22)` za vsako pozicijo.
- Izhodne `.fem` in `.ans` datoteke shranjene v `outputs/fem/`.

**Izhodna merila**:
- Za vsako od 5 rešitev se model uspešno izriše brez "node overlap" napak v FEMM.
- Mreženje (`smartmesh(0)`) konča v < 30 s na rešitev.
- Inducirana napetost U_ind ≤ 0.95 · U_grid = 0.95 · 200 V·√2 ≈ 269 V (amplituda) oz. 190 V (RMS) za vse 5 rešitev.
- Geometrija vizualno ustreza slikam `stroj.JPG` in `utor.JPG`.

**Odvisnosti**: Faza 2.

### 3.4 Faza 4 — Post-procesiranje in poročilo
**Obseg**: Analiza rezultatov in priprava gradiva za poročilo.
**Rezultati**:
- `post.py`:
  - Izračun U_ind iz Ψ(θ): U_ind = dΨ/dt; preračun na 7000 rpm.
  - FFT navorne krivulje (`numpy.fft`): osnovna harmonska + 2. harmonska + THD.
  - Izračun B_max v statorskem zobu in jarmu (`mo_getpointvalues` v zanki po Q_s).
  - Izračun izgub v železu iz B_max in nazivne frekvence (f = n·p/60 = 7000/60·3 = 350 Hz).
  - Primerjava analitično vs. FEMM (η, U_ind, B_max) v tabeli.
- Grafi (`figures/<design_id>/`):
  - U_ind(θ), spekter U_ind (FFT).
  - Navor(θ), spekter navora (FFT, 1. in 2. harmonska).
  - Skupna Pareto fronta z označenimi 5 rešitvami.
- `results.csv` z geometrijo, izgubami in η za vseh 5 rešitev.

**Izhodna merila**:
- Vse slike so v dovolj visoki ločljivosti (≥ 300 dpi) za vključitev v poročilo.
- THD navora < 15 % za vse 5 rešitev *(assumed cilj — ni eksplicitno zahtevan v opisu)*.
- Tabela primerjave analitično–FEMM pokaže relativno odstopanje < 10 % za U_ind in η.

**Odvisnosti**: Faza 3.

---

## 4. Zahteve

### 4.1 Funkcionalne zahteve (FR)

#### Analitični izračun (točke A–F naloge)
- **FR-1.1** [Must]: Program **shall** iz vhodnih podatkov (m = 3, p = 3, U_f = 200 V, P_c = 50 kW, n_c = 7000 vrt/min, n_max = 14000 vrt/min, material M330-35A) analitično izračunati popolno geometrijo stroja: D_r, l, zračna reža δ, Q_s, q, N_s, Z_q, b_ds, h_ds, h_ys, b_dr, h_yr, b_1s, b_1r.
- **FR-1.2** [Must]: Program **shall** določiti q ∈ {1, 1.5, 2, 2.5, 3} in odgovarjajoč k_w1 iz EMETOR korelacij; q sme biti optimizacijska spremenljivka.
- **FR-1.3** [Must]: Program **shall** kot konstrukcijsko omejitev upoštevati max gostoto toka J_Cu_s ≤ 10 A/mm² in J_Cu_r ≤ 5 A/mm².
- **FR-1.4** [Must]: Program **shall** določiti funkcijo B_δ(I_m, θ_m) — amplitudna vrednost gostote magnetnega pretoka v zračni reži v odvisnosti od vzbujalnega toka in vzbujalne magnetne napetosti.
- **FR-1.5** [Must]: Program **shall** s polinomsko aproksimacijo `poly23` določiti izgube v železu P_Fe(B, f) za material M330-35A iz priložene tabele.
- **FR-1.6** [Must]: Program **shall** izračunati izgube v železu ločeno v statorskih zobeh in jarmu, z empiričnimi koeficienti k_Fe_zob = 2 in k_Fe_jarem = 1.6 *(prevzeto iz `izgube_fem.txt`)*.
- **FR-1.7** [Must]: Program **shall** izračunati izgube v bakru P_Cu_s in P_Cu_r ob delovni temperaturi T_del = 80 °C, s temperaturnim koeficientom α_Cu = 0.00381 K⁻¹ in σ_Cu = 34·10⁶ S/m; dolžina čel navitja **shall** biti privzeto enaka aktivni dolžini paketa (L_r).
- **FR-1.8** [Must]: Program **shall** v vogalni točki (P_c, n_c) izračunati izkoristek η = P_c / (P_c + P_Fe + P_Cu).

#### Optimizacija (točke G, H, K naloge)
- **FR-2.1** [Must]: Program **shall** ovrednotiti analitični model znotraj GA z dvokriterijsko funkcijo `[minimize(1/η), minimize(V)]`, kjer V = π · (D_stator/2)² · l.
- **FR-2.2** [Must]: Program **shall** generirati Pareto fronto v 2D prostoru (η vs. V) z NSGA-II.
- **FR-2.3** [Must]: Program **shall** kot trde omejitve uporabiti: J_Cu_s ≤ 10, J_Cu_r ≤ 5, U_ind_analitično ≤ 0.95 · U_grid, izvedljivost navitij (zapolnitev utora ≤ K_Cu_max), izogibanje preploščatim geometrijam (l/D_r ≥ 0.4 *(assumed mejnik)*).
- **FR-2.4** [Must]: Program **shall** iz Pareto fronte izbrati 5 reprezentativnih rešitev.
- **FR-2.5** [Should]: Program **shall** izrisati Pareto fronto z označenimi 5 izbranimi rešitvami v PNG.

#### FEMM izris in simulacije (točke I, J, L, M naloge)
- **FR-3.1** [Must]: Program **shall** za vsako od 5 rešitev parametrično zgraditi 2D FEMM model (`.fem`) preko `pyFEMM`.
- **FR-3.2** [Must]: Rotorski zobje **shall** imeti sinusno aproksimirano zračno režo (zr(α) = zr_min / cos(κ·α)) z zaokroženim koncem, tako da je v sredini pola najmanjša reža.
- **FR-3.3** [Must]: Statorski zobje **shall** imeti razširjeno glavo ob zračni reži (cevelj), tako da se delno raztezajo v sosednje utorne odprtine; uporabi se `mi_createradius`.
- **FR-3.4** [Must]: Vzbujalno navitje **shall** biti modelirano kot dva ločena bloka po rotorskem zobu (en blok ovojev v pozitivni, en v negativni smeri), po en par na vsakem od 6 polov.
- **FR-3.5** [Must]: Statorsko navitje **shall** biti razporejeno v 36 utorih (q = 2) z razporeditvijo `A c c B B a a C C b b A …`; za drugi q se razporeditev generira iz EMETOR korelacij.
- **FR-3.6** [Must]: Materiali **shall** biti definirani: `M330-35A` z B-H krivuljo iz `BH.txt`, `Copper`, `Air`. Robni pogoj `A=0` na zunanjem statorskem robu.
- **FR-3.7** [Must]: Program **shall** za vsako rešitev izvesti **prosti tek**: I_A = I_B = I_C = 0, vzbujanje = I_m_analitično, vrtenje rotorja za 360°/p s korakom ≤ 5°, izračun Ψ_A(θ).
- **FR-3.8** [Must]: Program **shall** izračunati amplitudo U_ind = max|dΨ_A/dt| · (n_c/60)/(360°/p) za vsako rešitev in preveriti pogoj U_ind ≤ 0.95·U_grid.
- **FR-3.9** [Must]: Program **shall** za vsako rešitev izvesti **simulacijo navora**: I_A = √2·I_n, I_B = -I_A/2, I_C = -I_A/2, vzbujanje = I_m, vrtenje za 360°/p s korakom ≤ 5°, integral `mo_blockintegral(22)` za vsak korak.
- **FR-3.10** [Must]: Program **shall** s Fourierovo transformacijo razgraditi navorno krivuljo na osnovno in višje harmonske komponente; izrecno izpisati 2. harmonsko.

#### Izvoz rezultatov (povezuje s točko N in poročilom)
- **FR-4.1** [Must]: Program **shall** za vsako od 5 rešitev shraniti: `.fem`, `.ans`, sliko geometrije, U_ind(θ) graf, navor(θ) graf, FFT spekter navora.
- **FR-4.2** [Must]: Program **shall** v `results.csv` izvoziti za vsako rešitev: geometrijo, P_Fe, P_Cu, η, U_ind, M_avg, M_ripple, B_max_zob, B_max_jarem.
- **FR-4.3** [Should]: Program **shall** v ločeni datoteki natisniti tabelo vseh spremenljivk: kratico + polno ime + enoto + vir (za poročilo).

#### Komentarji in koda
- **FR-5.1** [Must]: Vsaka funkcija v kodi **shall** vsebovati komentarje, ki pojasnijo fizikalno ozadje, in **shall** vsebovati referenco na vir (npr. *Pyrhönen, Design of rotating electrical machines, 2nd ed., enačba X.Y*).
- **FR-5.2** [Should]: Koda **shall** biti pregledna in modulirana — ena odgovornost na modul/funkcijo.

### 4.2 Nefunkcionalne zahteve (NFR)

- **NFR-1.1** [Must, učinkovitost]: Eno ovrednotenje analitičnega modela (potreba GA) **shall** trajati < 100 ms na sodobnem PC.
- **NFR-1.2** [Should, učinkovitost]: Celoten GA zagon (≥ 100 osebkov × 50 generacij) **shall** se zaključiti v < 15 min na sodobnem PC.
- **NFR-1.3** [Should, učinkovitost]: Ena FEMM simulacija (prosti tek ali navor, 24 pozicij rotorja) **shall** se zaključiti v < 5 min.
- **NFR-2.1** [Must, ponovljivost]: Z istimi vhodnimi parametri in istim seed za GA **shall** program vrniti enako Pareto fronto (deterministično).
- **NFR-2.2** [Must, ponovljivost]: Vmesni rezultati (Pareto fronta, izbranih 5 rešitev) **shall** biti serializirani v `outputs/` za ponovno uporabo brez ponovne optimizacije.
- **NFR-3.1** [Must, pravilnost]: Analitično izračunan U_ind **shall** odstopati od FEMM rezultata za < 10 % za izbranih 5 rešitev.
- **NFR-3.2** [Must, pravilnost]: Analitično izračunan η **shall** odstopati od η izračunanega iz FEMM (P_Fe + P_Cu) za < 10 %.
- **NFR-4.1** [Should, berljivost]: Vse javne funkcije **shall** imeti docstring (Slovensko ali Angleško, dosledno).
- **NFR-4.2** [Should, prenosljivost]: Koda **shall** delovati na Windows + native FEMM; macOS/Linux + Wine je *nice-to-have* (kot v `primer1/python/KES.py`).
- **NFR-5.1** [Must, sledljivost virov]: Vsaka enačba v kodi in poročilu **shall** imeti citat (knjiga + poglavje/enačba ali znanstveni članek).

### 4.3 Omejitve

- **Tehnologija**: FEMM 4.2 podpira le 2D (planarno ali aksisimetrično) analizo; 3D efekti (čela navitij) so v analitičnem modelu zajeti preko empiričnih koeficientov, ne pa direktno simulirani.
- **Material**: B-H krivulja in tabela izgub M330-35A sta fiksni (priložen `BH.txt` in `izgube_fem.txt`); drugih materialov v okviru te naloge **ne** preverjamo.
- **Številčne omejitve**: q ∈ {1, 1.5, 2, 2.5, 3} (diskretno) — vpliva na število utorov Q_s = 2·p·m·q.
- **Geometrijske omejitve**: l/D_r ≥ 0.4 (izogibanje preploščatim strojem) *(assumed iz teksta naloge)*.
- **Predmetne omejitve**: Naloga zahteva, da je dolžina čela navitij enaka L_r (FR-1.7) — to je poenostavitev, ki dvigne P_Cu vrednost in pomaga z optimizacijo robov.

---

## 5. Tveganja, predpostavke in odvisnosti

### 5.1 Tehnična tveganja

| ID | Tveganje | Verjetnost | Vpliv | Ublažitev |
|----|----------|-----------|-------|-----------|
| R-1 | FEMM model ne mreži zaradi prekrivanja vozlišč ali napačne sinusne krivulje rotorja | Sredna | Visok | Robusten kod za risanje s preverjanjem razdalj med vozlišči (kot v `primer1`), `smartmesh(0)`, ročna preverba geometrije na enem primeru pred poganjanjem GA. |
| R-2 | GA ne konvergira ali Pareto fronta je prazna zaradi pretrdih omejitev | Sredna | Visok | Najprej zagon brez U_ind omejitve, nato dodajanje omejitev iterativno; ohlapne začetne meje optimizacijskih spremenljivk. |
| R-3 | Analitično U_ind močno odstopa od FEMM zaradi neupoštevanja nasičenja | Visoka | Sredna | Uporaba Carterjevega koeficienta in faktorja nasičenja k_sat (≈ 1.05–1.2) iz Pyrhönen knjige; primerjava na enem primeru pred polnim zagonom. |
| R-4 | pyFEMM klic se obesi (npr. FEMM odpre dialog) ali pade na specifični geometriji | Sredna | Sredna | Try/except okrog vsakega `mi_*` klica, logiranje v datoteko, samodejni preskok problematične rešitve in nadaljevanje. |
| R-5 | FEMM neuspešno najde rešitev pri ekstremnih GA gen (npr. zelo tanek zob) | Visoka | Sredna | Trde omejitve v fitness funkciji + kazenska funkcija (penalty) za nedopustne osebke; FEMM klicati le na končnih 5 rešitvah, ne za vsak GA osebek. |
| R-6 | Polinom `poly23` izgub se ne ujema dobro z visoko-frekvenčnim območjem (npr. 350 Hz, blizu zg. meje tabele 400 Hz) | Nizka | Sredna | Preveriti residuale; po potrebi uporabiti 2D interpolacijo (`scipy.interpolate.RBFInterpolator`) namesto polinoma. |
| R-7 | Avtor ne razume vseh izračunov (težave pri zagovoru) | Visoka | Visok | Vsi koraki v kodi z izčrpnimi komentarji + sklici na vir; ročno preverjanje rezultatov za en primer "na papirju". |

### 5.2 Predpostavke (vse označene `(assumed)`)

- A1: Uporabljen bo Python (ne MATLAB) — usklajeno z navodilom *"Priporočam python"*. (Tudi `primer1` je Python.)
- A2: GA knjižnica = `pymoo` (NSGA-II). MATLAB Toolbox `gamultiobj` ni na voljo.
- A3: K_Cu (faktor zapolnitve utora) = 0.38 za stator, 0.6 za rotor (vzeto iz `izgube_fem.txt` in primer1; za rotor predpostavljeno).
- A4: Sinusna zaokrožitev rotorja se modelira kot zr(α) = zr_min / cos(κ·α), κ ≈ 80/90 (kot v `primer1/python/KES_geometrija.py`, vrstica 178).
- A5: Število ovojev na pol rotorja = 20 (kot v primer1, vrstica 277); to bo verjetno optimizacijska spremenljivka.
- A6: Začetni interval za D_r = [100, 300] mm, L = [50, 250] mm. Omejitev l/D_r ≥ 0.4 zaradi tehnoloških razlogov.
- A7: Optimizacijske spremenljivke (gen): {D_r, l/D_r, B_δ, B_ds, B_sy, J_Cu_s, J_Cu_r, N_s_rotor, q}. (Skupaj 9 spremenljivk; q je diskretna.)
- A8: Mreženje FEMM na privzeti `smartmesh(0)`.
- A9: Frekvenca v vogalni točki: f_c = n_c/60 · p = 7000/60 · 3 = **350 Hz**.

### 5.3 Zunanje odvisnosti

- **FEMM 4.2** namestitev (Windows: privzeto `C:\femm42\`).
- **Python ≥ 3.10** z naslednjimi paketi (zaklenjenimi v `requirements.txt`):
  - `pyfemm`, `numpy`, `scipy`, `matplotlib`, `pymoo`, `pandas`, `pyyaml`.
- **Priloženi viri**:
  - `BH.txt` (B-H za M330-35A — iz `primer1/python/`).
  - `izgube_fem.txt` (tabela izgub M330-35A — koren projekta).
  - `SURA M330-35A.pdf` (referenca).
  - `Design of rotating electrical machines SECOND EDITION.pdf` (učbenik).
  - `manual.pdf`, `octavefemm.pdf` (FEMM dokumentacija).
  - `stroj.JPG`, `utor.JPG` (referenčni sliki za geometrijo).

### 5.4 Okoljske / regulativne omejitve
Ni relevantnih regulativnih zahtev (akademski projekt).

---

## 6. Specifikacije vmesnikov

### 6.1 Zunanji vmesniki

#### 6.1.1 FEMM
- **Tip**: krmiljenje preko `pyFEMM` Python klicev (notranja ovojnica za FEMM ALC ukaze).
- **Ključne funkcije**: `openfemm()`, `newdocument(0)`, `mi_probdef`, `smartmesh`, `mi_addnode`, `mi_addsegment`, `mi_drawarc`, `mi_createradius`, `mi_addmaterial`, `mi_addbhpoint`, `mi_addcircprop`, `mi_setblockprop`, `mi_addboundprop`, `mi_setarcsegmentprop`, `mi_moverotate`, `mi_analyze`, `mi_loadsolution`, `mo_getb`, `mo_getpointvalues`, `mo_blockintegral`, `mo_getcircuitproperties`, `mi_saveas`, `closefemm()`.
- **Datoteke**: `.fem` (vhod), `.ans` (rezultat). Imenovanje: `<design_id>_<scenario>.fem`, kjer scenario ∈ {`unbuilt`, `excited`, `noload`, `torque`, `Ld`, `Lq`}.

#### 6.1.2 CLI
- **Vstopna točka**: `python -m src.main [opcije]`.
- **Argumenti** *(assumed obseg)*:
  - `--config <pot>`: pot do YAML konfiguracije.
  - `--skip-ga`: preskoči GA, naloži zadnjo Pareto fronto.
  - `--design-id <ID>`: zagon FEMM-a le za eno rešitev.
  - `--seed <int>`: deterministična ponovljivost.

### 6.2 Notranji vmesniki

#### 6.2.1 Razred `MotorDesign`
```python
@dataclass
class MotorDesign:
    # vhod (gen)
    D_r: float          # premer rotorja [mm]
    l_to_D: float       # razmerje L_r / D_r [-]
    B_delta: float      # ciljna gostota mag. pretoka v zr. reži [T]
    B_ds: float         # B v statorskem zobu [T]
    B_sy: float         # B v statorskem jarmu [T]
    J_cu_s: float       # tokovna gostota stator [A/mm²]
    J_cu_r: float       # tokovna gostota rotor [A/mm²]
    N_r: int            # ovoji na rotorski pol
    q: float            # utori/(pol·fazo), ∈ {1, 1.5, 2, 2.5, 3}

    # izračunano
    L_r: float
    Q_s: int
    N_s: int            # ovoji statorja
    Z_q: int
    b_ds: float
    h_ds: float
    h_ys: float
    delta: float
    b_dr: float
    h_yr: float
    I_m: float          # vzbujalni tok [A]
    I_n: float          # nazivni statorski tok (fazni) [A]
    P_Fe: float
    P_Cu_s: float
    P_Cu_r: float
    eta: float
    V_active: float

    # validacijska zastavica
    feasible: bool
```

#### 6.2.2 GA <-> Analitični model
- GA podaja vektor `x` ∈ ℝ⁹ (zadnja komponenta razdeljena na 5 košev za diskretno q).
- `evaluate(x) -> (f1, f2, g)` kjer:
  - `f1 = 1/eta - 1` (minimizacija),
  - `f2 = V_active` [m³] (minimizacija),
  - `g` = vektor kršitev trdih omejitev (pozitivne vrednosti = kršitev).

### 6.3 Podatkovni modeli

#### 6.3.1 `results.csv`
| Stolpec | Tip | Enota | Opis |
|---------|-----|-------|------|
| `design_id` | str | — | Identifikator rešitve |
| `D_r` | float | mm | Premer rotorja |
| `L_r` | float | mm | Aktivna dolžina paketa |
| `q` | float | — | Utori/(pol·fazo) |
| `Q_s` | int | — | Število statorskih utorov |
| `N_s` | int | — | Ovoji statorja |
| `N_r` | int | — | Ovoji na rotorski pol |
| `I_m` | float | A | Vzbujalni tok |
| `I_n` | float | A | Nazivni statorski (fazni) tok |
| `P_Fe` | float | W | Skupne izgube v železu |
| `P_Cu` | float | W | Skupne izgube v bakru |
| `eta_analit` | float | — | Analitično izračunan η |
| `eta_femm` | float | — | η iz FEMM (P_Fe + P_Cu_femm) |
| `U_ind_analit` | float | V_RMS | Analitično U_ind pri 7000 rpm |
| `U_ind_femm` | float | V_RMS | FEMM U_ind |
| `M_avg` | float | Nm | Povprečni navor |
| `M_ripple` | float | % | Valovitost navora |
| `B_max_zob` | float | T | Maks. B v statorskem zobu |
| `B_max_jarem` | float | T | Maks. B v statorskem jarmu |

#### 6.3.2 `inputs.yaml` (privzeti vhodni podatki)
```yaml
m: 3
p: 3
U_f: 200          # V_RMS, fazna
P_c: 50000        # W
n_c: 7000         # rpm
n_max: 14000      # rpm
material: M330-35A
cos_phi: 0.95
K_Cu_s: 0.38
K_Cu_r: 0.6
T_del: 80         # °C
J_cu_s_max: 10    # A/mm²
J_cu_r_max: 5     # A/mm²
ga:
  pop_size: 100
  n_gen: 50
  seed: 42
femm:
  smartmesh: 0
  rotor_step_deg: 5
```

### 6.4 Ukazi / opkodi
Ni relevantno (ni vgrajeni sistem niti lasten protokol).

---

## 7. Operativni postopki

### 7.1 Namestitev
1. Namestiti **FEMM 4.2** (Windows): https://www.femm.info/wiki/Download (privzeta pot `C:\femm42\`).
2. Namestiti **Python ≥ 3.10** (https://python.org).
3. Iz korena projekta:
   ```powershell
   python -m venv .venv
   .venv\Scripts\Activate.ps1
   pip install -r requirements.txt
   ```
4. Preveriti, da `python -c "import femm; femm.openfemm()"` odpre FEMM okno.

### 7.2 Konfiguracija
- Privzeti vhodni podatki so v `src/inputs.py` ali `inputs.yaml`.
- Za spreminjanje GA parametrov (populacija, generacije, seed): urediti `inputs.yaml`.
- Za spremembo materiala: priložiti novo `BH.txt` in tabelo izgub; popraviti `material:` polje.

### 7.3 Normalno delovanje
1. Zagon polnega cevovoda:
   ```powershell
   python -m src.main --config inputs.yaml
   ```
2. Program v vrstnem redu:
   - Naloži vhodne podatke.
   - Zaganja NSGA-II GA (logiranje generacij v `outputs/ga_log.csv`).
   - Shrani Pareto fronto v `outputs/pareto.npz` in graf v `outputs/pareto.png`.
   - Izbere 5 rešitev → `outputs/selected5.json`.
   - Za vsako od 5 rešitev:
     - Odpre FEMM, izriše geometrijo, shrani `<design_id>_excited.fem`.
     - Izvede prosti tek, shrani `<design_id>_noload.fem` + slike.
     - Izvede simulacijo navora, shrani `<design_id>_torque.fem` + slike.
   - Zapolni `results.csv` in zapre FEMM.
3. Pregled rezultatov v `outputs/`.

### 7.4 Delni zagoni
- Samo GA (brez FEMM): `python -m src.main --skip-femm`.
- Samo FEMM za eno rešitev: `python -m src.main --skip-ga --design-id D03`.
- Samo izris geometrije (vizualna verifikacija): `python -m src.main --skip-ga --design-id D03 --draw-only`.

### 7.5 Vzdrževanje
- Izhodne datoteke (`outputs/`) **niso** verzionirane (vključno v `.gitignore`).
- Za reset stanja: `Remove-Item outputs/ -Recurse -Force`.

### 7.6 Postopki obnove
- **FEMM se obesi**: ubij `femm.exe` v Task Managerju, zaženi ponovno z `--skip-ga --design-id <ID>` od zadnje uspešne rešitve.
- **GA pade na pol poti**: `outputs/ga_log.csv` vsebuje zadnjo populacijo; pomožno skripto za nadaljevanje *(out of scope za to fazo)*.
- **Napaka pri risanju**: dodaj `--debug` zastavico za interaktivno odpiranje FEMM (brez `closefemm()`), preveri vizualno.

---

## 8. Verifikacija in validacija

### 8.1 Verifikacija faze 1 (analitični izračun)

| Test ID | Funkcionalnost | Postopek | Pričakovani rezultat |
|---------|----------------|----------|----------------------|
| TC-A.1 | Geometrija | Zaženi `analytical.compute_geometry()` za vhod 200 V / 50 kW / 7000 rpm / q = 2 | D_r, l, Q_s, N_s, b_ds, h_ds, h_ys, delta vrnjeni v razumnih mejah (npr. D_r ∈ [100, 300] mm); brez negativnih vrednosti. |
| TC-A.2 | EMETOR k_w1 | Klic `winding_factor(q)` za q ∈ {1, 1.5, 2, 2.5, 3} | Vrne k_w1 ∈ [0.86, 1.0]; za q = 2 vrne 0.966. |
| TC-A.3 | Gostota toka | Preveri J_Cu_s in J_Cu_r v izračunani geometriji | J_Cu_s ≤ 10, J_Cu_r ≤ 5. |
| TC-B.1 | B_δ(I_m) | Plot B_δ vs. I_m za fiksen θ_m | Linearno do nasičenja, nato saturacijska krivulja. |
| TC-C.1 | Polinom izgub | Preveri f(B = 1, f = 50) vs. tabela `izgube_fem.txt` | Relativno odstopanje < 5 %. |
| TC-D.1 | P_Fe | P_Fe statorja pri B_max = 1.5 T, f = 350 Hz | V razumnem območju (200–1500 W za 50 kW stroj); ločeno za zob in jarem. |
| TC-E.1 | P_Cu | P_Cu pri T_del = 80 °C | Skladno z R(80°C) = R(20°C) · (1 + α·60). |
| TC-F.1 | η | η = P_c / (P_c + P_Fe + P_Cu) | η > 0.9 za ciljno točko. |

### 8.2 Verifikacija faze 2 (optimizacija)

| Test ID | Funkcionalnost | Postopek | Pričakovani rezultat |
|---------|----------------|----------|----------------------|
| TC-G.1 | GA konvergenca | Zaženi GA s seed = 42, pop = 100, gen = 50 | Hipervolumen Pareto fronte se v zadnjih 10 generacijah ne spreminja > 1 %. |
| TC-G.2 | Determinizem | Dvakratni zagon z istim seedom | Identična Pareto fronta. |
| TC-H.1 | Pareto fronta | Pregled `pareto.png` | Vidna negativno-konveksna krivulja v ravnini (V, 1/η). |
| TC-H.2 | Trde omejitve | Preveri vsakega osebka na fronti | Vsi imajo `feasible = True`. |
| TC-K.1 | Izbor 5 rešitev | `pareto.select_five()` | Vrne 5 različnih `MotorDesign`; 1× max η, 1× min V, 3× vmesni. |

### 8.3 Verifikacija faze 3 (FEMM)

| Test ID | Funkcionalnost | Postopek | Pričakovani rezultat |
|---------|----------------|----------|----------------------|
| TC-I.1 | Izris geometrije | Zaženi `femm_model.build(design)` za D03 | `.fem` datoteka, vizualno enaka `stroj.JPG` (rotor s sinusno zaokroženimi zobmi, 6 polov, vzbuj. navitje na vsakem polu); FEMM brez `node overlap` napak. |
| TC-J.1 | Materiali in mreženje | Po izrisu zaženi `mi_analyze()` | Konča v < 30 s, brez napak; rešitev shranjena v `.ans`. |
| TC-L.1 | Prosti tek — U_ind | Zaženi `femm_sim.no_load(design)` | Ψ_A(θ) graf ima en period čez 360° / p = 120°; U_ind_amp ≤ 0.95 · √2 · 200 ≈ 269 V (amp.). |
| TC-L.2 | U_ind vs. analitično | Primerjaj FEMM U_ind z analitičnim | |Δ| < 10 % (NFR-3.1). |
| TC-M.1 | Navor | Zaženi `femm_sim.torque(design)` | navor(θ) graf z osnovno harmonsko; povprečni navor blizu nazivnemu M = P_c/ω_meh = 50000/(2π·7000/60) ≈ 68.2 Nm. |
| TC-M.2 | FFT navora | Spekter navora | 1. harmonska dominantna; 2. harmonska < 30 % osnovne *(assumed)*. |

### 8.4 Sledljivostna matrika

| Zahteva | Prioriteta | Test(i) | Status |
|---------|-----------|---------|--------|
| FR-1.1 | Must | TC-A.1 | Covered |
| FR-1.2 | Must | TC-A.2 | Covered |
| FR-1.3 | Must | TC-A.3 | Covered |
| FR-1.4 | Must | TC-B.1 | Covered |
| FR-1.5 | Must | TC-C.1 | Covered |
| FR-1.6 | Must | TC-D.1 | Covered |
| FR-1.7 | Must | TC-E.1 | Covered |
| FR-1.8 | Must | TC-F.1 | Covered |
| FR-2.1 | Must | TC-G.1 | Covered |
| FR-2.2 | Must | TC-H.1 | Covered |
| FR-2.3 | Must | TC-H.2 | Covered |
| FR-2.4 | Must | TC-K.1 | Covered |
| FR-2.5 | Should | TC-H.1 | Covered |
| FR-3.1 | Must | TC-I.1 | Covered |
| FR-3.2 | Must | TC-I.1 | Covered |
| FR-3.3 | Must | TC-I.1 | Covered |
| FR-3.4 | Must | TC-I.1 | Covered |
| FR-3.5 | Must | TC-I.1 | Covered |
| FR-3.6 | Must | TC-J.1 | Covered |
| FR-3.7 | Must | TC-L.1 | Covered |
| FR-3.8 | Must | TC-L.1, TC-L.2 | Covered |
| FR-3.9 | Must | TC-M.1 | Covered |
| FR-3.10 | Must | TC-M.2 | Covered |
| FR-4.1 | Must | (artifakti generirani v TC-L.1, TC-M.1) | Covered |
| FR-4.2 | Must | (`results.csv` po polnem zagonu) | Covered |
| FR-4.3 | Should | (ročna verifikacija) | Covered |
| FR-5.1 | Must | (code review) | Covered |
| FR-5.2 | Should | (code review) | Covered |
| NFR-1.1 | Must | benchmark analytical.evaluate | Covered |
| NFR-1.2 | Should | čas GA zagona | Covered |
| NFR-1.3 | Should | čas FEMM simulacije | Covered |
| NFR-2.1 | Must | TC-G.2 | Covered |
| NFR-2.2 | Must | preveri prisotnost `pareto.npz` | Covered |
| NFR-3.1 | Must | TC-L.2 | Covered |
| NFR-3.2 | Must | (TC-D.1 + TC-E.1 vs. FEMM) | Covered |
| NFR-4.1 | Should | code review | Covered |
| NFR-4.2 | Should | (ročno na macOS preko Wine) | Covered |
| NFR-5.1 | Must | code review + poročilo | Covered |

---

## 9. Vodnik za odpravljanje težav

| Simptom | Verjetni vzrok | Diagnostični koraki | Korektivni ukrep |
|---------|---------------|---------------------|------------------|
| FEMM odpre okno, a `mi_analyze()` vrne napako "node overlap" | Iterativni izris utora se ni ustavil v pravem trenutku → točke se prekrivajo | Izriši `<design_id>_excited.fem` ročno v FEMM GUI; preveri vse `mi_addnode` klice | Zmanjšaj inkrement v `while spodnja_bds > bds: utor[0] += 0.001` (npr. 0.0005); dodaj toleranco. |
| FEMM mreženje traja > 5 min | `smartmesh(0)` + drobna geometrija (npr. zelo tanek zračna reža) | Zaženi FEMM ročno, preveri velikost mreže | Povečaj `smartmesh` parameter (npr. `smartmesh(1)`); dodaj eksplicitne `mi_setblockprop` `meshsize` parametre. |
| `pyfemm` se obesi po `femm.openfemm()` | FEMM dialog (npr. licenčni warning) blokira | Glej procesno listo, ali se `femm.exe` izvaja | Zapri vse FEMM instance v Task Managerju, zaženi FEMM ročno enkrat, sprejmi dialog. |
| GA vrne prazno Pareto fronto | Trde omejitve preveč stroge | Zaženi GA brez U_ind omejitve, izpis kršitev po populaciji | Mehkejše začetne meje; uvedi penalty namesto hard cutoff. |
| U_ind v FEMM je 2× večji od analitičnega | Robni pogoj `A=0` ni nastavljen ali pa krožni obseg ne meji na statorski jarem | Preveri `mi_setarcsegmentprop(1, "A=0", 0, 0)` klic; v FEMM GUI preveri "boundary" oznake | Eksplicitno selektuj zunanji statorski rob in nastavi A=0. |
| Navor ima velike višje harmonske, 2. harmonska > 50 % | Statorski zobje nimajo razširjene glave (cevlja) | Vizualno preveri `<design_id>_torque.fem` | Aktiviraj `mi_createradius` korak v `femm_model.py`. |
| Polinom izgub vrne negativne vrednosti pri majhnih B | `poly23` ne dobro modelira majhnih vrednosti | Plot polinoma vs. tabele | Spremeni v `polyfit` z višjim redom ali uporabi `RBFInterpolator`. |
| Program konča, a `results.csv` manjka stolpec | En od FEMM klicev je padel "tiho" | Preveri `outputs/run.log` | Popravi try/except, da ne požira napak. |

---

## 10. Priloga

### 10.1 Konstante in privzete vrednosti

| Konstanta | Vrednost | Enota | Vir |
|-----------|----------|-------|-----|
| Specifična prevodnost bakra σ_Cu | 34·10⁶ | S/m | `izgube_fem.txt` |
| Temp. koeficient bakra α_Cu | 0.00381 | K⁻¹ | `izgube_fem.txt` |
| Delovna temperatura T_del | 80 | °C | `izgube_fem.txt` |
| Faktor zapolnitve utora K_Cu (stator) | 0.38 | — | `izgube_fem.txt` |
| Faktor zapolnitve K_Cu (rotor) | 0.6 *(assumed)* | — | predpostavka |
| Gostota železa ρ_Fe | 7650 | kg/m³ | `izgube_fem.txt` |
| k_Fe_zob (povečanje izgub v zobu) | 2 | — | `izgube_fem.txt` |
| k_Fe_jarem | 1.6 | — | `izgube_fem.txt` |
| k_fe (polnilni f. lameliranja) | 0.97 | — | primer1 |
| Carterjev koef. K_c (začetna ocena) | 1.17 | — | primer1 |
| Max J_Cu stator | 10 | A/mm² | naloga, dodatno navodilo |
| Max J_Cu rotor | 5 | A/mm² | naloga, dodatno navodilo |

### 10.2 Tabela simbolov (izvleček — polna verzija v poročilu)

| Simbol | Polno ime | Enota |
|--------|-----------|-------|
| D_r | Premer rotorja (zunanji) | mm |
| L_r | Aktivna dolžina paketa | mm |
| δ | Zračna reža (minimalna, na sredini pola) | mm |
| Q_s | Število statorskih utorov | — |
| q | Število utorov na pol in fazo | — |
| k_w1 | Faktor navitja osnovne harmonske | — |
| N_s | Število ovojev statorja na fazo | — |
| Z_q | Število ovojev v enem utoru | — |
| b_ds | Širina statorskega zoba | mm |
| h_ds | Višina statorskega zoba | mm |
| h_ys | Višina statorskega jarma | mm |
| b_dr | Širina rotorskega zoba (osnova) | mm |
| h_yr | Višina rotorskega jarma | mm |
| B_δ | Amplituda B v zračni reži | T |
| B_ds | B v statorskem zobu | T |
| B_sy | B v statorskem jarmu | T |
| I_m | Vzbujalni tok (rotor) | A |
| I_n | Nazivni statorski tok (faz., RMS) | A |
| J_Cu_s | Tokovna gostota stator | A/mm² |
| J_Cu_r | Tokovna gostota rotor | A/mm² |
| P_Fe | Izgube v železu | W |
| P_Cu | Izgube v bakru | W |
| η | Izkoristek | — |
| V_active | Volumen aktivnega dela | m³ |
| f_c | Frekvenca v vogalni točki | Hz |
| ω_meh | Mehanska kotna hitrost | rad/s |

### 10.3 Ključne enačbe (kratek izvleček)

```
ω_meh   = 2π·n_c / 60
ω_el    = ω_meh · p
f_c     = n_c · p / 60                                = 350 Hz
M       = P_c / ω_meh                                 ≈ 68.2 Nm
τ_p     = π · D_r / (2·p)                              ← polov korak
δ_min   ≈ 1/K_c · 4·1e-7 · τ_p · A_s / B_δ            ← zr. reža (Pyrhönen §3)
N_s     = √2 · U_f / (ω_el · α · B_δ · τ_p · L_r · k_w1)
A_s     = (M / 2 / V_r) · √2 / (B_δ · cosφ)            ← tokovna obloga
V_r     = π · (D_r/2)² · L_r
R(T)    = R₀ · (1 + α_Cu · (T - 20))
P_Cu    = 3 · I_n² · R_s(T)  +  I_m² · R_r(T)
P_Fe,zob   = k_Fe_zob   · m_Fe,zob   · f(B_zob, f_c)
P_Fe,jarem = k_Fe_jarem · m_Fe,jarem · f(B_jarem, f_c)
η       = P_c / (P_c + P_Fe + P_Cu)
U_ind   = max|dΨ_A/dθ| · ω_meh
```

### 10.4 Imenovanje datotek

- `outputs/fem/D01_excited.fem` … `D05_excited.fem` — model s preračunanim vzbujanjem.
- `outputs/fem/D01_noload.fem` — prosti tek (po vrtenju).
- `outputs/fem/D01_torque.fem` — simulacija navora.
- `outputs/figures/D01/uind.png`, `navor.png`, `fft.png`.
- `outputs/pareto.png`, `outputs/pareto.npz`, `outputs/selected5.json`, `outputs/results.csv`.

### 10.5 Reference

- **Pyrhönen, J., Jokinen, T., Hrabovcová, V.**: *Design of Rotating Electrical Machines*, 2nd ed., Wiley. *(priloženo: `Design of rotating electrical machines SECOND EDITION.pdf`)*
- **Meeker, D.**: *FEMM 4.2 User's Manual*. *(priloženo: `manual.pdf`)*
- **Meeker, D.**: *Octave-FEMM Interface Documentation*. *(priloženo: `octavefemm.pdf`)*
- **EMETOR**: korelacije q ↔ Q_s ↔ k_w1, https://www.emetor.com/ (online tabela; uporabljeno za določitev q in k_w1).
- **Sura/Cogent**: *Material data sheet M330-35A*. *(priloženo: `SURA M330-35A.pdf`)*
- **Deb, K., Pratap, A., Agarwal, S., Meyarivan, T.**: *A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II*, IEEE Trans. Evol. Comp., 2002. *(referenca za GA)*
- **primer1** (interna referenca): `primer1/python/KES.py`, `primer1/python/KES_geometrija.py`, `primer1/porocilo.pdf` — predloga za FEMM izris in arhitekturo Python kode.
- **primer2** (interna referenca): `primer2/koda/*.m` — MATLAB ekvivalent (za navzkrižno preverjanje analitičnih enačb in izračuna izgub).
- **primer3** (interna referenca): `primer3/Seminar KES - Vukovic.pdf` — letošnja predloga poročila (8-polni stroj).

---

## 11. Povezani dokumenti

- `[[primer1-porocilo]]` — referenčno poročilo za 6-polni stroj iz prejšnjega leta.
- `[[primer3-porocilo]]` — referenčno poročilo za 8-polni stroj iz tega leta (struktura poročila).
- `[[Sinhronc_opis]]` — izvorni opis projekta in seznam nalog A–N.
