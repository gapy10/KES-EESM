# KES v7 — 6-polni sinhronski motor z vzbujalnim navitjem

Programsko orodje za **analitično dimenzioniranje**, **večkriterijsko
optimizacijo** (NSGA-II: maks. izkoristek, min. volumen) in **FEMM
verifikacijo** 6-polnega sinhronskega stroja. Seminarska naloga pri
predmetu *Konstruiranje električnih strojev*, FE UL.

## Glavne značilnosti

- Analitičen preračun geometrije, magnetnih količin in izgub po Pyrhönen et al.
- Polinomska aproksimacija specifičnih izgub za jeklo **M330-35A**
  ($p_{Fe}(B, f)$, 12 koeficientov, RMSE ≈ 0,76 W/kg).
- **NSGA-II GA** (preko `pymoo`): 9 spremenljivk, 2 cilja
  ($1/\eta - 1$, $V_{\text{active}}$), 3 omejitve.
- Izbor **5 reprezentativnih rešitev** s Pareto fronte z enakomernim
  vzorčenjem po loku.
- **Parametričen FEMM model** (rotor s sinusno zaokrožitvijo, 36 utorov,
  tri faze, vzbujalno navitje, B-H krivulja).
- **Avtomatska fallback rutina**: če FEMM gradnja kakega stroja spodleti,
  `run_femm_pareto.py` samodejno preide na najbližjega soseda na Pareto
  fronti — brez ročnega posega.
- Generira Pareto grafe (η in P_loss), primerjalne grafe (η analitika ↔
  FEMM, V po designih), grafa $U_{\text{ind}}(\theta)$ in $M(\theta)$,
  FFT spektre.

## Zahteve

- **Python 3.11+**
- **FEMM 4.2** ([femm.info](https://www.femm.info/)) nameščen lokalno —
  potreben za korake 3 in 4 (faze FEMM verifikacije).
- Knjižnice iz `requirements.txt`:

```powershell
pip install -r requirements.txt
```

Glavne odvisnosti: `numpy`, `scipy`, `pandas`, `matplotlib`, `pyyaml`,
`pymoo`, `pyfemm`.

## Zagon

### 1. Faze 1+2 — analitika + GA + Pareto + izbor 5 rešitev

```powershell
python -m src.main --config inputs.yaml
```

Privzeti parametri NSGA-II: `--pop 100 --gen 50 --seed 42`. Brez argumenta
`--config` se uporabijo privzete vrednosti iz `src/inputs.py`.

Izhodi v `outputs/`:
- `pareto.npz` — celotna Pareto fronta (geni in cilji).
- `pareto.png` — graf fronte v ravnini (V, η).
- `selected5.json` — 5 izbranih rešitev.
- `results.csv` — vsi Pareto kandidati.

### 2. Faze 3+4 — FEMM simulacije in post-procesiranje

```powershell
python -m src.run_femm_pareto --step 5
```

Za vsako od 5 izbranih rešitev izvede **prosti tek** (izračun $U_{\text{ind}}$,
B_max v zobu in jarmu) in **navorno analizo** (M_1.harm, FFT, THD). Tipičen
celokupen zagon: ~18 min.

Argumenti:
- `--step 5` — kotni korak vrtenja v stopinjah (manjše = natančneje, počasneje).
- `--only mid` — samo ena oznaka (`max_eta`, `high_eta_mid`, `mid`,
  `low_eta_mid`, `min_V`); rezultati se zlijejo z obstoječim povzetkom.
- `--max-retries 5` — največje število alternativ pri napaki FEMM gradnje.
- `--skip-noload`, `--skip-torque` — preskoči ustrezno simulacijo.

Izhodi v `outputs/`:
- `femm_summary.csv` — povzetek vseh 5 rešitev (s stolpcema `pareto_idx`,
  `replaced`).
- `femm_pareto_results.json` — polni FEMM rezultati.
- `fem/*.fem` — FEMM modeli.
- `figures/D01_max_eta/`, …, `figures/D05_min_V/` — grafe v 300 dpi.

### 3. Pareto grafe z 5 FEMM-validiranimi rešitvami

```powershell
python -m src.plot_pareto_femm
```

Generira:
- `outputs/pareto_femm.png` — Pareto fronta (η na y-osi) z 5 označenimi
  kandidati.
- `outputs/pareto_femm_loss.png` — minimizacijski prikaz ($P_{\text{loss}}$
  na y-osi).
- `outputs/eta_comparison.png` — stolpčni graf η analitično vs FEMM.
- `outputs/V_comparison.png` — stolpčni graf $V_{\text{active}}$.

## Struktura projekta

```
KES v7/
├── data/
│   └── BH_M330-35A.txt              # B-H krivulja za FEMM
├── src/
│   ├── inputs.py                     # privzeti podatki, EMETOR tabela
│   ├── losses.py                     # polinom $p_{Fe}(B, f)$
│   ├── analytical.py                 # `analyze()` — analitični izračun
│   ├── optimization.py               # NSGA-II ovojnica preko pymoo
│   ├── pareto.py                     # izbor 5 reprezentativnih rešitev
│   ├── femm_model.py                 # parametrični izris geometrije
│   ├── femm_sim.py                   # prosti tek, navor
│   ├── post.py                       # P_Fe iz FEMM, primerjava, THD
│   ├── plot.py                       # vsi grafi @ 300 dpi
│   ├── main.py                       # faze 1+2 CLI
│   ├── run_femm_pareto.py            # faze 3+4 CLI (z retry-with-neighbor)
│   └── plot_pareto_femm.py           # zbirni Pareto + primerjalni grafe
├── tests/                            # pytest enote in smoke testi
├── inputs.yaml                       # privzeta konfiguracija
└── requirements.txt
```

## Testi

```powershell
python -m pytest tests/ --ignore=tests/smoke_sim.py --ignore=tests/smoke_femm.py
```

Smoke testi (`smoke_sim.py`, `smoke_femm.py`) zahtevajo nameščeno FEMM 4.2
in se zato preskakujejo v privzetem zagonu pytest.

## Reference

- **Pyrhönen, J., Jokinen, T., Hrabovcová, V.** *Design of Rotating Electrical
  Machines*, 2nd ed., Wiley, 2014.
- **Meeker, D.** *Finite Element Method Magnetics — User's Manual*, v. 4.2,
  [femm.info](https://www.femm.info/).
- **Deb et al.** "A Fast and Elitist Multiobjective Genetic Algorithm:
  NSGA-II", *IEEE TEC*, 2002.
- **EMETOR** — [emetor.com](https://www.emetor.com/) (faktorji navitja).
- **Sura/Cogent** — *M330-35A* data sheet.
