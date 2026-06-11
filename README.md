# 6-polni sinhronski motor z vzbujalnim navitjem

Orodje za **dimenzioniranje**, **optimizacijo** in **FEMM verifikacijo**
6-polnega sinhronskega stroja (50 kW @ 7000 rpm). Seminarska naloga pri
predmetu *Konstruiranje električnih strojev*, FE UL.

Potek dela:

1. **Analitika** — preračun geometrije, magnetnih količin in izgub (po Pyrhönen et al.).
2. **Optimizacija** — NSGA-II (`pymoo`): maks. izkoristek vs. min. volumen; iz Pareto fronte izbere 5 reprezentativnih strojev.
3. **FEMM verifikacija** — za vsak stroj prosti tek ($U_\text{ind}$, $B_\max$) in navorna analiza (FFT, THD).
4. **Grafi** — Pareto fronta, primerjava analitika ↔ FEMM, $U_\text{ind}(\theta)$, $M(\theta)$.

## Namestitev

- **Python 3.11+**
- **FEMM 4.2** ([femm.info](https://www.femm.info/)) — potreben le za FEMM verifikacijo.

```powershell
pip install -r requirements.txt
```

## Zagon

Najlažje — celoten potek (vse tri korake) z enim klikom:

```powershell
.\run_all.bat
```

Ali ročno, korak za korakom:

```powershell
python -m src.main --config inputs.yaml      # analitika + GA + izbor 5 rešitev
python -m src.run_femm_pareto --step 5       # FEMM simulacije (~18 min)
python -m src.plot_pareto_femm               # zbirni in primerjalni grafi
```

Vsi rezultati (`.csv`, `.json`, `.npz`, grafi, FEMM modeli) se zapišejo v `outputs/`.

Koristni argumenti za `run_femm_pareto`:
`--step` (kotni korak vrtenja), `--only <oznaka>` (samo en stroj),
`--skip-noload` / `--skip-torque` (preskoči simulacijo).

## Struktura

```
src/
├── inputs.py / losses.py        # podatki, polinom izgub p_Fe(B,f)
├── analytical.py                # analitični izračun
├── optimization.py / pareto.py  # NSGA-II + izbor 5 rešitev
├── femm_model.py / femm_sim.py  # FEMM geometrija + simulacije
├── post.py / plot.py            # post-procesiranje + grafi
├── main.py                      # koraka 1+2
├── run_femm_pareto.py           # koraka 3+4 (z auto-fallback ob napaki FEMM)
└── plot_pareto_femm.py          # zbirni grafi
data/   inputs.yaml   tests/   requirements.txt
```

## Testi

```powershell
python -m pytest tests/ --ignore=tests/smoke_sim.py --ignore=tests/smoke_femm.py
```

Smoke testi zahtevajo nameščeno FEMM 4.2 in se zato privzeto preskočijo.

## Reference

- Pyrhönen, Jokinen, Hrabovcová: *Design of Rotating Electrical Machines*, 2. izd., Wiley, 2014.
- Meeker, D.: *FEMM — User's Manual*, v. 4.2, [femm.info](https://www.femm.info/).
- Deb et al.: "NSGA-II", *IEEE TEC*, 2002.
- EMETOR ([emetor.com](https://www.emetor.com/)); Sura/Cogent *M330-35A* data sheet.
