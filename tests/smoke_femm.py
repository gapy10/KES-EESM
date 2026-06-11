"""Smoke FEMM test - zaženi v foreground, izpis na konec."""
import sys, math, json, time, femm
from src.femm_model import kill_stale_femm

I_m_default = 21.0
femm_path = "outputs/fem/smoke_test.fem"

try:
    data = json.loads(open('outputs/selected5.json', encoding='utf-8').read())
    mid = next(d for d in data if d['label'] == 'mid')
    I_m = mid['electrical']['I_m_A']
except Exception:
    I_m = I_m_default
    mid = None

print(f'I_m = {I_m:.2f} A', flush=True)
kill_stale_femm()
femm.openfemm()
femm.opendocument(femm_path)
femm.mi_setcurrent('DC', I_m)
t0 = time.time()
femm.mi_analyze()
femm.mi_loadsolution()
print(f'Analyzed in {time.time()-t0:.1f}s', flush=True)

if mid:
    delta = mid['geometry_mm']['delta']
    D_si = mid['geometry_mm']['D_si']
    Rr = D_si / 2 - delta
    R_test = Rr + delta / 2
else:
    R_test = 85.0

bvals = []
for k in range(60):
    ang = k * 2 * math.pi / 60
    bx, by = femm.mo_getb(R_test * math.cos(ang), R_test * math.sin(ang))
    bvals.append((bx ** 2 + by ** 2) ** 0.5)

print(f'B v zr. rezi (R={R_test:.2f} mm):', flush=True)
print(f'  min = {min(bvals):.3f} T', flush=True)
print(f'  max = {max(bvals):.3f} T', flush=True)
print(f'  avg = {sum(bvals)/len(bvals):.3f} T', flush=True)
print(f'  variansa = {sum((b-sum(bvals)/len(bvals))**2 for b in bvals)/len(bvals):.4f}', flush=True)
femm.closefemm()
print('done', flush=True)
