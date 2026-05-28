import femm
import math
import KES_geometrija as KES
import plot as plt
import gui


def postavi_model():
    femm.openfemm()
    femm.newdocument(0)
    femm.mi_probdef(0,'millimeters','planar',1E-8, 110, 30, 0)
    femm.smartmesh(0)
    KES.narisi_geometrijo()
    KES.koncaj_geometrijo()
    femm.mi_zoomnatural()
    femm.mi_saveas("KES_nevzbujan.fem")
    #femm.closefemm()

# Funkcije za analizo

# Doloci vzbujalni tok - le za dolocitev
def vzbujalni_tok(Im, Bamp):
    femm.mi_addcircprop("DC", Im, 1)
    I_m = Im
    tocka = [8, 95.2]
    femm.mi_analyze()
    femm.mi_loadsolution()
    bx, by = femm.mo_getb(tocka[0], tocka[1])
    b_amp = (bx**2 + by**2)**0.5

    while math.fabs(b_amp - Bamp) > 0.001:
        if (b_amp - Bamp) < 0:
            I_m = I_m + 0.1
        else:
            I_m = I_m - 0.1
        femm.mi_setcurrent("DC", I_m)
        femm.mi_analyze()
        femm.mi_loadsolution()
        bx, by = femm.mo_getb(tocka[0], tocka[1])
        b_amp = (bx**2 + by**2)**0.5 
        print(I_m, b_amp)
    return I_m

def nastavi_vzbujanje(I_m):
    #femm.openfemm()
    #femm.opendocument("KES_JS.fem")
    femm.mi_saveas("KES_vzbujanje.fem")
    femm.mi_addcircprop("DC", I_m, 1)
    femm.mi_saveas("KES_vzbujanje.fem")
    femm.mi_analyze()
    femm.mi_loadsolution()
    return I_m
    #femm.closefemm()

def inducirana_napetost(id):
    #femm.openfemm()
    #femm.opendocument("KES_vzbujanje.fem")
    if id != 1:
        return
    print("Inducirana napetost")
    psi_t = []
    uind_t = []
    u_ind = []
    t = []
    inc = 5
    for i in range(0, 121, inc):
        femm.mi_selectgroup(1)
        femm.mi_moverotate(0,0,inc)
        femm.mi_analyze()
        femm.mi_loadsolution()
        tok, nap, psi = femm.mo_getcircuitproperties("A")
        psi_t.append(psi)
        t.append(i / 120 / 400)
    for i in range(len(psi_t) - 1):
        uind = (psi_t[i+1] - psi_t[i]) / (t[i+1] - t[i])
        uind_t.append([uind, t[i+1]/2 + t[i]])
        u_ind.append(uind)
    ret = "Model shranjen kot KES_vzbujanje.fem\n"
    global ret_str
    ret_str = ret
    plt.plt_uind(u_ind)

def prosti_tek(id):
    if id != 2:
        return
    print("Prosti tek")
    i_m = []
    u_ind = []
    inc = 5
    femm.mi_saveas("KES_prosti_tek.fem")
    femm.mi_selectgroup(1)
    femm.mi_moverotate(0,0, 25)
    femm.mi_clearselected()
    femm.mi_saveas("KES_prosti_tek.fem")
    for i in range(0, 241, inc):
        femm.mi_setcurrent("DC", i)
        femm.mi_analyze()
        femm.mi_loadsolution()
        tok, nap, psi1 = femm.mo_getcircuitproperties("A")
        femm.mi_selectgroup(1)
        femm.mi_moverotate(0,0, 5)
        femm.mi_analyze()
        femm.mi_loadsolution()
        tok, nap, psi2 = femm.mo_getcircuitproperties("A")
        uind = (psi2 - psi1) / (5/120/400)
        print(uind, i)
        u_ind.append(uind/math.sqrt(2))
        i_m.append(i)
        femm.mi_selectgroup(1)
        femm.mi_moverotate(0,0, -5)
    ret = "Model shranjen kot KES_prosti_tek.fem\n"
    global ret_str
    ret_str = ret
    plt.plt_kpt(u_ind)


def navor(id):
    if id != 4:
        return
    print("Navor")

    femm.mi_saveas("KES_navor.fem")
    Ia = 142.85 * math.sqrt(2)
    femm.mi_setcurrent("A", Ia)
    femm.mi_setcurrent("B", -Ia/2)
    femm.mi_setcurrent("C", -Ia/2)
    navor = []
    inc=5
    for i in range(0, 121, inc):
        if i:
            femm.mi_selectgroup(1)  
            femm.mi_moverotate(0,0,inc)
            femm.mi_clearselected()
        femm.mi_analyze()
        femm.mi_loadsolution()
        femm.mo_groupselectblock(1)
        torque = femm.mo_blockintegral(22)
        femm.mo_clearblock()
        navor.append(torque)
    femm.mi_saveas("KES_navor.fem")
    femm.mi_analyze()
    femm.mi_loadsolution()
    ret = "Model shranjen kot KES_navor.fem\n"
    global ret_str
    ret_str = ret
    plt.plt_M(navor)

def stator(id):
    if id != 3:
        return
    print("Statorsko napajanje")
    femm.mi_saveas(f"KES_statorsko_napajanje_0.fem")
    Ia = 142.85 * math.sqrt(2)
    femm.mi_setcurrent("A", Ia)
    femm.mi_setcurrent("B", -Ia/2)
    femm.mi_setcurrent("C", -Ia/2)

    femm.mi_saveas(f"KES_statorsko_napajanje_0.fem")
    femm.mi_analyze()
    femm.mi_loadsolution()

    femm.mi_saveas(f"KES_statorsko_napajanje_60.fem")
    femm.mi_selectgroup(1)
    femm.mi_moverotate(0,0, 60)
    femm.mi_clearselected()
    femm.mi_saveas(f"KES_statorsko_napajanje_60.fem")

    femm.mi_analyze()
    femm.mi_loadsolution()

    ret = "Modela shranjena kot KES_statorsko_napajanje_0.fem in KES_statorsko_napajanje_60.fem"
    global ret_str
    ret_str = ret

def induktivnost_d(id):
    if id != 5:
        return
    print("Magnetna induktivnost d")
    femm.mi_saveas("KES_Ld.fem")
    femm.mi_setcurrent("DC", 0)
    femm.mi_setprevious('KES_vzbujanje.ans', 2)

    # Lp
    Ia = 1
    femm.mi_setcurrent("A", Ia)
    femm.mi_setcurrent("B", -Ia/2)
    femm.mi_setcurrent("C", -Ia/2)

    femm.mi_saveas("KES_Ld.fem")

    femm.mi_analyze()
    femm.mi_loadsolution()

    tok, nap, psiA0 = femm.mo_getcircuitproperties("A")
    tok, nap, psiB0 = femm.mo_getcircuitproperties("B")
    tok, nap, psiC0 = femm.mo_getcircuitproperties("C")
    psi_d = 2/3*(psiA0 - psiB0/2 - psiC0/2)

    ret = "Ld = " + str(1000*psi_d) + " [mH]\n Model shranjen kot KES_Ld.fem"
    global ret_str
    ret_str = ret

def induktivnost_q(id):
    if id != 6:
        return
    print("Prečna induktivnost q")
    femm.mi_saveas("KES_Lq.fem")
    femm.mi_setcurrent("DC", 0)
    femm.mi_setprevious('KES_vzbujanje.ans', 2)
    Ib = math.sqrt(3)/2
    femm.mi_setcurrent("A", 0)
    femm.mi_setcurrent("B", Ib)
    femm.mi_setcurrent("C", -Ib)

    femm.mi_saveas("KES_Lq.fem")

    femm.mi_analyze()
    femm.mi_loadsolution()
    tok, nap, psiA = femm.mo_getcircuitproperties("A")
    tok, nap, psiB = femm.mo_getcircuitproperties("B")
    tok, nap, psiC = femm.mo_getcircuitproperties("C")
    psi_q = math.sqrt(3)/3*(psiB - psiC)

    ret = "Lq = " + str(1000*psi_q) + " [mH]\n Model shranjen kot KES_Lq.fem"
    global ret_str
    ret_str = ret




# macOS / Linux using Wine only - comment on Windows
winepath = '/opt/homebrew/bin/wine' 
femmpath = '/drive_c/femm42'

id, naloga = gui.run_gui()
#id = 1
#naloga = "Inducirana napetost"
ret_str = None

# Open FEMM
postavi_model()
Im = nastavi_vzbujanje(57)

inducirana_napetost(id)
prosti_tek(id)    
stator(id)
navor(id)
induktivnost_d(id)
induktivnost_q(id)

gui.show_value_window(naloga, ret_str)
femm.closefemm()

while True:
    pass

