import femm
import math

def zrcali_rotiraj(x, y, kot):
    kot = kot*math.pi/180
    x = -x
    x1 = x*math.cos(kot) - y*math.sin(kot)
    y1 = x*math.sin(kot) + y*math.cos(kot)
    return [x1, y1]

def narisi_geometrijo():
    # Nazivni podatki
    Pn = 120000
    pp = 3
    Udc = 800
    nn = 8000
    nm = 20000
    m = 3

    # Število utorov
    q = 2 # Qs = 18 utorov * q
    Qs = 2*pp*m*q
    Kws1 = 0.966 # emetor @ 36, 1@18

    # Dimenzije 
    Dr = 190
    l = 110
    ld = l/Dr
    Rr = Dr / 2
    Vr = math.pi*math.pow(Rr, 2)*l/1000000000
    omega_meh = (2*math.pi*nn/60)
    M = Pn / omega_meh
    tangencialna_obremenitev = M / 2 / Vr
    cosfi = 0.95
    Bzr = 1
    A_tokovna_obloga = tangencialna_obremenitev/Bzr*math.sqrt(2)/cosfi
    obseg_rotorja = math.pi*Rr*2
    tau_p = obseg_rotorja / 6
    Kc = 1.172
    zr_min = 1/Kc*4*(1e-7)*tau_p*A_tokovna_obloga/Bzr
    U1f = Udc/math.sqrt(6)*0.95
    #Emf = 1.2*U1f
    alfa = 2/math.pi
    omega_el = omega_meh * pp
    Ns = math.sqrt(2)*U1f/omega_el/alfa/Bzr/tau_p/l*1000000/Kws1
    Zqs = 2*m*Ns/Qs
    Zqs = math.ceil(Zqs)
    Ns = Zqs*Qs/2/m
    bzr_ = math.sqrt(2)*U1f/omega_el/alfa/Ns/tau_p/l*1000000/Kws1
    izkoristek = 0.95
    Is = Pn/izkoristek/m/cosfi/U1f
    Jcu = 6
    SCus = Is / Jcu
    dcu = math.sqrt(SCus)
    b1s = 4
    b1r = tau_p - 90
    K_bs = b1s/zr_min/(5+b1s/zr_min)

    b1r = 9.48
    K_br = b1r/zr_min/(5+b1r/zr_min)

    bes = K_bs*b1s
    ber = K_br*b1r
    tau_u = obseg_rotorja / Qs
    Kcs = tau_u/(tau_u-bes)
    Kcr = tau_p/(tau_p-ber)
    Kc = Kcr * Kcr

    zr_min = 1/Kcs*4*(1e-7)*tau_p*A_tokovna_obloga/Bzr

    Bds = 1.5
    kfe = 0.97
    Bzr_homogeno = Bzr / Kc
    bds = Bzr_homogeno*tau_u/Bds/kfe
    bds = 9.52
    KCUs = 0.38
    Sus = Zqs*SCus/KCUs
    bs = tau_u - bds
    hs = Sus/bs
    hs = 36.25

    Bsy = 1.25
    fluks_m = alfa*Bzr*(tau_p/1000)*(l/1000)
    hys = fluks_m/2/kfe/(l/1000)/(Bsy/1000)
    Rsn = Rr + zr_min
    obseg_statorja_notranji = math.pi*Rsn*2
    l_reza_utora = 1 
    zagozda = 2
    hs = 36.25
    Rsz = Rsn + hys + hs + l_reza_utora + zagozda 
    Dstroja = 2*Rsz

    #print("Premer zunanji, notranji:",2*Rsz, 2*Rr)

    # Statorski utori
    stator_nodes = []
    stator_nodes.append([Rsn*math.cos(85/180*math.pi), Rsn*math.sin(85/180*math.pi)])
    utorska_odprtina = [b1s/2, math.sqrt(math.pow(Rsn,2) - math.pow(b1s/2,2))]
    stator_nodes.append(utorska_odprtina)
    utorska_odprtina = [utorska_odprtina[0], utorska_odprtina[1]+l_reza_utora]
    stator_nodes.append(utorska_odprtina)
    utorska_odprtina = [utorska_odprtina[0], utorska_odprtina[1]+zagozda]


    #### Višina utora ####
    for i in range(30):
        # Spodnja širina utora - iterativen postopek, dokler ne dosežemo prave širine zoba
        utor = [utorska_odprtina[0], utorska_odprtina[1]]
        spodnja_bds = 1000
        while(spodnja_bds > bds):
            utor[0] = utor[0] + 0.001
            kopija_utora = zrcali_rotiraj(utor[0], utor[1], -360/Qs) 
            spodnja_bds = math.sqrt(math.pow(utor[0] - kopija_utora[0],2) + math.pow(utor[1] - kopija_utora[1],2))


        # Zgornja širina utora - iterativen postopek, dokler ne dosežemo prave širine zoba
        utor_z = [utor[0], utor[1] + hs]
        zgornja_bds = 1000
        while(zgornja_bds > bds):
            utor_z[0] = utor_z[0] +  0.001
            kopija_utora = zrcali_rotiraj(utor_z[0], utor_z[1], -360/Qs)
            zgornja_bds = math.sqrt(math.pow(utor_z[0] - kopija_utora[0],2) + math.pow(utor_z[1] - kopija_utora[1],2))

        hs = Sus / (bs/2 + utor_z[0])
    #print("Visina statorskih zob, povrsina utora:",hs, Sus)


    stator_nodes.append(utor)
    stator_nodes.append(utor_z)
    utor_sredina = [0, utor_z[1]]
    stator_nodes.append(utor_sredina)
    stator_nodes.append([0, Rsn+l_reza_utora])
    stator_nodes.append([0, utor[1]])

    for node in stator_nodes:
        femm.mi_addnode(node[0], node[1])

    femm.mi_drawarc(stator_nodes[0][0], stator_nodes[0][1], stator_nodes[1][0], stator_nodes[1][1], 8, 4/360*2*math.pi*Rsn)
    femm.mi_addsegment(stator_nodes[1][0], stator_nodes[1][1], stator_nodes[2][0], stator_nodes[2][1])
    femm.mi_addsegment(stator_nodes[3][0], stator_nodes[3][1], stator_nodes[2][0], stator_nodes[2][1])
    femm.mi_addsegment(stator_nodes[3][0], stator_nodes[3][1], stator_nodes[4][0], stator_nodes[4][1])
    femm.mi_addsegment(stator_nodes[5][0], stator_nodes[5][1], stator_nodes[4][0], stator_nodes[4][1])
    femm.mi_addsegment(stator_nodes[3][0], stator_nodes[3][1], stator_nodes[7][0], stator_nodes[7][1])


    ###### NAREDI RADIJ ZOBA
    femm.mi_createradius(stator_nodes[4][0], stator_nodes[4][1], 1)

    femm.mi_selectcircle(0, 0,Rsz,1)
    femm.mi_mirror2(0,0,Rsz*math.cos(math.pi/2), Rsz*math.sin(math.pi/2),1)
    femm.mi_selectcircle(0, 0,Rsz,3)
    femm.mi_mirror2(0,0,Rsz*math.cos(math.pi/2), Rsz*math.sin(math.pi/2),3)

    femm.mi_selectcircle(0, 0,Rsz,1)
    femm.mi_copyrotate2(0,0, 360/Qs, Qs-1, 1)
    femm.mi_selectcircle(0, 0,Rsz,3)
    femm.mi_copyrotate2(0,0, -360/Qs, Qs-1, 3)

    Rsz = Rsn + hys + hs + l_reza_utora + zagozda 
    #print(Rsz, Rr/Rsz)
    # Stator
    femm.mi_addnode(-Rsz, 0)
    femm.mi_addnode(Rsz, 0)
    femm.mi_drawarc(Rsz, 0, -Rsz, 0, 180, 1/2*2*math.pi*Rsz)
    femm.mi_drawarc(-Rsz, 0, Rsz, 0, 180, 1/2*2*math.pi*Rsz)


    ### Rotor
    bdr = 45
    tau_cevlja_r = 90
    b1r = tau_p - tau_cevlja_r

    ################################################
    # Sinusna porazdelitev zoba
    tocke_a = []
    cevelj_rad = (tau_cevlja_r/2 / obseg_rotorja) *2 * math.pi
    krivulja_rad = 80/90*math.pi/2
    rad_a = cevelj_rad/(tau_cevlja_r/2)
    krivulja_a_rad = krivulja_rad/(tau_cevlja_r/2)

    for a in range(int(tau_cevlja_r/2)+1):
        zr_a = zr_min / math.cos(a*krivulja_a_rad) 
        R_zoba = Rsn - zr_a
        tocke_a.append([R_zoba*math.cos(math.pi/2 - a*rad_a), R_zoba*math.sin(math.pi/2 - a*rad_a)])

    tocke_a.append([R_zoba*math.cos(math.pi/2 - a*rad_a), Dr/2 - 16])

    for tocka in tocke_a:
        femm.mi_addnode(tocka[0], tocka[1])

    for i in range(len(tocke_a)-1):
        femm.mi_addsegment(tocke_a[i][0],tocke_a[i][1], tocke_a[i+1][0], tocke_a[i+1][1])

    ################################################

    rotor_nodes = []
    rotor_nodes.append([tocka[0],tocka[1]])
    rotor_nodes.append([bdr/2,tocka[1]])
    rotor_nodes.append([bdr/2,bdr*math.cos(math.pi/6)])

    for node in rotor_nodes:
        femm.mi_addnode(node[0], node[1])

    femm.mi_addsegment(rotor_nodes[0][0],rotor_nodes[0][1], rotor_nodes[1][0], rotor_nodes[1][1])
    femm.mi_addsegment(rotor_nodes[0][0],rotor_nodes[0][1], tocka[0], tocka[1])
    femm.mi_addsegment(rotor_nodes[2][0],rotor_nodes[2][1], rotor_nodes[1][0], rotor_nodes[1][1])
    femm.mi_addsegment(rotor_nodes[2][0],rotor_nodes[2][1], rotor_nodes[0][0], rotor_nodes[0][1])

    femm.mi_selectcircle(0, 0,Rr+1,4)
    femm.mi_mirror2(0,0,0, Rr, 1)

    femm.mi_selectcircle(0, 0,Rr+1,4)
    femm.mi_copyrotate2(0,0, 60, 5, 1)


    # Rotorski jarm
    hyr = 22.27
    for i in range(0,6):
        femm.mi_addnode((bdr - hyr)*math.sin(i*math.pi/3+math.pi/6),(bdr - hyr)*math.cos(i*math.pi/3 + math.pi/6))


def koncaj_geometrijo():
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    femm.mi_addcircprop("A", 0, 1)
    femm.mi_addcircprop("B", 0, 1)
    femm.mi_addcircprop("C", 0, 1)


    navitja = ["A", "c", "c", "B", "B", "a", "a", "C", "C", "b", "b", "A", "A", "c", "c", "B", "B", "a", "a", "C", "C", "b", "b", "A", "A","c", "c", "B", "B", "a", "a", "C", "C", "b", "b", "A"]
    #print(len(navitja))
    #femm.mi_getmaterial('M330-35A')

    # Load B-H data from a file
    bh_points = []
    with open("BH.txt", "r") as f:
        for line in f:
            if line.strip():
                B, H = map(float, line.strip().split())
                bh_points.append((B, H))

    # Add a placeholder material with default properties
    femm.mi_addmaterial("M330-35A", 1, 1, 0, 0, 0, 0)

    # Add B-H curve to the material
    for B, H in bh_points:
        femm.mi_addbhpoint("M330-35A", B, H)


    femm.mi_getmaterial('Air')
    femm.mi_addmaterial('Copper', 1, 1, 0)

    block_label = [10, 150]
    femm.mi_addblocklabel(block_label[0], block_label[1])
    femm.mi_selectcircle(block_label[0], block_label[1], 5, 2)
    femm.mi_setblockprop('M330-35A', 1, 0, "<None>", 0, 0, 1)
    femm.mi_clearselected()

    block_label = [0, 0]
    femm.mi_addblocklabel(block_label[0], block_label[1])
    femm.mi_selectcircle(block_label[0], block_label[1], 5, 2)
    femm.mi_setblockprop('M330-35A', 1, 0, "<None>", 0, 1, 1)
    femm.mi_clearselected()

    block_label = [0, 97]
    femm.mi_addblocklabel(block_label[0], block_label[1])
    femm.mi_selectcircle(block_label[0], block_label[1], 5, 2)
    femm.mi_setblockprop('Air', 1, 0, "<None>", 0, 0, 1)
    femm.mi_clearselected()


    for i in range(6):
        block_label = [70*math.cos(3*math.pi/8 + i*60*math.pi/180), 70*math.sin(3*math.pi/8 + i*60*math.pi/180)]
        femm.mi_addblocklabel(block_label[0], block_label[1])
        femm.mi_selectcircle(block_label[0], block_label[1], 5, 2)
        if(i%2 == 0):
            ovoji = 20
        else:
            ovoji= -20
        femm.mi_setblockprop('Copper', 1, 0, "DC", 0, 1, ovoji)
        femm.mi_clearselected()

        block_label = [70*math.cos(5*math.pi/8 + i*60*math.pi/180), 70*math.sin(5*math.pi/8 + i*60*math.pi/180)]
        femm.mi_addblocklabel(block_label[0], block_label[1])
        femm.mi_selectcircle(block_label[0], block_label[1], 5, 2)
        if(i%2 == 0):
            ovoji = -20
        else:
            ovoji= 20
        femm.mi_setblockprop('Copper', 1, 0, "DC", 0, 1, ovoji)
        femm.mi_clearselected()


    for i in range(36):
        block_label = [120*math.cos(math.pi/2 + i*10*math.pi/180), 120*math.sin(math.pi/2 + i*10*math.pi/180)]
        femm.mi_addblocklabel(block_label[0], block_label[1])
        femm.mi_selectcircle(block_label[0], block_label[1], 5, 2)
        ovoji = 5
        navitje = navitja[i]
        if(navitja[i] == "a" or navitja[i] == "b" or navitja[i] == "c"):
            ovoji = -ovoji
            if(navitja[i] == "a"):
                navitje = "A"
            elif(navitja[i] == "b"):
                navitje = "B"
            elif(navitja[i] == "c"):
                navitje = "C"
        femm.mi_setblockprop('Copper', 1, 0, navitje, 0, 0, ovoji)
        #femm.mi_setblockprop('Copper', 1, 0, "<None>", 0, 0, 1)
        femm.mi_clearselected()


    femm.mi_selectarcsegment(0, 163)
    femm.mi_selectarcsegment(0, -163)
    femm.mi_setarcsegmentprop(1, "A=0", 0, 0)
    femm.mi_clearselected()

    # Vrtenje
    femm.mi_selectcircle(0, 0, 96, 1)
    femm.mi_setsegmentprop("", 0, 1, 0, 1)
    

    # Nastavi kolesni kot na 0 stopinj
    femm.mi_selectgroup(1)
    femm.mi_moverotate(0,0, -35)
    femm.mi_clearselected()

