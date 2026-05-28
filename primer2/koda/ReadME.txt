%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pozdravljeni, tukaj so razložene vse funkcije v programu in kako med seboj delujejo
 
V grobem so trije tipi funkcij-glavne, podporne in za analizo. Podpore so naprimer funkcija za izris stroja v FEEM, za izračun navora, rotacijo rotorja,....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delitev in opis funkcij:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Glavne:-

Zahtevke: Tukaj se definira zahtevane parametre stroja, konstante ter izhodišče parametre za kasnejšo optimizacijo. Naredi datoteki konst.mat in result.mat

OTP3: Optimizacijska funkcija, ki vključuje tudi kriterijsko funkcijo. Po optimizaciji shrani končne podatke v result.mat in fval, Kriterihska funkcija kliče design12, izrisstroja in izgube.

design12: To je funckcija za analitični izračun stroja

izrisstroja: Izriše navorne karakteristiko stroja in preko FFT pridobi vzbujalni in samodržni navor. Ta funkcija tudi najde kolsni kot pri katerm doseže stroj največji navor in ga shrani v kotMAX.mat

Izgube: Kliče kotMAX.am in postavi stroj v najbolj obremenjeno stanje, ter izračuna izgube. Pri izračunu upošteva skin efekt, inducirane tokove v vodniki zaradi mag.polja v utorih, iz FEMM uvozi magnetno stanje stroja in izračuna izgube v železu,.... Pri izračunu se zanemarijo izgube  rotorskem železu

Rezultati: Ta funckija kliče vse funkcije za analozo stroja in se avtomatsko izvede po končani optimizaciji

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Podporne:-

Izris2: Izriše stroj v FEMM

Navor: Izračuna navor preko vgrajenega orodja v FEMM

Navorkota: Izračuna navorno karakteristiko stroja

zavrtiR: zavrti rotor za nek mehanski kot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Za Analizo:-

NavornaKarak: Izriše surovo navorno karakteristiko, karkaterisko samo vzbujalnega in reluktančnega anvora ter FFT speketer navora

RelekNavor: Izriše karakteristiko reluktančnega navora

SamodržniNav: Izriše karakteristiko samodržnega navora

IndNap: Izračuna magentne slopitve, inducirano napetost in FFT spekter inducirane napetosti

KPT: Izračuna in izriše karakteristiko prostega teka

DQind: Izračun linearne in nasičene DQ induktivnosti

Breza: Izriše magnento gototo v reži

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Če želite le verificirate rezultate lahko poženete le Rezultati.m, ki bo izrisal vse grafe in podatke (izračuni lahko trajajo nekaj časa saj koda ni optimizirana za hitrost,odvisno od računalnika 30min-2h) 




