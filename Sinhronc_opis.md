# Sinhronski elektromotor z vzbujalnim navitjem

To je moj zaključni projekt za predmet Konstruiranje električnih strojev na Fakulteti za Elektrotehniko Univerze v Ljubljani. Uporaba umetne inteligence je dovoljena.

Narediti je potrebno MATLAB ali Python program, ki bo analitično izračunal dimenzije sinhronskega elektromotorja z navitim rotorjem glede na vhodne parametra, optimiziral stroj in ga simuliral v programu FEMM Finite Element Method Magnetics in analiziral rezultate. Priporočam python, da ga boš lažje zagnal in analiziral rezultate, če pa znaš zagnati MATLAB in predbrati razultate simulacije v matlabu pa naredi program v matlabu. V kodi naj bodo komentarji, tako da bo kodo lahko razumel kdorkoli, ki se spozna na električne stroje. 

Da ti bo lažje sem priložil sem tri mape s podobnimi projekti pri tem predmetu:
PRIMER1: tukaj je pdf poročilo, ki opisujeta potek analitičnega izračuna stroja in ostalih delov programa in python koda za 6 polni sinhronski stroj z vzbujalnim navitjem od lanskega leta.
PRIMER2: tukaj je pdf poročilo, ki opisujeta potek analitičnega izračuna stroja in ostalih delov programa in matlab koda za 4 polni sinhronski stroj z vzbujalnim navitjem od lanskega leta.
PRIMER3: tukaj je pdf poročilo, ki opisujeta potek analitičnega izračuna stroja in ostalih delov programa za 8 polni sinhronski stroj z vzbujalnim navitjem od letos.

podrobno preberi vse 3 pdf dokumente in preglej matlab in python kodo.

NUJNO: Najbolj pomebna je python koda v primer1. Njegov motor je prav tako kot moj šestpolni z vzbujalnim navitjem. Nujno zaženi python program znotraj mape primer1 in analiziraj delovanje python kode. Sploh tisti del, ki v programu FEMM izriše cel motor. Natančno si poglej kako je motor narisan. Moj mora biti praktično dobesedno enako narisan, samo malenkost drugačne dimenzije.

NUJNO: preberi vsa 3 pdf poročila v mapah primer1 primer2 in primer3, da boš vedel na kakšen način je treba analitično izračunati vse parametre sinhronskega motorja.



## Vhodni podatki

Število faz: 3
Število polov: 6 (torej število polovih parov: 3)
Napetost 200 V (fazna RMS vrednost napetosti)
P<sub>c</sub>: 50 kW (nazivna moč v vogalni točki)
n<sub>c</sub>: 7000 vrt/min (mehanska hitrost vrtenja v nazivni/vogalni točki)
N<sub>max</sub>: 14000 vrt/min (najvišja obratovalna hitrost)
Material: M330-35A (lamelirana elektropločevina)



### Dodatna navodila


Sinhronski motor z vzbujalnim navitjem naj bo 3 fazni. Za vrednost q (število utorov na pol in fazo) izbirajte med vrednostmi 1, 1.5, 2, 2.5, 3. To vam bo definiralo število utorov na statorju ter faktor navitja. Korelacije q (število utorov statorja, najvišji faktor navitja k<sub>w1</sub>) povzemite iz spletne strani EMETOR. Parameter q je lahko optimizacijska spremenljivka, ni pa nujno.

Najprej je potrebno analitično izračunati začetno geometrijo stroja, ki mora biti že v štartu čim boljša. Preberi si oba pdf poročila od mojih kolegov, kjer je na začetku obrazloženo kako so izračunali stroj. Preberi si tudi obsežno knjigo "Design of rotating electrical machines SECOND EDITION.pdf"., kjer je vse to razloženo v potankosti. Upoštevaj, da mora biti narisan model v programu FEMM točno tak kot je v python kodi od kolega. Rotorski zobje za vsak pol morajo biti na koncu široki in se zaključiti zaokroženo po sinusu, tako da je na sredini najmanjša zračna reža, na konceh pola pa širša. To da lepo sinusno inducirano napetost. Vzbujalno rotorsko navitje mora biti na rotorskem zobu za vsak pol posebej, ločeno od drugega pola. Tako kot je na sliki "stroj.jpg" Statorski zobje morajo biti tudi na koncu pri zračni reži širši, da se raztegnejo malo v sosednji vdolbini za statorska navitja. Tako kot je na sliki statorski utor "utor.jpg". Pazi na to, da že bo analitičen izračun tak, da ne pridejo rezultati preveč ploščati motorji ki imajo velik radij in majhno dolžino, ker to ni dobro za tehnološko izvedbo.

![Zajeta slika](C:\Users\gaspe\OneDrive\Namizje\KES v7\stroj.JPG)

![Zajeta slika1](C:\Users\gaspe\OneDrive\Namizje\KES v7\utor.JPG)

B-H diagram materiala M330-35A je v pdf dokumentu. Matlab koda za izračun izgub v železu je b tekstovni datoteki izgube_fem.txt. Če boš uporabil python jo konvertiraj v python in jo uporabi v izračunu izgub. Če boš uporbil MATLAB pa jo obdrži in vključi v izračun izgub.

Pri optimizaciji je ključna maksimizacija izkosristka in minimizacija volumna. Pazi pa na to, da bo motor še vedno možno fizično izdelati, torej, da se da namestiti vsa navitja v svoje utore.

Program na koncu izbere še nekaj optimiziranih modelov in jih "D izriše v programu za simulacijo FEMM Finite Element Method Magnetics. Navodila za uporabo programam FEMM si preberi v pdfju "manual.pdf" in za interakcijo s programom preko kode si preberi pdf dokument "octavefemm.pdf". Tukaj si res pomagaj z že priloženo kodo od kolega, da jo dobro pregledaš da boš znal pravilno upravljati FEEM skozi kodo.

Za izdelavo seminarske naloge lahko uporabite vsa razpoložljiva programska orodja kot tudi umetno inteligenco (UI). V primeru, da boste za tvorjenje matematične kode uporabili UI morate zahtevati, da bo kota napisana z dodanimi komentarji h kodi kot tudi z referencami na znanstvene publikacije (iz kje je UI črpala vsebino za predlagano kodo). Koda naj bo pregledna tako, da jo boste razumeli in mi jo znali razložiti/interpretirati v fazi zagovora. Če si se odločil za MATLAB potem naj bo koda v obliki »Matlab Live Script«, če pa bo python potem pa python.

NUJNO: Na koncu napiši še podrobno obsežno poročilo. Poročila ne rabi generirati program. Program naj naredi analitičen izračun, optimizacijo z genetskimi algoritmi, pareto fronto, simulacijo 5 motorjev v programu FEMM in naj vrne rezultate simulacije. V poročilu pa ti kot umetna inteligenca predstavi celoten potek delovanja programa. Vključi vse enačbe in postopek analitičnega izračuna. Dodaj tudi vire na znanstvene reference iz kje so dobljene enačbe oziroma priporočila za vrednosti parametrov. Napiši tudi tabelo vseh spremenljivk njihove kratice in polno ime spremenljivke. Razloži cel postopek optimizacije, kako deluje genetski algoritem, katere so optimizacijske spremenljivke. Razloži kako deluje risanje motorja in analiza motorja znotraj programa FEMM. Na koncu pa še za vseh 5 analiziranih motorjev opiši in komentiraj. Poročilo je lahko v markdown lahko pa je v latex, obvezno da ima notri vse enačbe za vse izračune. Obvezno naj bo v poročilu pojasnjeno čisto vse kar dela program, da če me profesor kaj vpraša glede programa bom vedel kaj program dela.



#### Podrobne naloge - obvezno:

###### A:

Glede na zahtevane lastnosti sinhronskega stroja z električnim vzbujanjem analitično preračunajte vrednosti geometrijskih spremenljivk (celotno geometrijo). Nasvet: v statorskem navitju privzemite maksimalno gostoto toka 10 A/mm2 , v rotorskem vzbujalnem navitju pa 5 A/mm2

###### B:

Analitično določite funkcijsko odvisnost amplitudne vrednosti gostote magnetnega pretoka v zračni reži 𝐵𝛿 od vzbujalnega toka in od vzbujalne magnetne napetost 𝜃m.

###### C:

S polinomsko aproksimacijo določite izgube v železu v funkciji amplitude gostote magnetnega pretoka in frekvence za neorientirano elektro pločevino tipa M230-35A (Sura Cogent podatkovni list je v priponki, kot tudi Matlab koda za polinomsko aproksimacijo izgub v odvisnosti od amplitude gostote magnetnega pretoka in frekvence).

###### D:

V vsakem delu statorskega paketa (zob, jarem) izračunajte izgube v železu Pfe ter nato Pfe v celotnem statorskem paketu. Pfe izračunajte v dani obratovalni točki (pozor na vrednost B-ja in frekvence v obratovalni točki).

###### E:

Izračunajte izgube v navitju statorja in v vzbujalnem navitju na rotorju ter skupne izgube PCu. Dolžino glav navitji prevzemite enako kot je dolžina statorskega/rotorskega feromagnetnega paketa. PCu izračunajte v dani obratovalni točki.

###### F:

V dani obratovalni točki izračunajte izkoristek motorja.

###### G:

Analitičen preračun vključite v optimizacijski postopek na podlagi genetskih algoritmov (GA). 

###### H:

Ciljni funkciji sta minimizacija izgub (maksimizacija izkoristka) in minimizacija volumna aktivnega dela motorja (kot aktivni del je mišljen premer statorja in dolžina paketa). Vse preračune diktirane s strani genetskih algoritmov predstavite v obliki 2D Pareto optimalnih rešitev. Iz množice rešitev določite Pareto fronto.

###### I:

Za dobljeni analitičen model pod točko A zgradite v programskem okolju FEMM parametrično zasnovan numeričen model (s pomočjo LUA skript-a).

###### J:

Modelu stroja v programskem okolju FEMM definirajte materiale, vstavite statorsko in rotorsko navitje ter ga pravilno omrežite s končnimi elementi (metoda končnih elementov).

###### K:

 Iz Pareto fronte po lastni izbiri izberite 5 Pareto optimalnih rešitev (rezultat dela pod točko H). 

###### L:

Za vsako od petih Pareto optimalnih rešitev opravite preračun v okolju FEMM in sicer tako da, rotorsko vzbujalno navitje napajajte z analitično izračunanim tokom (analitično preračunan s pomočjo GA dobljeno pod točko H in K ) in pri tem vrtite rotor s kotnim inkrementom tako, da dosežete skupen premik za en polov par. Tok v statorskih navitjih je nič, saj stroj deluje v prostem teku v generatorskem režimu. Na podlagi spremembe magnetnega pretoka (glede na inkremente rotacije rotorja) skozi fazno navitje izračunajte inducirano napetost pri hitrosti v vogalni točki (definirana kot vhodni podatek). Opozorilo: Inducirana napetost naj bi bila U<sub>ind</sub> <= (U<sub>grid</sub>*0.95)  (primerno za motorsko delovanje).

###### M:

Za vsako od petih Pareto optimalnih rešitev opravite preračun v okolju FEMM in sicer tako da, rotorsko navitje napajajte z vrednost vzbujalnega toka (analitično preračunan s pomočjo GA dobljeno pod točko H in K ) ter statorsko navitje napajajte z analitično dobljenim statorskim tokom (analitično preračunan s pomočjo GA dobljeno pod točko H in K ) in sicer tako da je: Ia -> amplitudna vrednost nazivnega toka, Ib= - Ia/2 in Ic= -Ia/2 in pri tem vrtite rotor s kotnim inkrementom tako, da dosežete skupen premik za en polov par. Za vsako točko premika izračunajte vrednost navora. Rezultate prikažite v grafu: Navor proti kolesnemu kotu (kotu premika). Dobljeno karakteristiko analizirajte s Fourierovo transformacijo (osnovna in višje harmonske komponente predvsem 2. harmonska komponenta).

###### N:

Kritično ovrednotite in komentirajte pet izbranih Pareto optimalnih rešitev.