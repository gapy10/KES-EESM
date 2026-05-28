%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Rezultati

load("result") %inicializacija spremenljivk za preizkuse stroja
AdP=result(1);
Bd=result(2);
kE=result(3);
ni=result(4);
cosfi=result(5);
kLok=result(6);
Ad=25000*AdP;
design12(Ad,Bd,kE,ni,cosfi,kLok)

M=NavornaKarak()
RelekNavor()

SamodrzniNav()

IndNap()
KPT()

load('mere.mat')
load('konst.mat')
Lm=(m/pi^2)*mi0*Tp*l*kfe*(kw*N)^2/(2*Pp*d)
Ldq=DQind()

Breza()
Izgube(M)

V=l*(Dr/2)^2*pi
l/D