function f=SamodrzniNav() %izračun samodržnega novra novra

%clc;clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("mere.mat")
load('konst.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rpc=Dr/2;
Rr=Dr/2;

IDC=Idc;
IS=Is;
ar=ar;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Samodržni navor

t=0; %-4*pi/Q;  %-(2*pi/Pp)-(Q/pi)
Idc=IDC;
Is=0;
Is1=sqrt(2)*Is*sin(t)/a;
Is2=sqrt(2)*Is*sin(t+2*pi/3)/a;
Is3=sqrt(2)*Is*sin(t+4*pi/3)/a;

ar=ar;

openfemm;
newdocument(0); % 0 for magnetics problem
smartmesh(0);
rezolucija=1e-8; % original 1e-8
minangle=30;
mi_probdef(0, 'meters', 'planar', rezolucija, l, minangle);

izris2(Pp,Q,D,Dr,d,Rs,Rrz,hpc,bdr,ar,br,Ndc,hsz,hs,hys,b1,b2,bdss,bds,fizob,Zq,Idc,m,Is1,Is2,Is3,alfa,Rpc,angle1)
mi_selectcircle(0,0,Rr,4); %group 1 je cel rotor
mi_setgroup(1);
mi_clearselected();

%Izracun navorne karakteristike
fiel=2; %inkrement vrtenja (električnih stopinje)
kotkarakteristike=180; %željen pregled kota (električne stopinje)
T=Navorkota(fiel,Pp,kotkarakteristike);
F=0:fiel:kotkarakteristike;

figure()
hold on
grid on
title('Samodržni Navor')
plot(F,T)
xlabel('Električni kot [°]')
ylabel('Samodržni Navor [Nm]')
xlim([0,kotkarakteristike])

closefemm;











