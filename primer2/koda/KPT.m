function f=KPT() %karakteristika prostega teka
%večkratni izračun ind nap

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Karakteristika prostega teka

inkTOK=IDC/5;

Idc=0:inkTOK:IDC*2;
Is1=0;
Is2=0;
Is3=0;

% Idc(end+1)=IDC*5;
% Idc(end+1)=IDC*10;


ii=1;
while ii<=length(Idc)

openfemm;
newdocument(0); % 0 for magnetics problem
smartmesh(0);
rezolucija=1e-8; % original 1e-8
minangle=30;
mi_probdef(0, 'meters', 'planar', rezolucija, l, minangle);

izris2(Pp,Q,D,Dr,d,Rs,Rrz,hpc,bdr,ar,br,Ndc,hsz,hs,hys,b1,b2,bdss,bds,fizob,Zq,Idc(ii),m,Is1,Is2,Is3,alfa,Rpc,angle1)
mi_selectcircle(0,0,Rr,4); %group 1 je cel rotor
mi_setgroup(1);
mi_clearselected();

mi_saveas('design1.FEM');
mi_createmesh;
mi_analyze(1);

mi_loadsolution;

A1=mo_getcircuitproperties('NavitjeAC1');
psi1=A1(3)/a;
A2=mo_getcircuitproperties('NavitjeAC2');
psi2=A2(3)/a;
A3=mo_getcircuitproperties('NavitjeAC3');
psi3=A3(3)/a;

theta=180; %360; % gledan kot (el)
theta1=90; % št izračunov na električno periodo
ink=(theta/theta1)/Pp; %vsako meansko stopinje
fi=0:ink:theta/Pp; % ° meh

i=1;
while i<length(fi)

zavrtiROT(ink);

mi_saveas('design1.FEM');
mi_createmesh;
mi_analyze(1);

mi_loadsolution;

A1=mo_getcircuitproperties('NavitjeAC1');
psi1(i+1)=A1(3)/a;
A2=mo_getcircuitproperties('NavitjeAC2');
psi2(i+1)=A2(3)/a;
A3=mo_getcircuitproperties('NavitjeAC3');
psi3(i+1)=A3(3)/a;

dt=(fi(i+1)-fi(i))*60/RPMk/360;
E1(i)=(psi1(i+1)-psi1(i))/dt;

i=i+1;
end


E1=[E1,fliplr(E1)];
fi=0:ink:(length(E1)-1)*ink;

% figure()
% grid on
% hold on
% title('Inducirana napetost')
% xlabel('Električni kot [°]')
% ylabel('Inducirana napetost [V]')
% plot(fi(1:end),E1)
% xlim([5,theta])


E1=[E1,E1,E1,E1,E1,E1,E1,E1,E1,E1];
fi=0:ink:(length(E1)-1)*ink;


dt=dt(1);
timedom=0:dt:length(E1)*dt;
L=length(timedom);

warning('off','all');

Y=fft(E1);

% Compute the two-sided spectrum P2 and single-sided spectrum P1
P2 = abs(Y/L);       % Normalize the amplitude
P1 = P2(1:L/2+1);    % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double amplitudes (except DC and Nyquist)

% Frequency domain vector
f1 = 1/dt*(0:(L/2))/L;

warning('on','all');

% figure;
% plot(f1, P1);
% title('Harmonske komponente inducirane napetosti');
% xlabel('Frekvenca [Hz]');
% ylabel('Amplitude harmonskih komponent [V]');
% %xlim([0,1000])

%[peaks,loc]=findpeaks(P1,f1,'SortStr','descend');

Uind(ii)=max(P1)/sqrt(2);

closefemm;

ii=ii+1;
clear E1
end

figure()
grid on
hold on
title('Karakteristika prostega teka')
xlabel('Vzbujalni tok [A]')
ylabel('Inducirana napetost [V]')
plot(Idc,Uind*1.045)
xlim([0,14])


