function f=IndNap() %izračun induciane napetosti in FFT ind nap

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
%Inducirana napetos

Idc=IDC;
Is1=0;
Is2=0;
Is3=0;

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

theta=360; %360; % gledan kot (el)
theta1=180; % št izračunov na električno periodo
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
E1(i)=(psi1(i+1)-psi1(i))/(dt);
E2(i)=(psi2(i+1)-psi2(i))/(dt);
E3(i)=(psi3(i+1)-psi3(i))/(dt);

i=i+1;
end

figure()
grid on
hold on
title('Inducirana napetost')
xlabel('Električni kot [°]')
ylabel('Inducirana napetost [V]')
plot(fi(1:end-1)*2,E1)
plot(fi(1:end-1)*2,E2)
plot(fi(1:end-1)*2,E3)
xlim([5,theta])
ylim([-700,700])

figure()
grid on
hold on
title('Magnetni sklep')
xlabel('Električni kot [°]')
ylabel('Magnetna sklopitev [Wb]')
plot(fi(1:end)*2,psi1)
plot(fi(1:end)*2,psi2)
plot(fi(1:end)*2,psi3)
xlim([5,theta])

E1=[E1,E1,E1,E1,E1,E1,E1,E1,E1,E1];

dt=dt(1);
timedom=0:dt:length(E1)*dt;
L=length(timedom);

Y=fft(E1);

warning('off','all');
% Compute the two-sided spectrum P2 and single-sided spectrum P1
P2 = abs(Y/L);       % Normalize the amplitude
P1 = P2(1:L/2+1);    % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double amplitudes (except DC and Nyquist)

% Frequency domain vector
f1 = 1/dt*(0:(L/2))/L;

warning('on','all');

% Plot
figure;
hold on
grid on
plot(f1, P1)
title('Harmonske komponente inducirane napetosti')
xlabel('Frekvenca [Hz]')
ylabel('Amplitude harmonskih komponent [V]')

% Plot
f2=f1/f1(11);
P2=P1/P1(11)*100;

figure;
hold on
grid on
plot(f2, P2)
title('Harmonske komponente inducirane napetosti')
xlabel('Harmonske komponente');
ylabel('Amplitude harmonskih komponent [%]');
xlim([0,50])

%closefemm;