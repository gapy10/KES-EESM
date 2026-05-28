%function f=NavornaKarak()
%izris navorne karaktersitike in FFT navora
clc;clear all
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
%Celotna Navorna karakteristike

t=0; %-4*pi/Q;  %-(2*pi/Pp)-(Q/pi)
Idc=IDC;
Is=IS;
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
kotkarakteristike=360; %360; %željen pregled kota (električne stopinje)
T=Navorkota(fiel,Pp,kotkarakteristike);
F=0:fiel:11*kotkarakteristike+10*fiel;

T=[T,T,T,T,T,T,T,T,T,T,T];

figure()
hold on
grid on
title('Navor')
plot(F,T)
xlabel('električni kot [°]')
ylabel('Navor [Nm]')
xlim([0,kotkarakteristike])

dt=fiel*60/RPMk/360/Pp;
timedom=0:dt:length(F)*dt;
L=length(timedom);

Y=fft(T);

% Compute the two-sided spectrum P2 and single-sided spectrum P1
P2 = abs(Y/L);       % Normalize the amplitude
P1 = P2(1:L/2+1);    % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double amplitudes (except DC and Nyquist)

% Frequency domain vector
f = 1/dt*(0:(L/2))/L;

% Plot
figure;
hold on
grid on
plot(f, P1);
title('Harmonske komponente Navorne karakteristike');
xlabel('Frekvenca [Hz]');
ylabel('Amplitude harmonskih komponent [Nm]');
%xlim([0,1000])

[peaks,loc]=findpeaks(P1,f,'SortStr','descend');
t=0:pi/180:pi;

figure()   
hold on
grid on
plot(t*180/pi,peaks(1)*sin(t)+peaks(2)*sin(2*t))
xlabel('električni kot [°]')
ylabel('Navor [Nm]')

f2=f/f(12);
P2=P1/P1(12)*100;

figure;
hold on
grid on
plot(f2, P2)
title('Harmonske komponente Navorne karakteristike');
xlabel('Harmonske komponente');
ylabel('Amplitude harmonskih komponent [%]');
xlim([0,50])

%closefemm;

f=max(peaks(1)*sin(t)+peaks(2)*sin(2*t));

closefemm;

