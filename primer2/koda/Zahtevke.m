clc;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vhodni podatki

Udc=1000; % V
Pm=150000; % kW kornr point
RPMk=7000; %RPM
RPMm=14000; %RPM
RPMmehanskimax=16000;

Umf=Udc/sqrt(2);  %3*sqrt(6)/pi^2; % V
Uf=Umf/sqrt(3); %V
Ufa=Uf*sqrt(2); % V

Mn=Pm/(2*pi*(RPMk/60)); %Nm
knav=1.05; %dodaten faktor za navor, da ga lahko kasneje poševljamo
Mn=Mn*knav;

m=3;
Pp=2;
kut=2;
Q=3*2*Pp*kut; %Utori
a=2;

fk=Pp*RPMk/60;
fm=Pp*RPMm/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Konstante / parametri storja

mi0=1.256637*10^-6; % mi 0 

b1=0.003; % razdalja med zobi 2mm-4mm
b2=0.0025; % visina zoba
fizob=30*pi/180; %kot med povezovalnim delom zgornejga in spodnjega dela zoba
kfe=0.95; %polnilni faktor pločevine

alfa=0.63662;
ep=4*10^-7; % konstanta  za zračno režo za SS

kw=0.96;

Jr=5000000; % A/m2  tokovna gostota rotor
Js=6000000; %A/m2 tokovna gostota stator

kvar=RPMmehanskimax/RPMm; %varnostni faktor vrtilne hitrosti
kv=2; % faktor varnosti moči gredi

kpr=0.55; % faktor polenja vzb navitja
kp=0.45; % koeficient polenja satora

diz=0.0002; % debelina izolacij žic
dizz=0.0005; %debelina izolacije rotorja
dl=0.00035; %debelina lamele

% Bds=1.3; % B v zobu statorja
% Bdr=1.4; % B v stebru pola rotorja
% Bdj=1.35; % B v jarmu

Bds=1.25; % B v zobu statorja
Bdr=1.35; % B v stebru pola rotorja
Bdj=1.3; % B v jarmu

kg=0.4; %delež pola, ki je polov cevel 40%

save('konst.mat',"Udc","Pm","RPMk","RPMm","RPMmehanskimax","Umf","Uf","Ufa","Mn","knav","m","Pp","kut","Q","a","fk","fm","mi0","b1","b2","fizob","kfe","kfe","alfa","ep","kw","Jr","Js","kvar","kv","kp","kpr","diz","dizz","dl","Bds","Bdr","Bdj","kg")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Izhhodiščni parametri za optimizacijo

% AdP=1.1250 ;
% Bd=0.9625;
% kE=1.2 ;
% ni=0.9625 ;
% cosfi=0.9313;
% kLok=-0.0156;
% Ad=25000*AdP;

AdP=1.1250 ;
Bd=0.9625;
kE=1.1750 ;
ni=0.9625 ;
cosfi=0.9313;
kLok=-0.0156;
Ad=25000*AdP;

% result=[AdP,Bd,kE,ni,cosfi,kLok];
% 
% save("result.mat","result")

%1.1250    0.9625    1.1750    0.9625    0.9313   -0.0156

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




