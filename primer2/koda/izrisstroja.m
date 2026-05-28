function f=izrisstroja()

load("mere.mat")
load('konst.mat')

openfemm;
newdocument(0); % 0 for magnetics problem
smartmesh(0);
rezolucija=1e-8; % original 1e-8
minangle=30;
mi_probdef(0, 'meters', 'planar', rezolucija, l, minangle);

Rpc=Dr/2;
Rr=Dr/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Navorna karakteristike

t=0; % z spreminjanjem št. utotorv je potrebno spremniti magnetno polje stroja, da se karakteristika začne pri 0 Nm
Is1=sqrt(2)*Is*sin(t)/a;
Is2=sqrt(2)*Is*sin(t+2*pi/3)/a;
Is3=sqrt(2)*Is*sin(t+4*pi/3)/a;

ar=ar;

%Izris stroja
izris2(Pp,Q,D,Dr,d,Rs,Rrz,hpc,bdr,ar,br,Ndc,hsz,hs,hys,b1,b2,bdss,bds,fizob,Zq,Idc,m,Is1,Is2,Is3,alfa,Rpc,angle1)
%Zagon izracuna
mi_selectcircle(0,0,Rr,4); %group 1 je cel rotor
mi_setgroup(1);
mi_clearselected();

%Izracun navorne karakteristike
fiel=2; %inkrement vrtenja električnih stopinj
kotkarakteristike=180; %360; (el stopinje)
T=Navorkota(fiel,Pp,kotkarakteristike);

T=[T,-fliplr(T)]; %izračunamo samo polovico navrone karakteristike in prezrcalimo za drugo (180° el.)
T=[T,T,T,T,T,T,T,T,T,T]; %karakteristiko kooepramo da pridobimo rezolucijo v FFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ločimo valovitost navorne karakteristike in seštejemo samodržni in vzbujalni navor

fi=0:fiel:10*(2*kotkarakteristike+fiel)+9*fiel; % v elek °
firot=fi/Pp;  % v °

% figure()    
% grid on
% plot(firot,T)
% xlabel('mehanski kot [°]')
% ylabel('Navor [Nm]')
% xlim([0,180])

warning('off','all');

dt=fiel*60/RPMk/360/Pp;
timedom=0:dt:length(firot)*dt;
L=length(timedom); %prezrcalimo karakteristiko zato x2

Y=fft(T);

% Compute the two-sided spectrum P2 and single-sided spectrum P1
P2 = abs(Y/L);       % Normalize the amplitude
P1 = P2(1:L/2+1);    % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double amplitudes (except DC and Nyquist)

f = 1/dt*(0:(L/2))/L;

warning('on','all');
 
% figure;
% plot(f, P1);
% title('Single-Sided Amplitude Spectrum of X(t)');
% xlabel('Frequency (Hz)');
% ylabel('|P1(f)|');
% %xlim([0,1000])

[peaks,loc]=findpeaks(P1,f,'SortStr','descend'); %izluščimo vzb in relek navaor ter ju seštejemo
t=0:pi/180:pi; %izračunomo maksimalni navor in kolesni kot maksimalnega navora

kotMAX=t(find(peaks(1)*sin(t)+peaks(2)*sin(2*t)==max(peaks(1)*sin(t)+peaks(2)*sin(2*t))));
save("kotMAX","kotMAX");

% figure()    
% grid on
% plot(t*180/pi,peaks(1)*sin(t)+peaks(2)*sin(2*t))
% xlabel('mehanski kot [°]')
% ylabel('Navor [Nm]')

f=max(peaks(1)*sin(t)+peaks(2)*sin(2*t));

%f=max(P1);

