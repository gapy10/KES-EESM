function f=Izgube(M) %izračun izgub

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% analiza v vogalni toki
load("mere.mat")
load('konst.mat')
load("kotMAX");  %elek kot

k1=1.8;
ke=0.14;
kh=256;

Rpc=Dr/2;
Rr=Dr/2;

IDC=Idc;
IS=Is;
ar=ar;

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

zavrtiROT(kotMAX/Pp*180/pi);

mi_saveas('design1.FEM');
mi_createmesh;
mi_analyze(1);

mi_loadsolution;

R=Rr+d+(hs+b2)/2;
alfa2=2*pi/Q;

y=R*cos(alfa2/2);
x=R*sin(alfa2/2);

alfa1=alfa2;
for i = 1:Q
    A1(i,:) = mo_getpointvalues(x, y);
    y=R*cos(alfa1-alfa2/2);
    x=R*sin(alfa1-alfa2/2);
    alfa1=alfa1+alfa2;
    i=i+1;
end


alfa1=0:2*pi/Q:alfa1-2*pi/Q;
B1=A1(:,3);
B2=A1(:,2);

B=sqrt((B1.^2)+(B2.^2));
clear B1 B2

% figure()
% grid on
% plot(alfa1(1:end-1),B)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Določanje izgub in izkoristka

ROcu=1.724*10^-8; %prevodnost bakra
alfacu=3.93*10^-3; % temperaturni koeficient bakra
alfafe=3.93*10^-3; % temperaturni koeficient pločevine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Izgube v železu

mo_selectblock(0, Rr+d+hsz)
Vs =callfemm('mo_blockintegral(10)');  % blockintegral(1) = Area
mo_selectblock(0, Rr+d+hsz)

n=10000;
R=D/2-hys/2;
alfa1=2*pi;  %/Pp; % celotna krožnica
ink=alfa1/n;

alfa1=0;
for i = 1:n
    x=R*cos(alfa1);
    y=R*sin(alfa1);
    A1(i,:) = mo_getpointvalues(x, y);
    alfa1=alfa1+ink;
    i=i+1;
end

alfa1=(0:ink:alfa1-ink)*180/pi;
B1=A1(:,3);
B2=A1(:,2);

B=sqrt((B1.^2)+(B2.^2));
clear B1 B2

% figure()
% grid on
% hold on
% plot(alfa1,B)

Bsmax=max(B);  % average B^2 over area

Ph=kh*fk*Bsmax^(k1);
Pe=ke*fk^2*Bsmax;
Pfe=(Ph+Pe)*Vs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Izgube v navitju Rot

Rdc=ROcu*Ndc*(2*(bdr+ar)+2*(l+ar))/r^2;
Pdc=(Idc^2*Rdc*(1+alfacu*(60)))*(2*Pp); % 60 °C nad temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Izgube v navitju stat

Lz=2*l*1.2+2*Tp;  %dolžina zanke, 1.2 za glavo navitja
Rac1=ROcu*Lz*N/(pi*rs^2); % DC izgube

glob=sqrt(2/(2*pi*fk*mi0*ROcu));  %izgube skin efekta
Kskin=1+rs/glob;

Rac2=Rac1*Kskin; % dodatne AC izgube

Vcu=Lz*N*pi*rs^2*m*Pp;
Peddy=(pi^2+max(B)^2*(2*rs)^2*fk^2*ROcu/6)*Vcu;

Pac=m*Is^2*Rac2*(1+alfacu*(60))+Peddy;  % 2*bs je dodatna dolžina glave navitja (približek) % 60 °C nad temp + dodatne izgube zaradi izmeničnega bolja v utorih

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Meh Izgube

P0=50;
k=0.5;
Pmeh=P0*(RPMk/6000)^(k+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Izkoristek

Pizg=Pac+Pdc+Pfe+Pmeh;
P=M*2*pi*RPMk/60;
Izk=P/(Pizg+P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gostota moči

mcu=Vcu*9860; %(Vcu+;(2*(l+ar)+2*(bdr+ar))*Ndc*pi*r^2)*9860;

mi_selectcircle(0,0,Rr,4);
mo_selectblock(0,Rr+d+hs);

Vfe=callfemm('mo_blockintegral(10)');
mfe=Vfe*7850;

m=mfe+mcu;
Pd=Pm/m/1000; %kW/kg 

f=[Izk,Pd];

