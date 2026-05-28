function f=Breza() %potek gosote magnetnega polja v zraćni reži v neobreenjenem in obremenjenem stanju

load("mere.mat")
load('konst.mat')
load("kotMAX");  %elek kot

Rpc=Dr/2;
Rr=Dr/2;

IDC=Idc;
IS=Is;
ar=ar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%Brez obremenitve

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

mi_saveas('design1.FEM');
mi_createmesh;
mi_analyze(1);

mi_loadsolution;

%oblika B-ja
n=10000;
R=Rr+d/2;
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

figure()
grid on
hold on
plot(alfa1-45,B)
title('Gostota magnetnega polja v zračni reži')
xlabel('Mehanski kot [°]')
ylabel('Gostota magnetnega polja  [T]')
xlim([0,180])

closefemm;

%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Obremenjen

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

%oblika B-ja
n=10000;
R=Rr+d/2;
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

figure()
grid on
hold on
plot(alfa1-45,B)
title('Gostota magnetnega polja v zračni reži')
xlabel('Mehanski kot [°]')
ylabel('Gostota magnetnega polja  [T]')
xlim([0,180])

closefemm;