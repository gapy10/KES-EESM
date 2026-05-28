function f=DQind()
%izračun prečne in vzdolžne induktivnosti v lineanrem in obremenjen stanju
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

Idc=0;

Is1=sqrt(2)*Is;
Is2=0;
Is3=0;

ar=ar;

%Izris stroja
izris2(Pp,Q,D,Dr,d,Rs,Rrz,hpc,bdr,ar,br,Ndc,hsz,hs,hys,b1,b2,bdss,bds,fizob,Zq,Idc,m,Is1,Is2,Is3,alfa,Rpc,angle1)
%Zagon izracuna
mi_selectcircle(0,0,Rr,4); %group 1 je cel rotor
mi_setgroup(1);
mi_clearselected();

theta=180/Pp; % gledan kot
theta1=theta; % št izračunov na električno periodo
ink=theta/(theta1); %20-vsake 4.5 stopinje (meh)
fi=0:ink:theta; % ° ugibam da meh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nasičen del
i=1;
while i<=length(fi)

mi_saveas('design1.FEM');
mi_createmesh;
mi_analyze(1);

mi_loadsolution;

flux1 = callfemm('mo_getcircuitproperties("NavitjeAC1")');
flux(i)=flux1(3);

zavrtiROT(ink);
i=i+1;
end

Ldn=max(flux)/Is;
Lqn=min(flux)/Is;

figure()
grid on
hold on
plot(fi,flux)
title('Magnetna sklopitev')
xlabel('Električni kot [°]')
ylabel('Magnetni sklop  [Wb]')

%f=[Ldn,Lqn];

closefemm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%linearni del

openfemm;
newdocument(0); % 0 for magnetics problem
smartmesh(0);
rezolucija=1e-8; % original 1e-8
minangle=30;
mi_probdef(0, 'meters', 'planar', rezolucija, l, minangle);

Rpc=Dr/2;
Rr=Dr/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Navorna karakteristike

Idc=0;

Is=Is/10;

Is1=sqrt(2)*Is;
Is2=0;
Is3=0;

%Izris stroja
izris2(Pp,Q,D,Dr,d,Rs,Rrz,hpc,bdr,ar,br,Ndc,hsz,hs,hys,b1,b2,bdss,bds,fizob,Zq,Idc,m,Is1,Is2,Is3,alfa,Rpc,angle1)
%Zagon izracuna
mi_selectcircle(0,0,Rr,4); %group 1 je cel rotor
mi_setgroup(1);
mi_clearselected();

i=1;
while i<=length(fi)

mi_saveas('design1.FEM');
mi_createmesh;
mi_analyze(1);

mi_loadsolution;

flux1 = callfemm('mo_getcircuitproperties("NavitjeAC1")');
flux(i)=flux1(3);

zavrtiROT(ink);
i=i+1;
end

Ldl=max(flux)/Is;
Lql=min(flux)/Is;

figure()
grid on
hold on
plot(fi,flux)
title('Magnetna sklopitev')
xlabel('Električni kot [°]')
ylabel('Magnetni sklop  [Wb]')


f=[Ldn,Lqn,Ldl,Lql];

closefemm;