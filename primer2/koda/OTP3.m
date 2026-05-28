clc; clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("result.mat","result")

fun= @(vars) OPT(vars(1), vars(2), vars(3), vars(4), vars(5), vars(6));
%options = optimoptions('patternsearch', 'Display', 'iter', 'MeshTolerance', 0.0005, 'StepTolerance', 0.0005,'MaxIterations',50,'PollOrderAlgorithm','Success');
options = optimoptions('patternsearch', 'Display', 'iter', 'MeshTolerance', 0.0005, 'StepTolerance', 0.0005,'MaxIterations',35);
[result, fval] = patternsearch(fun, result, [], [], [], [], [0.6,0.8,1.1,0.8,0.8,-0.05], [1.5,1.1,1.3,0.98,0.98,1], options);

save("result.mat","result")
save("fval.mat","fval")

% load("result") %inicializacija spremenljivk za preizkuse stroja
% result=num2cell(result);
% design12(result{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preizkusi storja

Rezultati();

% M=NavornaKarak()
% RelekNavor()
% SamodrzniNav()
% 
% IndNap()
% KPT()
% 
% load('mere.mat')
% load('konst.mat')
% Lm=(m/pi^2)*mi0*Tp*l*kfe*(kw*N)^2/(2*Pp*d); %izračun magnetilne induktivnosti
% Ldq=DQind()
% 
% Breza()
% Izgube(M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=OPT(AdP,Bd,kE,ni,cosfi,kLok)

load('konst.mat')

Ad=25000*AdP; % A/m tokovna obloga (RMS) 24748-42326

[AdP,Bd,kE,ni,cosfi,kLok]  %izpis trenutnih spremenljik

A=design12(Ad,Bd,kE,ni,cosfi,kLok); %izračun dimanzij stroja
Eror=A(1); 
Rr=A(2);
l=A(3);
V=l*(pi*Rr^2);

if Eror==0  % ce ni erorja nared izris;
    M=izrisstroja(); %izris navorne karakterisitke stroja
else
    M=1;
end

if M<Mn && Eror==0 %če stroj ne proizvede dovolj navora korigiramo dolžino stroja in navitje, da dobimo željen navor

    load('mere.mat')

    l=(1+(Mn-M)/M)*l;
    
    Bda=IN/(d/(mi0)); 
    Ef=Uf*kE;
    N=sqrt(2)*Ef/(2*pi*fk*2/pi*kw*Bda*Tp*l);
    Zq=round(2*a*m*N/Q);

    save('mere.mat')

    M=izrisstroja(); %ponovni izračun storja

    if M>=Mn
        V=l*(pi*Rr^2);
    end

end

K=0;
if M<Mn-1
    K=1;
    Izk=0.8;
else
    Izg=Izgube(M);
    Izk=Izg(1);
    Pd=Izg(2);
end

Izk
M
f=3500*(1-Izk)+8760*V+100000*K/M+10000*Eror; %kriterijska funckija

end
