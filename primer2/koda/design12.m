function f=design12(Ad,Bd,kE,ni,cosfi,kLok)

% clc;clear all
% load("result.mat")
% AdP=result(1) ;
% Bd=result(2);
% kE=result(3);
% ni=result(4);
% cosfi=result(5);
% kLok=result(6);
% Ad=25000*AdP;

kreza=1;
kvzb=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Konstante
load('konst.mat')

UE=1/kE;
%groba ocena volumna rotorja
Sigma=Ad*Bd*cosfi;
Vr=Mn/(2*Sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MEHANIKA-max dimenzije Rrotorja 
ro=8760; % železo
c=1; % velike izvrtine
omega=2*pi*RPMm/60;
sigma=300000000; % MPa

Drm=2*sqrt(sigma/(c*ro*(omega*kvar)^2));
Rrm=Drm/2;  %maksimalni varni polmer rotorja
    
% Mehanika Dimenzija gredi
sigmaM=550000000;
Ds=2*nthroot(kv*16*Mn/(pi*sigmaM),3);
Rs=Ds/2;  %minimalni varni polmer gredi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zuzanja zanka preveri če je b1 prevelik glede na rotor

ppp=0;
while ppp==0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Oblikovanje rotorja
    Rr=Rs;
    inkrement=0.0001;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % inkrementira cez Rr
    Ads=Ad;
    k1=1;
    count1=1;    
    p1=0;
    while p1==0
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Tokovna obloga --izhodiščna tokovna obloga in realna nista enaki
        %na začetku, zato ju je potrebno skonvergirati
        Ad=Ads;
        delta=100;
        count4=1;
        Tp=pi*(Rr*2)/(2*Pp);
        while delta>1
    
            Ad1=Ad;
    
            l=Vr/(pi*(Rr)^2);
            Is=Pm/(cosfi*ni*Uf*m);
            Ef=Uf*kE;
            N=round(sqrt(2)*Ef/(2*pi*fk*2/pi*kw*Bd*Tp*l));
            %N=sqrt(2)*Ef/(2*pi*fk*kw*Bd*Tp*l);

            Ad=Is*N/(Tp); %/sqrt(2);
    
            count4=count4+1;
            delta=abs(Ad-Ad1);
    
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Zračna reža
    
        kc=1.1;
        d=(1/kc)*Tp*ep*Ad/Bd;
        
        if d<5*10^-4
            d=5*10^-4;
        end
        d=d*kreza;

        % Tp=pi*2*Rr/(2*Pp);
        % delta=1;
        % count3=1;   
        % d1=0.0001; % izhodiščna zračna reža
        % while delta>0.0001 && count3<500 % convergira zračno režo
        % 
        %     kapa=(b1/d1)/(5+(b1/d1));
        %     d=(1/kc)*Tp*ep*Ad/Bd;
        %     delta=abs(d-d1);
        % 
        %     d1=d1+inkrement;
        %     count3=count3+1;
        % end
        % d=d*kreza;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Izhodiščni Rrz -notranji magnento prevodni del rotorja
        bdr=alfa*Bd*pi*Rr/(kc*Bdr*Pp);
        Rrz=bdr/2+Rs;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Vzbujalno navitje
        IN=Bd*d*kvzb/(mi0); 
        Scur=1/Jr;
        Snr=IN*Scur/kpr;  %površina navitja 
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inkrementira cez Rrz
        k2=1;
        count2=1;
        p2=0;
        while p2==0
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Geometrija rotorja
            Rp=Rr-Rrz;
            hpc=kg*Rp; % višina polovega cevlja  %upoštevamo da čevelj drži navitje in mora biti vaj 15% dolzine pola
            br=Rp-hpc;
    
            ar=Rrz*tan(pi/(2*Pp))-bdr/2;
    
            if (ar*br>=Snr && bdr/2+a>=Rrz*tan(pi/(2*Pp)) && (Rr>=Rrz+br+hpc) && ar>0 && br>0 && hpc>0 && ar/br<0.5) || Rrz>=Rr  % zagotavljanje normalnih razmerij navitja
                    p2=1;
            end
            Rrz=Rrz+inkrement;
            count2=count2+1;
        end   
    
        if (ar*br>=Snr && bdr/2+a>=Rrz*tan(pi/(2*Pp)) && ar>0 && br>0 && hpc>0  && ar/br<0.5) || Rr>=Rrm 
            p1=1;
        end
    
        Rr=Rr+inkrement;
        count1=count1+1;
    
    end
    Dr=2*Rr;
    
    if ~(Rrz>=Rr || Rr>=Rrm) %preverimo če je izračunan rotor mogoč

    Eror=0; %če ni damo error

    l1=-1; %konvergiramo l
    delta=0;
    while delta<inkrement
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Korekcija dolžine
        Tu=2*pi*(Dr/2)/Q;
        Tp=pi*Dr/(2*Pp);
        
        Sigma=Ad*Bd*cosfi;
        Vr=Mn/(2*Sigma);
    
        l=Vr/(pi*(Dr/2)^2);
        
        C=pi^2/sqrt(2)*kw*Ad*sqrt(2)*Bd;
        P=ni*cosfi*UE*C*(RPMk/60)*Dr^2*l;
        l=(1+(Pm-P)/P)*l;
        
        P=ni*cosfi*UE*C*(RPMk/60)*Dr^2*l;
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %STATOR
        
        kapa=(b1/d)/(5+(b1/d));
        kc=Tu/(Tu-kapa*b1);
    
        Bdss=Bd;
        bds=(Bd/kc)*Tu/(Bds*kfe);   %debelia zoba zgoraj
        %bdss=(Bd/kc)*Tu/(Bdss*kfe);  %debelia zoba spodaj
        bdss=2*pi*(Rr+d)/Q-b1;
    
        %Ind nap in N
        Ef=Uf*kE;

        Bda=IN/(d/(mi0)); 
        N=sqrt(2)*Ef/(2*pi*fk*kw*2/pi*Bda*Tp*l);
        Zq=round(2*a*m*N/Q);
        %Zq1=2*a*m*N/Q;

        %dimenzije statorskih tokovodnikov
        Is=P/(cosfi*ni*Uf*m);
        Scus=Is/(a*Js);
        Sutor=Zq*Scus/kp;
    
        rs=sqrt(Scus/pi);
        FFs=pi*rs^2*Zq/Sutor;
        
        %dimenzije utorov in jarma
        F=Rr+d+b2;
        bs1=2*pi*(F)/Q-bds; %debelia utora spodaj
        F1=(bs1/2)+(pi*F/Q)-(bds/2); %upoštevamo da je utor trapezen
        Strikot=((bs1-b1)/2)^2*tan(fizob); %upoštevamo dodaten material nad zobom
        hs=(-F1+sqrt(F1^2+4*(pi/Q)*(Sutor+Strikot)))/(2*pi/Q);
        bs2=2*pi*(F+hs)/Q-bds; %debelina utora zgoraj
         
        hsz=(b2+hs);
     
        hys=alfa*Bd*pi*Rr/(2*kc*Bdj*Pp);  %debelina jarma
        
        delta=abs(l-l1);
        l1=l;
    
        D=2*(Rr+d+hsz+hys);
        LD=l/D;
        
        end
        
        if bdss<bds*1.05
            b1=b1-0.0001; %preverimo da je dovolj razdalje med zobi
        else
            ppp=1;
        end
        
        if b1<=0.001
            ppp=1;
            Eror=1;
        end

else
    Eror=1;
    ppp=1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Navitje rotorja
Ink=0.01;
Idc=Ink;
p3=0;
%kvadratni vodniki
while p3==0
    r=sqrt(Idc/Jr);
    Ndc=IN/Idc;
    
    if ((Ndc*(r+2*diz)^2))<=Snr && Idc<=Jr*r^2
        p3=1;
    end

    if r>ar || r>br
        p3=1;
        Eror=1;
    end
    Idc=Idc+Ink;
end
FFr=(Ndc*r^2)/Snr;
Ndc=round(Ndc);

%  %okrogli vodniki
% while p3==0;
%     Nb=(br-2*dizz)/(2*r+2*diz); %dizz izolacija med navitjem in rotorjem
%     Na=(ar-2*dizz)/((r+diz)*sqrt(3));
%     Ndc=floor(Na*Nb/2+(Na-1)*Nb/2);
%     Idc=Jr*r^2*pi;
% 
%     if ((Ndc*pi*(r+diz)^2))<Snr && Idc*Ndc>=IN || r>ar || r>br
%         p3=1;
%     end
%     r=r+inkrement/100;
% end
% FFr=pi*r^2*Ndc/Snr;

if Eror==0
    if tan(fizob)*((bdss-bds)/2)>=bds/2  %preden izrišemo stroj preverimo nekaj parametrov da so možni
        Eror=1;
    end
end

if Eror==0
    if hs<2*rs
        Eror=1;
    end
end

if Eror==0
    angle1=kLok*0.05;  %stvar polovega cevlja
    save("mere.mat","Pp","Q","D","Dr","l","d","Rs","Rrz","hpc","bdr","ar","br","Ndc","r","Jr","hsz","hs","hys","b1","b2","bdss","bds","fizob","Zq","rs","Is","Idc","m","a","Tp","alfa","kLok","angle1","IN","bs1","bs2","N")
end

f=[Eror,Rr,l];
