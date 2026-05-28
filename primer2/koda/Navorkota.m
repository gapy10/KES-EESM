function f=Navorkota(fiel,Pp,kotkarakteristike)
fimeh=0;

i=1;
while i-1<=kotkarakteristike/fiel 
    
    zavrtiROT(fimeh)

    mi_saveas('design1.FEM');
    mi_createmesh;
    mi_analyze(1);
    
    T(i)=-Navor();
    
    fimeh=fiel/Pp;
    i=i+1;
end

f=T;
