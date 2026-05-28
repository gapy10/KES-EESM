function f=Navor()
%izracun navora
mi_loadsolution;

mo_groupselectblock(1);
Torq = mo_blockintegral(22); 
%Torq=mo_gapintegral('Torque',1);
mo_clearblock();

f=Torq;