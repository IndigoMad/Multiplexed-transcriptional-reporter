%obtaining parameters
load('..\ABC\ABC.mat')
kon1A=ABC.fitting(1).kon1;
konA=ABC.fitting(1).kon;
koffA=ABC.fitting(1).koff;
ktA=ABC.fitting(1).kT;
kpA=10;
kdmA=0.00116;
kdpA=0.000566;
kon1B=ABC.fitting(4).kon1_max;
kaB=ABC.fitting(4).KA;
nB=ABC.fitting(4).n;
konB=ABC.fitting(4).kon;
koffB=ABC.fitting(4).koff;
ktB=ABC.fitting(4).kT;
kpB=10;
kdmB=0.00116;
kdpB=0.000566;
kon1C=ABC.fitting(3).kon1_max;
kaC=ABC.fitting(3).KA;
nC=ABC.fitting(3).n;
konC=ABC.fitting(3).kon;
koffC=ABC.fitting(3).koff;
ktC=ABC.fitting(3).kT;
kpC=10;
kdmC=0.00116;
kdpC=0.000566;
T=3000;
tau=10;
time=(0:10:3000)';

lenc=100;
ABC_simulation_degradation.decription="A transcription simulation of the ABC dataset from the fitted parameters, in addition, a simulation that multiply all kdp and kdm by 100, and divid all ka by 100";
ABC_simulation_degradation.cell_sets(1,1).name='ABC_simulation_degradation';
ABC_simulation_degradation.cell_sets(2,1).name='ABC_random_control';
ABC_simulation_degradation.cell_sets(1).cells(lenc,1).name=[];
ABC_simulation_degradation.cell_sets(2).cells(lenc,1).name=[];
for c=1:lenc
    c
    str=sprintf("simulation_%d",c);
    ABC_simulation_degradation.cell_sets(1).cells(c,1).name=str;
    [A,B,C,Ac,Bc,Cc,Apc,Bpc,Cpc,btA,btB,btC,bsA,bsB,bsC]=TRS(kon1A,konA,koffA,ktA,kpA,kdmA,kdpA,kon1B,kaB,nB,konB,koffB,ktB,kpB,kdmB,kdpB,kon1C,kaC,nC,konC,koffC,ktC,kpC,kdmC,kdpC,T,tau);
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(1,1).name='Gene A';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.intensity=A;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.mRNA=Ac;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.protein=Apc;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(1,1).name='Gene A';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(1,1).time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(1,1).intensity=A;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(1,1).name='Gene A';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(1,1).time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(1,1).intensity=cumsum(A);
    
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(2,1).name='Gene B';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.intensity=B;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.mRNA=Bc;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.protein=Bpc;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(2,1).name='Gene B';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(2,1).time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(2,1).intensity=B;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(2,1).name='Gene B';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(2,1).time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(2,1).intensity=cumsum(B);
    
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(3,1).name='Gene C';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.intensity=C;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.mRNA=Cc;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.protein=Cpc;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(3,1).name='Gene C';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(3,1).time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).sum_intensity(3,1).intensity=C;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(3,1).name='Gene C';
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(3,1).time=time;
    ABC_simulation_degradation.cell_sets(1).cells(c,1).cumulative_intensity(3,1).intensity=cumsum(C);
end

for c=1:lenc
    c
    str=sprintf("degradation_%d",c);
    ABC_simulation_degradation.cell_sets(2).cells(c,1).name=str;
    [A,B,C,Ac,Bc,Cc,Apc,Bpc,Cpc,btA,btB,btC,bsA,bsB,bsC]=TRS(kon1A,konA,koffA,ktA,kpA,kdmA*100,kdpA*100,kon1B,kaB/100,nB,konB,koffB,ktB,kpB,kdmB*100,kdpB*100,kon1C,kaC/100,nC,konC,koffC,ktC,kpC*100,kdmC*100,kdpC,T,tau);
    
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(1,1).name='Gene A';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.intensity=A;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.mRNA=Ac;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.protein=Apc;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(1,1).name='Gene A';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(1,1).time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(1,1).intensity=A;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(1,1).name='Gene A';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(1,1).time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(1,1).intensity=cumsum(A);
    
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(2,1).name='Gene B';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.intensity=B;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.mRNA=Bc;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.protein=Bpc;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(2,1).name='Gene B';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(2,1).time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(2,1).intensity=B;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(2,1).name='Gene B';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(2,1).time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(2,1).intensity=cumsum(B);
    
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(3,1).name='Gene C';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.intensity=C;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.mRNA=Cc;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.protein=Cpc;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(3,1).name='Gene C';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(3,1).time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).sum_intensity(3,1).intensity=C;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(3,1).name='Gene C';
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(3,1).time=time;
    ABC_simulation_degradation.cell_sets(2).cells(c,1).cumulative_intensity(3,1).intensity=cumsum(C);
end

save('ABC_simulation_degradation.mat','ABC_simulation_degradation')