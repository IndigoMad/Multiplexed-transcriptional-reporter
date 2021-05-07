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

lenc=1000;
ABC_simulation.decription="A transcription simulation of the ABC dataset from the fitted parameters, in addition, a simulation of the same parameter but without regulation relationship.";
ABC_simulation.cell_sets(1,1).name='ABC_simulation';
ABC_simulation.cell_sets(2,1).name='ABC_random_control';
ABC_simulation.cell_sets(1).cells(lenc,1).name=[];
ABC_simulation.cell_sets(2).cells(lenc,1).name=[];
for c=1:lenc
    c
    str=sprintf("simulation_%d",c);
    ABC_simulation.cell_sets(1).cells(c,1).name=str;
    [A,B,C,Ac,Bc,Cc,Apc,Bpc,Cpc,btA,btB,btC,bsA,bsB,bsC]=TRS(kon1A,konA,koffA,ktA,kpA,kdmA,kdpA,kon1B,kaB,nB,konB,koffB,ktB,kpB,kdmB,kdpB,kon1C,kaC,nC,konC,koffC,ktC,kpC,kdmC,kdpC,T,tau);
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(1,1).name='Gene A';
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.time=time;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.intensity=A;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.mRNA=Ac;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(1,1).spots.protein=Apc;
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(1,1).name='Gene A';
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(1,1).time=time;
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(1,1).intensity=A;
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(1,1).name='Gene A';
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(1,1).time=time;
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(1,1).intensity=cumsum(A);
    
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(2,1).name='Gene B';
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.time=time;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.intensity=B;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.mRNA=Bc;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(2,1).spots.protein=Bpc;
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(2,1).name='Gene B';
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(2,1).time=time;
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(2,1).intensity=B;
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(2,1).name='Gene B';
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(2,1).time=time;
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(2,1).intensity=cumsum(B);
    
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(3,1).name='Gene C';
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.time=time;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.intensity=C;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.mRNA=Cc;
    ABC_simulation.cell_sets(1).cells(c,1).spot_sets(3,1).spots.protein=Cpc;
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(3,1).name='Gene C';
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(3,1).time=time;
    ABC_simulation.cell_sets(1).cells(c,1).sum_intensity(3,1).intensity=C;
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(3,1).name='Gene C';
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(3,1).time=time;
    ABC_simulation.cell_sets(1).cells(c,1).cumulative_intensity(3,1).intensity=cumsum(C);
end

for c=1:lenc
    c
    str=sprintf("random_control_%d",c);
    ABC_simulation.cell_sets(2).cells(c,1).name=str;
    [A,B,C,Ac,Bc,Cc,Apc,Bpc,Cpc,btA,btB,btC,bsA,bsB,bsC]=TRSrandomcontrol_full(kon1A,konA,koffA,ktA,kpA,kdmA,kdpA,kon1B,konB,koffB,ktB,kpB,kdmB,kdpB,kon1C,konC,koffC,ktC,kpC,kdmC,kdpC,T,tau);
    
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(1,1).name='Gene A';
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.time=time;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.intensity=A;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.mRNA=Ac;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(1,1).spots.protein=Apc;
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(1,1).name='Gene A';
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(1,1).time=time;
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(1,1).intensity=A;
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(1,1).name='Gene A';
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(1,1).time=time;
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(1,1).intensity=cumsum(A);
    
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(2,1).name='Gene B';
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.time=time;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.intensity=B;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.mRNA=Bc;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(2,1).spots.protein=Bpc;
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(2,1).name='Gene B';
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(2,1).time=time;
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(2,1).intensity=B;
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(2,1).name='Gene B';
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(2,1).time=time;
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(2,1).intensity=cumsum(B);
    
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(3,1).name='Gene C';
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.time=time;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.intensity=C;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.mRNA=Cc;
    ABC_simulation.cell_sets(2).cells(c,1).spot_sets(3,1).spots.protein=Cpc;
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(3,1).name='Gene C';
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(3,1).time=time;
    ABC_simulation.cell_sets(2).cells(c,1).sum_intensity(3,1).intensity=C;
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(3,1).name='Gene C';
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(3,1).time=time;
    ABC_simulation.cell_sets(2).cells(c,1).cumulative_intensity(3,1).intensity=cumsum(C);
end

save('ABC_simulation.mat','ABC_simulation')
