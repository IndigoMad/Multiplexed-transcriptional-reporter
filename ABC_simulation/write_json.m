load('ABC_simulation.mat')
load('ABC_simulation_degradation.mat')

txt=jsonencode(ABC_simulation);
f=fopen('ABC_simulation.json','w','n','UTF-8');
fprintf(f,txt);
fclose(f);

txt=jsonencode(ABC_simulation_degradation);
f=fopen('ABC_simulation_degradation.json','w','n','UTF-8');
fprintf(f,txt);
fclose(f);