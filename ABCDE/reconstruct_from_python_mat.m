%load mat file and save as json file
load('ABCDE.mat')
txt=jsonencode(ABCDE);
f=fopen('ABCDE.json','w','n','UTF-8');
fprintf(f,txt);
fclose(f);

%update mat file from the json file
f=fopen('ABCDE.json','r','n','UTF-8');
txt=fgets(f);
fclose(f);
ABCDE=jsondecode(txt);
save('ABCDE.mat','ABCDE')