%load mat file and save as json file
load('ABC.mat')
txt=jsonencode(ABC);
f=fopen('ABC.json','w','n','UTF-8');
fprintf(f,txt);
fclose(f);

%update mat file from the json file
f=fopen('ABC.json','r','n','UTF-8');
txt=fgets(f);
fclose(f);
ABC=jsondecode(txt);
save('ABC.mat','ABC')
