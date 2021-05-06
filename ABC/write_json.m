txt=jsonencode(ABC);
f=fopen('ABC.json','w','n','UTF-8');
fprintf(f,txt);
fclose(f);
