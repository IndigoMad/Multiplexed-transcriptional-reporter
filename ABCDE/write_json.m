txt=jsonencode(ABCDE);
f=fopen('ABCDE.json','w','n','UTF-8');
fprintf(f,txt);
fclose(f);
