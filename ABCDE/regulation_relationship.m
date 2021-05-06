load('ABCDE.mat')

ssnamelist=["Gene A";"Gene B";"Gene C";"Gene D";"Gene E"];
ssupstream=["0";"Gene A";"BE";"0";"Gene D"];
%"0": ustream strating from time 0
%"BE": gene C is regulated by Gene B and C
[lenssnamelist,~]=size(ssnamelist);

[lencs,~]=size(ABCDE.cell_sets);
for cs=1:lencs %cs: cell set
    temcs=ABCDE.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    
    %initialization
    ABCDE.cell_sets(cs).regulation_relationship(2,1).name=0;
    ABCDE.cell_sets(cs).regulation_relationship(2).response_time=0;
    ABCDE.cell_sets(cs).regulation_relationship(2).response_time_upstream_rate=0;
    ABCDE.cell_sets(cs).regulation_relationship(2).response_time_upstream_cumulation=0;
    
    for c=1:lenc %c: cell
        temc=temcs.cells(c);
        [lenss,~]=size(temc.spot_sets);
        %==============================================================   
        %==========Response time upstream rate and cumulation===========
        %==============================================================
        %initialization
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3,1).name=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_upstream_rate=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_upstream_cumulation=0;
        ABCDE.cell_sets(cs).regulation_relationship(3).starting=0;
        ABCDE.cell_sets(cs).regulation_relationship(3).starting_upstream=0;
        
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_B=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_upstream_rate_B=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_upstream_cumulation_B=0;
        ABCDE.cell_sets(cs).regulation_relationship(3).starting_B=0;
        ABCDE.cell_sets(cs).regulation_relationship(3).starting_upstream_B=0;
        
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_E=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_upstream_rate_E=0;
        ABCDE.cell_sets(cs).cells(c).regulation_relationship(3).response_time_upstream_cumulation_E=0;
        ABCDE.cell_sets(cs).regulation_relationship(3).starting_E=0;
        ABCDE.cell_sets(cs).regulation_relationship(3).starting_upstream_E=0;
        %Upstream rate of response time, definded as the cumulative
        %intensity before response divided by reponse time, in another
        %word, the average intensity from the first upstream signal to 
        %downstream response. 
        
        %Response time: excluding response after mitosis, except there is
        %no response before mitosis.
        
        %Considering only gene B C and E
        sscounter=0;
        for ss=1:lenss
            if ssupstream(ss)=="0"
                1;
            elseif ssupstream(ss)=="BE"
                sscounter=sscounter+1;
                ci=temc.cumulative_intensity;
                ci_B=ci(2).intensity;
                ci_t_B=ci(2).time;
                ust_B=temc.sum_intensity(2).starting;
                
                ci_E=ci(5).intensity;
                ci_t_E=ci(5).time;
                ust_E=temc.sum_intensity(5).starting;
                
                temss=temc.spot_sets(ss);
                [lens,~]=size(temss.spots);
                
                temst_B=ones(lens,1)*999999999;
                temrt_B=ones(lens,1)*999999999;
                temrtur_B=ones(lens,1)*999999999;
                temrtuc_B=ones(lens,1)*999999999;
                temust_B=ones(lens,1)*999999999;

                temst_E=ones(lens,1)*999999999;
                temrt_E=ones(lens,1)*999999999;
                temrtur_E=ones(lens,1)*999999999;
                temrtuc_E=ones(lens,1)*999999999;
                temust_E=ones(lens,1)*999999999;
                
                for s=1:lens
                    tems=temss.spots(s);
                    rt_B=tems.response_time_B;%response time
                    rt_E=tems.response_time_E;%response time
                    if isempty(rt_B)
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate_B=[];
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation_B=[];
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream_B=[];
                    else
                        st=tems.starting;%starting time
                        rtur_B=ci_B(ci_t_B==st)/rt_B;%response time upstream rate
                        
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate_B=rtur_B;
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation_B=ci_B(ci_t_B==st);
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream_B=ust_B;
                        
                        temst_B(s)=st;
                        temrt_B(s)=rt_B;
                        temrtur_B(s)=rtur_B;
                        temrtuc_B(s)=ci_B(ci_t_B==st);
                        temust_B(s)=ust_B;
                    end
                    
                    if isempty(rt_E)
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate_E=[];
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation_E=[];
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream_E=[];
                    else
                        st=tems.starting;%starting time
                        rtur_E=ci_E(ci_t_E==st)/rt_E;%response time upstream rate
                        
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate_E=rtur_E;
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation_E=ci_E(ci_t_E==st);
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream_E=ust_E;
                        
                        temst_E(s)=st;
                        temrt_E(s)=rt_E;
                        temrtur_E(s)=rtur_E;
                        temrtuc_E(s)=ci_E(ci_t_E==st);
                        temust_E(s)=ust_E;
                    end
                end
                %removal of empty value
                temst_B=temst_B(temst_B~=999999999);
                temrt_B=temrt_B(temrt_B~=999999999);
                temrtur_B=temrtur_B(temrtur_B~=999999999);
                temrtuc_B=temrtuc_B(temrtuc_B~=999999999);
                temust_B=temust_B(temust_B~=999999999);

                temst_E=temst_E(temst_E~=999999999);
                temrt_E=temrt_E(temrt_E~=999999999);
                temrtur_E=temrtur_E(temrtur_E~=999999999);
                temrtuc_E=temrtuc_E(temrtuc_E~=999999999);
                temust_E=temust_E(temust_E~=999999999);
                %Response time: excluding response after mitosis, except there is
                %no response before mitosis.

                %the cell did not undergo mitosis
                if isempty(ABCDE.cell_sets(cs).cells(c).mitosis)
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).name=ssnamelist(ss);
                    
                    
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_B=temrt_B;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate_B=temrtur_B;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation_B=temrtuc_B;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_B=temst_B;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream_B=temust_B;
                    
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_E=temrt_E;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate_E=temrtur_E;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation_E=temrtuc_E;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_E=temst_E;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream_E=temust_E;

                %the cell underwent mitosis
                else
                    if min(temst)>ABCDE.cell_sets(cs).cells(c).mitosis
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).name=ssnamelist(ss);
                        
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_B=temrt_B;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate_B=temrtur_B;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation_B=temrtuc_B;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_B=temst_B;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream_B=temust_B;
                        
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_E=temrt_E;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate_E=temrtur_E;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation_E=temrtuc_E;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_E=temst_E;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream_E=temust_E;
                        
                    else
                        %Removing response after mitosis
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).name=ssnamelist(ss);
                        
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_B=temrt_B(temst_B<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate_B=temrtur_B(temst_B<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation_B=temrtuc_B(temst_B<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_B=temst_B(temst_B<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream_B=temust_B(temst_B<ABCDE.cell_sets(cs).cells(c).mitosis);
                        
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_E=temrt_E(temst_E<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate_E=temrtur_E(temst_E<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation_E=temrtuc_E(temst_E<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_E=temst_E(temst_E<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream_E=temust_E(temst_E<ABCDE.cell_sets(cs).cells(c).mitosis);
                    end
                end
            else
                sscounter=sscounter+1;
                ci=temc.cumulative_intensity;
                ci_up=ci(ssnamelist==ssupstream(ss)).intensity;
                ci_t_up=ci(ssnamelist==ssupstream(ss)).time;
                ust=temc.sum_intensity(ssnamelist==ssupstream(ss)).starting;

                temss=temc.spot_sets(ss);
                [lens,~]=size(temss.spots);
                temst=ones(lens,1)*999999999;
                temrt=ones(lens,1)*999999999;
                temrtur=ones(lens,1)*999999999;
                temrtuc=ones(lens,1)*999999999;
                temust=ones(lens,1)*999999999;
                for s=1:lens
                    tems=temss.spots(s);
                    rt=tems.response_time;%response time

                    if isempty(rt)
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate=[];
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation=[];
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream=[];
                    else
                        st=tems.starting;%starting time
                        rtur=ci_up(ci_t_up==st)/rt;%response time upstream rate
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate=rtur;
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation=ci_up(ci_t_up==st);
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream=ust;
                        temst(s)=st;
                        temrt(s)=rt;
                        temrtur(s)=rtur;
                        temrtuc(s)=ci_up(ci_t_up==st);
                        temust(s)=ust;
                    end
                end
                %removal of empty value
                temst=temst(temst~=999999999);
                temrt=temrt(temrt~=999999999);
                temrtur=temrtur(temrtur~=999999999);
                temrtuc=temrtuc(temrtuc~=999999999);
                temust=temust(temust~=999999999);

                %Response time: excluding response after mitosis, except there is
                %no response before mitosis.

                %the cell did not undergo mitosis
                if isempty(ABCDE.cell_sets(cs).cells(c).mitosis)
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).name=ssnamelist(ss);
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time=temrt;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate=temrtur;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation=temrtuc;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting=temst;
                    ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream=temust;

                %the cell underwent mitosis
                else
                    if min(temst)>ABCDE.cell_sets(cs).cells(c).mitosis
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).name=ssnamelist(ss);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time=temrt;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate=temrtur;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation=temrtuc;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting=temst;
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream=temust;
                    else
                        %Removing response after mitosis
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).name=ssnamelist(ss);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time=temrt(temst<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_rate=temrtur(temst<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).response_time_upstream_cumulation=temrtuc(temst<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting=temst(temst<ABCDE.cell_sets(cs).cells(c).mitosis);
                        ABCDE.cell_sets(cs).cells(c).regulation_relationship(sscounter).starting_upstream=temust(temst<ABCDE.cell_sets(cs).cells(c).mitosis);
                    end
                end
            end
        end
    end
    
    
    cell_sets_regulation_relationship=vertcat(ABCDE.cell_sets(cs).cells(:).regulation_relationship);
    sscounter=0;
    for ss=1:lenss
        if ssupstream(ss)=="0"
            1;
        elseif ssupstream(ss)=="BE"
            sscounter=sscounter+1;
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).name=ssnamelist(ss);
            
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_B=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_B);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_upstream_rate_B=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_rate_B);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_upstream_cumulation_B=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_cumulation_B);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).starting_B=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting_B);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).starting_upstream_B=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting_upstream_B);
            
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_E=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_E);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_upstream_rate_E=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_rate_E);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_upstream_cumulation_E=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_cumulation_E);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).starting_E=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting_E);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).starting_upstream_E=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting_upstream_E);
        else
            sscounter=sscounter+1;
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).name=ssnamelist(ss);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_upstream_rate=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_rate);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).response_time_upstream_cumulation=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_cumulation);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).starting=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting);
            ABCDE.cell_sets(cs).regulation_relationship(sscounter).starting_upstream=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting_upstream);
        end
    end
end
                
save('ABCDE.mat','ABCDE')










