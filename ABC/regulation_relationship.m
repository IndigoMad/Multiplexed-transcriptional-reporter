load('ABC.mat')

ssnamelist=["Gene A";"Gene B";"Gene C"];
ssupstream=["0";"Gene A";"Gene B";"NA"];
%"0": upstream strating from time 0
%"NA": could not defind upstream
[lenssnamelist,~]=size(ssnamelist);

[lencs,~]=size(ABC.cell_sets);
for cs=1:lencs %cs: cell set
    temcs=ABC.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    
    %initialization
    ABC.cell_sets(cs).regulation_relationship(2,1).name=0;
    ABC.cell_sets(cs).regulation_relationship(2).response_time=0;
    ABC.cell_sets(cs).regulation_relationship(2).response_time_upstream_rate=0;
    ABC.cell_sets(cs).regulation_relationship(2).response_time_upstream_cumulation=0;
    ABC.cell_sets(cs).regulation_relationship(2).starting=0;
    ABC.cell_sets(cs).regulation_relationship(2).starting_upstream=0;
    ABC.cell_sets(cs).regulation_relationship(2).burst_propogation=0;
    
    for c=1:lenc %c: cell
        temc=temcs.cells(c);
        [lenss,~]=size(temc.spot_sets);
        %==============================================================   
        %==========Response time upstream rate and cumulation===========
        %==============================================================
        %initialization
        ABC.cell_sets(cs).cells(c).regulation_relationship(2,1).name=0;
        ABC.cell_sets(cs).cells(c).regulation_relationship(2).response_time=0;
        ABC.cell_sets(cs).cells(c).regulation_relationship(2).response_time_upstream_rate=0;
        ABC.cell_sets(cs).cells(c).regulation_relationship(2).response_time_upstream_cumulation=0;
        ABC.cell_sets(cs).cells(c).regulation_relationship(2).starting=0;
        ABC.cell_sets(cs).cells(c).regulation_relationship(2).starting_upstream=0;
        ABC.cell_sets(cs).cells(c).regulation_relationship(2).burst_propogation=0;
        

        
        %Upstream rate of response time, definded as the cumulative
        %intensity before response divided by reponse time, in another
        %word, the average intensity from the first upstream signal to 
        %downstream response. 
        
        %Response time: excluding response after mitosis, except there is
        %no response before mitosis.
        
        %Considering only gene B and C
        for ss=2:3
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
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate=[];
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation=[];
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream=[];
                else
                    st=tems.starting;%starting time
                    rtur=ci_up(ci_t_up==st)/rt;%response time upstream rate
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_rate=rtur;
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_upstream_cumulation=ci_up(ci_t_up==st);
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting_upstream=ust;
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
            if isempty(ABC.cell_sets(cs).cells(c).mitosis)
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).name=ssnamelist(ss);
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time=temrt;
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_rate=temrtur;
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_cumulation=temrtuc;
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting=temst;
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting_upstream=temust;

            %the cell underwent mitosis
            else
                if min(temst)>ABC.cell_sets(cs).cells(c).mitosis
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).name=ssnamelist(ss);
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time=temrt;
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_rate=temrtur;
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_cumulation=temrtuc;
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting=temst;
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting_upstream=temust;
                else
                    %Removing response after mitosis
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).name=ssnamelist(ss);
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time=temrt(temst<ABC.cell_sets(cs).cells(c).mitosis);
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_rate=temrtur(temst<ABC.cell_sets(cs).cells(c).mitosis);
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_cumulation=temrtuc(temst<ABC.cell_sets(cs).cells(c).mitosis);
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting=temst(temst<ABC.cell_sets(cs).cells(c).mitosis);
                    ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting_upstream=temust(temst<ABC.cell_sets(cs).cells(c).mitosis);
                end
            end
        end
    end
    
    
    cell_sets_regulation_relationship=vertcat(ABC.cell_sets(cs).cells(:).regulation_relationship);
    for ss=2:3
        ABC.cell_sets(cs).regulation_relationship(ss-1).name=ssnamelist(ss);
        ABC.cell_sets(cs).regulation_relationship(ss-1).response_time=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time);
        ABC.cell_sets(cs).regulation_relationship(ss-1).response_time_upstream_rate=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_rate);
        ABC.cell_sets(cs).regulation_relationship(ss-1).response_time_upstream_cumulation=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).response_time_upstream_cumulation);
        ABC.cell_sets(cs).regulation_relationship(ss-1).starting=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting);
        ABC.cell_sets(cs).regulation_relationship(ss-1).starting_upstream=vertcat(cell_sets_regulation_relationship(vertcat(cell_sets_regulation_relationship.name)==ssnamelist(ss)).starting_upstream);
    end
end


save('ABC.mat','ABC')










