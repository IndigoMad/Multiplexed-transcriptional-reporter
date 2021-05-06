load('ABCDE.mat')

ssnamelist=["Gene A";"Gene B";"Gene C";"Gene D";"Gene E"];
ssupstream=["0";"Gene A";"BE";"0";"Gene D"];
%"0": ustream strating from time 0
%"BE": gene C is regulated by Gene B and C
[lenssnamelist,~]=size(ssnamelist);

[lencs,~]=size(ABCDE.cell_sets);
for cs=1:lencs %cs: cell set
    [lenc,~]=size(ABCDE.cell_sets(cs).cells);
    for c=1:lenc %c: cell
        [lenss,~]=size(ABCDE.cell_sets(cs).cells(c).spot_sets);

        
        %==============================================================   
        %====================Bursts information========================
        %==============================================================
        
        
        %Extracting information from the sum intensity of all spots of a 
        %gene in this cell.
        for ss=1:lenssnamelist
            seq=ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).intensity;
            time=ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).time;
            if isempty(seq)
                seq=ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).time*0;
                ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).intensity=seq;
                ABCDE.cell_sets(cs).cells(c).cumulative_intensity(ss).intensity=seq;
            end
            [bursts,intervals,truncated_seq,truncated_time]=extract_bursts(seq,time);
            ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).bursts=bursts;
            ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).intervals=intervals;
            ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).truncated_seq=truncated_seq;
            ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).truncated_time=truncated_time;
            %starting, the time when the first signal detected.
            if isempty(bursts)
                ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).starting=[];
            else
                ABCDE.cell_sets(cs).cells(c).sum_intensity(ss).starting=bursts(1).starting;
            end
        end
        
        
        %Extracting information from every spots.
        for ss=1:lenss %ss: spot set
            [lens,~]=size(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots);
            for s=1:lens %s: spot
                seq=ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).intensity;
                time=ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).time;
                [bursts,intervals,truncated_seq,truncated_time]=extract_bursts(seq,time);
                ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).bursts=bursts;
                ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).intervals=intervals;
                ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).truncated_seq=truncated_seq;
                ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).truncated_time=truncated_time;
                %starting, the time when the first signal detected.
                if isempty(bursts)
                    ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting=[];
                else
                    ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting=bursts(1).starting;
                end
                %Calculating response time.
                [lensb,~]=size(bursts);
                
                %"0": ustream strating from time 0, Response time: starting of 1st burst of the spot
                if ssupstream(ss)=="0"
                    if lensb>0
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=bursts(1).starting;
                    else
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=[];
                    end
                
                elseif ssupstream(ss)=="BE"
                    [lenBb,~]=size(ABCDE.cell_sets(cs).cells(c).sum_intensity(2).bursts);
                    [lenEb,~]=size(ABCDE.cell_sets(cs).cells(c).sum_intensity(5).bursts);
                    if (lensb>0 & lenBb>0)
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_B=bursts(1).starting-ABCDE.cell_sets(cs).cells(c).sum_intensity(2).bursts(1).starting;
                    else
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_B=[];
                    end
                    if (lensb>0 & lenEb>0)
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_E=bursts(1).starting-ABCDE.cell_sets(cs).cells(c).sum_intensity(5).bursts(1).starting;
                    else
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time_E=[];
                    end
                %Response time: (starting of 1st burst of the spot)-(starting of 1st burst of sum of upstream) 
                else
                    [lenUb,~]=size(ABCDE.cell_sets(cs).cells(c).sum_intensity(ssnamelist==ssupstream(ss)).bursts);
                    if (lensb>0 & lenUb>0)
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=bursts(1).starting-ABCDE.cell_sets(cs).cells(c).sum_intensity(ssnamelist==ssupstream(ss)).bursts(1).starting;
                    else
                        ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=[];
                    end
                end
            end
        end
        %==============================================================   
        %Time statistics: Response time, burst duration, burst interval
        %==============================================================
        %initialization
        ABCDE.cell_sets(cs).cells(c).time_statistics(lenssnamelist,1).name=0;
        ABCDE.cell_sets(cs).cells(c).time_statistics(lenssnamelist).response_time=0;
        ABCDE.cell_sets(cs).cells(c).time_statistics(lenssnamelist).burst_duration=0;
        ABCDE.cell_sets(cs).cells(c).time_statistics(lenssnamelist).burst_interval=0;
        
       
        for ss=1:lenssnamelist
            
            if ss~=3
                ABCDE.cell_sets(cs).cells(c).time_statistics(ss).name=ssnamelist(ss);
                %Response time: excluding response after mitosis, except there is
                %no response before mitosis.

                %A list of response time of this gene in this cell 
                if isempty(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots)
                    response_time_all=[];
                    starting_all=[];
                else
                    response_time_all=[ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).response_time]';
                    starting_all=[ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).starting]';
                end
                

                %the cell did not undergo mitosis
                if isempty(ABCDE.cell_sets(cs).cells(c).mitosis)
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time=response_time_all;

                %the cell underwent mitosis
                else
                    if min(starting_all)>ABCDE.cell_sets(cs).cells(c).mitosis
                        ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time=response_time_all;
                    else
                        %Removing response after mitosis
                        ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time=response_time_all(starting_all<ABCDE.cell_sets(cs).cells(c).mitosis);
                    end
                end
            else
                ABCDE.cell_sets(cs).cells(c).time_statistics(ss).name=ssnamelist(ss);
                %Response time: excluding response after mitosis, except there is
                %no response before mitosis.

                %A list of response time of this gene in this cell 
                if isempty(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots)
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_B=[];
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_E=[];
                else
                    response_time_B_all=[ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).response_time_B]';
                    response_time_E_all=[ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).response_time_E]';
                    starting_all=[ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).starting]';

                    %the cell did not undergo mitosis
                    if isempty(ABCDE.cell_sets(cs).cells(c).mitosis)
                        ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_B=response_time_B_all;
                        ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_E=response_time_E_all;


                    %the cell underwent mitosis
                    else
                        if min(starting_all)>ABCDE.cell_sets(cs).cells(c).mitosis
                            ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_B=response_time_B_all;
                            ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_E=response_time_E_all;
                        else
                            %Removing response after mitosis
                            ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_B=response_time_B_all(starting_all<ABCDE.cell_sets(cs).cells(c).mitosis);
                            ABCDE.cell_sets(cs).cells(c).time_statistics(ss).response_time_E=response_time_E_all(starting_all<ABCDE.cell_sets(cs).cells(c).mitosis);
                        end
                    end
                end
            end
            
            %Burst_duration
            if isempty(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots)
                ABCDE.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=[];
            else
                burst_duration=vertcat(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).bursts);
                if isempty(burst_duration)
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=[];
                else
                    burst_duration=vertcat(burst_duration(:).duration);
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=burst_duration;
                end
            end

            %Burst_interval
            if isempty(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots)
                ABCDE.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=[];
            else
                burst_interval=vertcat(ABCDE.cell_sets(cs).cells(c).spot_sets(ss).spots(:).intervals);
                if isempty(burst_interval)
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).burst_interval=[];
                else
                    burst_interval=vertcat(burst_interval(:).duration);
                    ABCDE.cell_sets(cs).cells(c).time_statistics(ss).burst_interval=burst_interval;
                end
            end
        end
    end
    
    %==============================================================   
    %===============Time statistics of cell sets===================
    %==============================================================
    
    %initialization
    ABCDE.cell_sets(cs).time_statistics(lenssnamelist,1).name=0;
    ABCDE.cell_sets(cs).time_statistics(lenssnamelist,1).response_time=0;
    ABCDE.cell_sets(cs).time_statistics(lenssnamelist,1).burst_duration=0;
    ABCDE.cell_sets(cs).time_statistics(lenssnamelist,1).burst_interval=0;
    
    spot_set_time_statistics=vertcat(ABCDE.cell_sets(cs).cells(:).time_statistics);
    for ss=1:lenssnamelist
        ABCDE.cell_sets(cs).time_statistics(ss,1).name=ssnamelist(ss);
        if ss~=3
            ABCDE.cell_sets(cs).time_statistics(ss,1).response_time=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).response_time);
        else
            ABCDE.cell_sets(cs).time_statistics(ss,1).response_time_B=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).response_time_B);
            ABCDE.cell_sets(cs).time_statistics(ss,1).response_time_E=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).response_time_E);
        end
        ABCDE.cell_sets(cs).time_statistics(ss,1).burst_duration=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).burst_duration);
        ABCDE.cell_sets(cs).time_statistics(ss,1).burst_interval=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).burst_interval);
    end
    
    
end
                
save('ABCDE.mat','ABCDE')






%given a intensity sequence, extracting bursts information.
function[bursts,intervals,truncated_seq,truncated_time]=extract_bursts(seq,time)
    bi_seq=(seq>0);%binarized sequence
    diff_bi_seq=[bi_seq(1);diff(bi_seq)];%differenced binarized sequence.
    
%     =============Here is a discarded algorithm===============
%     % 1 represents the starting of a burst while -1 represents ending of a burst.
%     dd_seq=diff(diff_bi_seq);% difference the differenced sequence again.
%     % 2 represents a 1-frame gap between two burst, which is belived to be a
%     % missed signal but not a real ending and starting of bursts.
%     gap_filter=([dd_seq;0]~=2).*([0;dd_seq]~=2); % gaps filter
%     filtered_burst=diff_bi_seq.*gap_filter;

%   =============Current algorithm===============
    filtered_burst=diff_bi_seq;
    burst_starting=time(filtered_burst>0);%1 represent the starting of a burst
    burst_ending=time(filtered_burst<0);%-1 represent the ending of a burst
    
    
    [bn,~]=size(burst_ending);
    burst_starting=burst_starting(1:bn);%remove unfinished bursts at the end of the time

    if bn>0
        %initialization
        bursts(bn,1).starting=0;
        bursts(bn).ending=0;
        bursts(bn).duration=0;
        
        for i=1:bn
            bursts(i).starting=burst_starting(i);
            bursts(i).ending=burst_ending(i);
            bursts(i).duration=burst_ending(i)-burst_starting(i);
        end
    else
        bursts=[];
    end
    if bn-1>0
        %initialization
        intervals(bn-1,1).starting=0;
        intervals(bn-1).ending=0;
        intervals(bn-1).duration=0;
        
        for i=1:bn-1
            intervals(i).starting=burst_ending(i);
            intervals(i).ending=burst_starting(i+1);
            intervals(i).duration=burst_starting(i+1)-burst_ending(i);
        end
    else
        intervals=[];
    end
    
    if bn>0
        time_end=find(time==burst_ending(bn),1)-1;
        truncated_time=time(1:time_end);
        truncated_seq=seq(1:time_end);
    else
        truncated_time=[];
        truncated_seq=[];
    end
end



