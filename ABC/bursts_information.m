ssnamelist=["Gene A";"Gene B";"Gene C"];
ssupstream=["0";"Gene A";"Gene B";"NA"];
%"0": ustream strating from time 0
%"NA": could not defind upstream
[lenssnamelist,~]=size(ssnamelist);

[lencs,~]=size(ABC.cell_sets);
for cs=1:lencs %cs: cell set
    [lenc,~]=size(ABC.cell_sets(cs).cells);
    for c=1:lenc %c: cell
        [lenss,~]=size(ABC.cell_sets(cs).cells(c).spot_sets);

        
        %==============================================================   
        %====================Bursts information========================
        %==============================================================
        
        
        %Extracting information from the sum intensity of all spots of a 
        %gene in this cell.
        for ss=1:lenssnamelist
            seq=ABC.cell_sets(cs).cells(c).sum_intensity(ss).intensity;
            time=ABC.cell_sets(cs).cells(c).sum_intensity(ss).time;
            [bursts,intervals,truncated_seq,truncated_time]=extract_bursts(seq,time);
            ABC.cell_sets(cs).cells(c).sum_intensity(ss).bursts=bursts;
            ABC.cell_sets(cs).cells(c).sum_intensity(ss).intervals=intervals;
            ABC.cell_sets(cs).cells(c).sum_intensity(ss).truncated_seq=truncated_seq;
            ABC.cell_sets(cs).cells(c).sum_intensity(ss).truncated_time=truncated_time;
            %starting, the time when the first signal detected.
            if isempty(bursts)
                ABC.cell_sets(cs).cells(c).sum_intensity(ss).starting=[];
            else
                ABC.cell_sets(cs).cells(c).sum_intensity(ss).starting=bursts(1).starting;
            end
        end
        
        
        %Extracting information from every spots.
        for ss=1:lenss %ss: spot set
            [lens,~]=size(ABC.cell_sets(cs).cells(c).spot_sets(ss).spots);
            for s=1:lens %s: spot
                seq=ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).intensity;
                time=ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).time;
                [bursts,intervals,truncated_seq,truncated_time]=extract_bursts(seq,time);
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).bursts=bursts;
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).intervals=intervals;
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).truncated_seq=truncated_seq;
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).truncated_time=truncated_time;
                %starting, the time when the first signal detected.
                if isempty(bursts)
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting=[];
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).burst_duration=[];
                else
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).starting=bursts(1).starting;
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).burst_duration=vertcat(bursts(:).duration);
                end
                if isempty(intervals)
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).burst_interval=[];
                else
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).burst_interval=vertcat(intervals(:).duration);
                end
                
                
                %Calculating response time.
                [lensb,~]=size(bursts);
                
                %"0": ustream strating from time 0, Response time: starting of 1st burst of the spot
                if ssupstream(ss)=="0"
                    if lensb>0
                        ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=bursts(1).starting;
                    else
                        ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=[];
                    end
                
                %"NA": could not defind upstream
                elseif ssupstream(ss)=="NA"
                    ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=[];
                    
                %Response time: (starting of 1st burst of the spot)-(starting of 1st burst of sum of upstream) 
                else
                    [lenUb,~]=size(ABC.cell_sets(cs).cells(c).sum_intensity(ssnamelist==ssupstream(ss)).bursts);
                    if (lensb>0 & lenUb>0)
                        ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=bursts(1).starting-ABC.cell_sets(cs).cells(c).sum_intensity(ssnamelist==ssupstream(ss)).bursts(1).starting;
                    else
                        ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).response_time=[];
                    end
                end
            end
        end
        %==============================================================   
        %Time statistics: Response time, burst duration, burst interval
        %==============================================================
        %initialization
        ABC.cell_sets(cs).cells(c).time_statistics(lenssnamelist,1).name=0;
        ABC.cell_sets(cs).cells(c).time_statistics(lenssnamelist).response_time=0;
        ABC.cell_sets(cs).cells(c).time_statistics(lenssnamelist).burst_duration=0;
        ABC.cell_sets(cs).cells(c).time_statistics(lenssnamelist).burst_interval=0;
        
       
        for ss=1:lenssnamelist
            
            ABC.cell_sets(cs).cells(c).time_statistics(ss).name=ssnamelist(ss);
            %Response time: excluding response after mitosis, except there is
            %no response before mitosis.
            
            %A list of response time of this gene in this cell 
            
            response_time_all=[ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(:).response_time]';
            starting_all=[ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(:).starting]';
            
            %the cell did not undergo mitosis
            if isempty(ABC.cell_sets(cs).cells(c).mitosis)
                ABC.cell_sets(cs).cells(c).time_statistics(ss).response_time=response_time_all;
            
            %the cell underwent mitosis
            else
                if min(starting_all)>ABC.cell_sets(cs).cells(c).mitosis
                    ABC.cell_sets(cs).cells(c).time_statistics(ss).response_time=response_time_all;
                else
                    %Removing response after mitosis
                    ABC.cell_sets(cs).cells(c).time_statistics(ss).response_time=response_time_all(starting_all<ABC.cell_sets(cs).cells(c).mitosis);
                end
            end
            
            
             %Burst_duration
             burst_duration=vertcat(ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(:).bursts);
             if isempty(burst_duration)
                 ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=[];
             else
                 burst_duration=vertcat(burst_duration(:).duration);
                 ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=burst_duration;
             end
             
             %Burst_interval
             burst_interval=vertcat(ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(:).intervals);
             if isempty(burst_interval)
                 ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_interval=[];
             else
                 burst_interval=vertcat(burst_interval(:).duration);
                 ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_interval=burst_interval;
             end
        end
    end
    
    %==============================================================   
    %===============Time statistics of cell sets===================
    %==============================================================
    
    %initialization
    ABC.cell_sets(cs).time_statistics(lenssnamelist,1).name=0;
    ABC.cell_sets(cs).time_statistics(lenssnamelist,1).response_time=0;
    ABC.cell_sets(cs).time_statistics(lenssnamelist,1).burst_duration=0;
    ABC.cell_sets(cs).time_statistics(lenssnamelist,1).burst_interval=0;
    
    spot_set_time_statistics=vertcat(ABC.cell_sets(cs).cells(:).time_statistics);
    for ss=1:lenssnamelist
        ABC.cell_sets(cs).time_statistics(ss,1).name=ssnamelist(ss);
        ABC.cell_sets(cs).time_statistics(ss,1).response_time=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).response_time);
        ABC.cell_sets(cs).time_statistics(ss,1).burst_duration=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).burst_duration);
        ABC.cell_sets(cs).time_statistics(ss,1).burst_interval=vertcat(spot_set_time_statistics(vertcat(spot_set_time_statistics.name)==ssnamelist(ss)).burst_interval);
    end
    
    
end
                
save('ABC.mat','ABC')






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




