load('ABC.mat')

ssnamelist=["Gene A";"Gene B";"Gene C";"Gene B fast";"Gene B slow"];
ssupstream=["0";"Gene A";"Gene B";"Gene A";"Gene A"];
ssupstreami=[0;1;2;1;1];
%"0": upstream strating from time 0
%"NA": could not defind upstream
[lenssnamelist,~]=size(ssnamelist);

%=========================================================
%=======correlation between downstream and upstream=======
%=========================================================
for cs=1:3 %cs: cell set
    temcs=ABC.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    for c=1:lenc
        temc=temcs.cells(c);
        for ss=2:5
            temss=temc.spot_sets(ss);
            useq=temc.sum_intensity(ssupstreami(ss)).intensity;
            useqc=temc.cumulative_intensity(ssupstreami(ss)).intensity;
            [lens,~]=size(temss.spots);
            for s=1:lens
                tems=temss.spots(s);
                seq=tems.intensity;
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).correlation_upstream_intensity=corr(seq,useq);
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).correlation_upstream_cumulative_intensity=corr(seq,useqc);
                xcss=Indigo_xcorr(seq,useq,(-150:150)');
                xccs=Indigo_xcorr(seq,useqc,(-150:150)');
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).cross_correlation_upstream_intensity=xcss;
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).cross_correlation_upstream_cumulative_intensity=xccs;
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).max_cross_correlation_upstream_intensity=max(xcss(~isnan(xcss)));
                ABC.cell_sets(cs).cells(c).spot_sets(ss).spots(s).max_cross_correlation_upstream_cumulative_intensity=max(xccs(~isnan(xccs)));
            end
        end
    end
end
            
        
for cs=1:3 %cs: cell set
    temcs=ABC.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    for c=1:lenc
        temc=temcs.cells(c);
        temcorrelation(4,1).name=[];
        temcorrelation(4,1).correlation_upstream_intensity=[];
        temcorrelation(4,1).correlation_upstream_cumulative_intensity=[];
        temcorrelation(4,1).max_cross_correlation_upstream_intensity=[];
        temcorrelation(4,1).max_cross_correlation_upstream_cumulative_intensity=[];
        for ss=2:5
            temss=temc.spot_sets(ss);
            temcorrelation(ss-1).name=ssnamelist(ss);
            if isempty(temss.spots)
                temcorrelation(ss-1).correlation_upstream_intensity=[];
                temcorrelation(ss-1).correlation_upstream_cumulative_intensity=[];
                temcorrelation(ss-1).max_cross_correlation_upstream_intensity=[];
                temcorrelation(ss-1).max_cross_correlation_upstream_cumulative_intensity=[];
            else
                temcorrelation(ss-1).correlation_upstream_intensity=vertcat(temss.spots(:).correlation_upstream_intensity);
                temcorrelation(ss-1).correlation_upstream_cumulative_intensity=vertcat(temss.spots(:).correlation_upstream_cumulative_intensity);
                temcorrelation(ss-1).max_cross_correlation_upstream_intensity=vertcat(temss.spots(:).max_cross_correlation_upstream_intensity);
                temcorrelation(ss-1).max_cross_correlation_upstream_cumulative_intensity=vertcat(temss.spots(:).max_cross_correlation_upstream_cumulative_intensity);
            end
        end
        ABC.cell_sets(cs).cells(c).correlation=temcorrelation;
    end
end 


for cs=1:3 %cs: cell set
    temcs=ABC.cell_sets(cs);
    temcorrelation(4,1).name=[];
    temcorrelation(4,1).correlation_upstream_intensity=[];
    temcorrelation(4,1).correlation_upstream_cumulative_intensity=[];
    temcorrelation(4,1).max_cross_correlation_upstream_intensity=[];
    temcorrelation(4,1).max_cross_correlation_upstream_cumulative_intensity=[];
    sumcor=vertcat(temcs.cells(:).correlation);
    for ss=2:5
        temss=sumcor(vertcat(sumcor(:).name)==ssnamelist(ss));
        temcorrelation(ss-1).name=ssnamelist(ss);
        temcorrelation(ss-1).correlation_upstream_intensity=vertcat(temss(:).correlation_upstream_intensity);
        temcorrelation(ss-1).correlation_upstream_cumulative_intensity=vertcat(temss(:).correlation_upstream_cumulative_intensity);
        temcorrelation(ss-1).max_cross_correlation_upstream_intensity=vertcat(temss(:).max_cross_correlation_upstream_intensity);
        temcorrelation(ss-1).max_cross_correlation_upstream_cumulative_intensity=vertcat(temss(:).max_cross_correlation_upstream_cumulative_intensity);
    end
    ABC.cell_sets(cs).correlation=temcorrelation;
end 


save('ABC.mat','ABC')
