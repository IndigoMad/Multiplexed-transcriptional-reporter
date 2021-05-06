load('ABC.mat')
load('..\figure_scripts\color.mat')

ssnamelist=["Gene A";"Gene B";"Gene C";"Gene B fast";"Gene B slow"];
ssupstream=["0";"Gene A";"Gene B";"Gene A";"Gene A"];
%"0": upstream strating from time 0
%"NA": could not defind upstream
[lenssnamelist,~]=size(ssnamelist);

%=========================================================
%==========Response time - upstream rate partition========
%=========================================================

%lists of information for gene B partition and locating the spot after
%partition.
totallist.is=[];%the list of the spot index in the spot set.
totallist.ic=[];%the list of the cell index in the cell set.
totallist.ics=[];%the list of the cell set index in the cell sets.
totallist.BRT=[];%the list of response time of gene B.
totallist.BRTUR=[];%the list of response time upstream rate of gene B.
totallist.Bstart=[];%the list of starting time of gene B.
for cs=1:3 %cs: cell set
    temcs=ABC.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    for c=1:lenc%c: cell
        temc=temcs.cells(c);
        %initialization
        ABC.cell_sets(cs).cells(c).spot_sets(4).name='Gene B fast';
        ABC.cell_sets(cs).cells(c).spot_sets(4).spots=[];
        ABC.cell_sets(cs).cells(c).spot_sets(4).description='mCherry:CFP=0:1, the spots in Gene B that response faster';
        ABC.cell_sets(cs).cells(c).spot_sets(5).name='Gene B slow';
        ABC.cell_sets(cs).cells(c).spot_sets(5).spots=[];
        ABC.cell_sets(cs).cells(c).spot_sets(5).description='mCherry:CFP=0:1, the spots in Gene B that response slower';

        %building the sublist in a cell 
        [lens,~]=size(temc.spot_sets(2).spots);
        is=(1:lens)';
        ic=repmat(c,lens,1);
        ics=repmat(cs,lens,1);
        BRT=vertcat(temc.spot_sets(2).spots(:).response_time);
        BRTUR=vertcat(temc.spot_sets(2).spots(:).response_time_upstream_rate);
        Bstart=vertcat(temc.spot_sets(2).spots(:).starting);
        mitosis=temc.mitosis;
        %Response time: excluding response after mitosis, except there is
        %no response before mitosis.

        %the cell did not undergo mitosis
        if isempty(mitosis)
            totallist.is=[totallist.is;is];
            totallist.ic=[totallist.ic;ic];
            totallist.ics=[totallist.ics;ics];
            totallist.BRT=[totallist.BRT;BRT];
            totallist.BRTUR=[totallist.BRTUR;BRTUR];
            totallist.Bstart=[totallist.Bstart;Bstart];
            
        %the cell underwent mitosis
        else
            if min(Bstart)>mitosis
                totallist.is=[totallist.is;is];
                totallist.ic=[totallist.ic;ic];
                totallist.ics=[totallist.ics;ics];
                totallist.BRT=[totallist.BRT;BRT];
                totallist.BRTUR=[totallist.BRTUR;BRTUR];
                totallist.Bstart=[totallist.Bstart;Bstart];
            else
                %Removing response after mitosis
                totallist.is=[totallist.is;is(Bstart<mitosis)];
                totallist.ic=[totallist.ic;ic(Bstart<mitosis)];
                totallist.ics=[totallist.ics;ics(Bstart<mitosis)];
                totallist.BRT=[totallist.BRT;BRT(Bstart<mitosis)];
                totallist.BRTUR=[totallist.BRTUR;BRTUR(Bstart<mitosis)];
                totallist.Bstart=[totallist.Bstart;Bstart(Bstart<mitosis)];
            end
        end
    end
end

%partition
BRT=totallist.BRT;
BRTUR=totallist.BRTUR;

BRTUR=log(BRTUR);% Upstream rate obey exponential distribution, though we 
%have no idea why. logarithm can transform this into a normal-like 
%distribution, in order to facillitate the performance of cluster algorithm. 

% Using Gaussian mixture distribution model with 2 components for
% clustering. 
GMModel=fitgmdist([BRT,BRTUR],2);
species=cluster(GMModel,[BRT,BRTUR]);


%======Check the results, for EM algorithm is local optimum, although======
%======the patition is stable (the results are the same in most time)======
BRT1=BRT(species==1);
BRTUR1=BRTUR(species==1);
BRT2=BRT(species==2);
BRTUR2=BRTUR(species==2);

figure()
[X,Y]=meshgrid(0:max(BRT)/49:max(BRT),min(BRTUR):(max(BRTUR)-min(BRTUR))/49:max(BRTUR));
X=reshape(X,2500,1);
Y=reshape(Y,2500,1);
f=ksdensity([BRT,BRTUR],[X,Y]);
X=reshape(X,50,50);
Y=reshape(Y,50,50);
f=reshape(f,50,50);
pcolor(X,Y,f);
shading interp
hold on;
scatter(BRT1,BRTUR1,'filled','markerfacecolor',colorr,'markeredgecolor','k','linewidth',0.75)
scatter(BRT2,BRTUR2,'filled','markerfacecolor',colorb,'markeredgecolor','k','linewidth',0.75)


[X1,Y1]=meshgrid(0:max(BRT1)/49:max(BRT1),min(BRTUR1):(max(BRTUR1)-min(BRTUR1))/49:max(BRTUR1));
X1=reshape(X1,2500,1);
Y1=reshape(Y1,2500,1);
f1=ksdensity([BRT1,BRTUR1],[X1,Y1]);
X1=reshape(X1,50,50);
Y1=reshape(Y1,50,50);
f1=reshape(f1,50,50);

[X2,Y2]=meshgrid(0:max(BRT2)/49:max(BRT2),min(BRTUR2):(max(BRTUR2)-min(BRTUR2))/49:max(BRTUR2));
X2=reshape(X2,2500,1);
Y2=reshape(Y2,2500,1);
f2=ksdensity([BRT2,BRTUR2],[X2,Y2]);
X2=reshape(X2,50,50);
Y2=reshape(Y2,50,50);
f2=reshape(f2,50,50);

plot(sum(X1.*(f1./sum(f1,2)),2),min(BRTUR1):(max(BRTUR1)-min(BRTUR1))/49:max(BRTUR1),'color',colorr,'linewidth',0.75);
plot(sum(X2.*(f2./sum(f2,2)),2),min(BRTUR2):(max(BRTUR2)-min(BRTUR2))/49:max(BRTUR2),'color',colorb,'linewidth',0.75);
c = colorbar;
c.Label.String = 'Density';
title("A-B")
xlabel("Response time (min)")
ylabel("log(Upstream expression rate)")
set(gca,'linewidth',0.75);
set(gca,'FontSize',7,'Fontname','Arial')
box on;

str="Input 1 if the results are normal, or input 0 to kill this program\n";
yon=input(str);
if yon~=1
    error('You decided to discard this program')
end


%======Building Gene B partition=======
m1=mean(BRT1);
m2=mean(BRT2);
if m1>m2
    fast=2;
    slow=1;
else
    fast=1;
    slow=2;
end

%====Gene B fast====
is=totallist.is(species==fast);
ic=totallist.ic(species==fast);
ics=totallist.ics(species==fast);
[lens,~]=size(is);
for i=1:lens
    [temlens,~]=size(ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(4).spots);
    if temlens==0
        ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(4).spots=ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(2).spots(is(i));
    else
        ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(4).spots(temlens+1,1)=ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(2).spots(is(i));
    end
end

%====Gene B slow====
is=totallist.is(species==slow);
ic=totallist.ic(species==slow);
ics=totallist.ics(species==slow);
[lens,~]=size(is);
for i=1:lens
    [temlens,~]=size(ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(5).spots);
    if temlens==0
        ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(5).spots=ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(2).spots(is(i));
    else
        ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(5).spots(temlens+1,1)=ABC.cell_sets(ics(i)).cells(ic(i)).spot_sets(2).spots(is(i));
    end
end



for cs=1:3 %cs: cell set
    temcs=ABC.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    for c=1:lenc %c: cell
        temc=temcs.cells(c);
        for ss=4:5 %ss: spot sets
            temss=temc.spot_sets(ss);
            if isempty(temss.spots)
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).name=temss.name;
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time=[];
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_rate=[];
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_cumulation=[];
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting=[];
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting_upstream=[];
                ABC.cell_sets(cs).cells(c).time_statistics(ss).name=temss.name;
                ABC.cell_sets(cs).cells(c).time_statistics(ss).response_time=[];
                ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_interval=[];
                ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=[];
            else
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).name=temss.name;
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time=vertcat(temss.spots(:).response_time);
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_rate=vertcat(temss.spots(:).response_time_upstream_rate);
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).response_time_upstream_cumulation=vertcat(temss.spots(:).response_time_upstream_cumulation);
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting=vertcat(temss.spots(:).starting);
                ABC.cell_sets(cs).cells(c).regulation_relationship(ss-1).starting_upstream=vertcat(temss.spots(:).starting_upstream);
                ABC.cell_sets(cs).cells(c).time_statistics(ss).name=temss.name;
                ABC.cell_sets(cs).cells(c).time_statistics(ss).response_time=vertcat(temss.spots(:).response_time);
                ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_interval=vertcat(temss.spots(:).burst_interval);
                ABC.cell_sets(cs).cells(c).time_statistics(ss).burst_duration=vertcat(temss.spots(:).burst_duration);
            end
        end
    end
end


for cs=1:3 %cs: cell set
    temcs=ABC.cell_sets(cs);
    [lenc,~]=size(temcs.cells);
    temcsrr=vertcat(temcs.cells(:).regulation_relationship);
    temcsts=vertcat(temcs.cells(:).time_statistics);
    for ss=4:5 %ss: spot sets
        temss=temcsrr(vertcat(temcsrr(:).name)==ssnamelist(ss));
        ABC.cell_sets(cs).regulation_relationship(ss-1).name=ssnamelist(ss);
        ABC.cell_sets(cs).regulation_relationship(ss-1).response_time=vertcat(temss(:).response_time);
        ABC.cell_sets(cs).regulation_relationship(ss-1).response_time_upstream_rate=vertcat(temss(:).response_time_upstream_rate);
        ABC.cell_sets(cs).regulation_relationship(ss-1).response_time_upstream_cumulation=vertcat(temss(:).response_time_upstream_cumulation);
        ABC.cell_sets(cs).regulation_relationship(ss-1).starting=vertcat(temss(:).starting);
        ABC.cell_sets(cs).regulation_relationship(ss-1).starting_upstream=vertcat(temss(:).starting_upstream);
        
        temss=temcsts(vertcat(temcsts(:).name)==ssnamelist(ss));
        ABC.cell_sets(cs).time_statistics(ss).name=ssnamelist(ss);
        ABC.cell_sets(cs).time_statistics(ss).response_time=vertcat(temss(:).response_time);
        ABC.cell_sets(cs).time_statistics(ss).burst_interval=vertcat(temss(:).burst_interval);
        ABC.cell_sets(cs).time_statistics(ss).burst_duration=vertcat(temss(:).burst_duration);
        
    end
end





save('ABC.mat','ABC')
