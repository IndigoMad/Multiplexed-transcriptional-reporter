
load('..\..\..\data_reconstruction\ABC\ABC.mat')
load('..\..\..\figure_scripts\color.mat')
ssnamelist=["Gene A","Gene B","Gene C"];




response_time{5,3}=0;
burst_interval{5,3}=0;
burst_duration{5,3}=0;
for cs=1:5
    for ss=1:3
        response_time{cs,ss}=ABC.cell_sets(cs).time_statistics(ss).response_time;
        burst_duration{cs,ss}=ABC.cell_sets(cs).time_statistics(ss).burst_duration;
        burst_interval{cs,ss}=ABC.cell_sets(cs).time_statistics(ss).burst_interval;
    end
end




%%
%=======================================================================
%Response time, burst duration and interval distribution of Gene A, B, C
%=======================================================================
figure('Position',[0,0,550,200])
sgtitle("Response time, burst duration and interval distribution of Gene A, B, C",'FontSize',9,'FontName','Arial','FontWeight','bold')
subplot(131)
[faclist]=Indigo_stacking_his([response_time{1,1};response_time{2,1};response_time{3,1}],[response_time{1,2};response_time{2,2};response_time{3,2}],[response_time{1,3};response_time{2,3};response_time{3,3}],'Color',colorlist,'LineWidth',0.75,'Edges',5:100:3005,'Normalization','pdf');
lmd=1/(mean([response_time{1,1};response_time{2,1};response_time{3,1}])-10);
x=0:3005;
y=lmd*exp(-lmd*x);
hold on;
plot(x,y*faclist(3)+2,'Color',colory,'LineWidth',1);
str=sprintf("kon1=%.4f/min",lmd);
text(0.2,0.95,str,'FontSize',8,'Units','normalized','FontName','Arial')
xlabel("Response time (min)")
xlim([0,3005])
set(gca,'linewidth',0.75);
set(gca,'FontSize',7,'Fontname','Arial')
box on;

subplot(132)
[faclist]=Indigo_stacking_his([burst_duration{1,1};burst_duration{2,1};burst_duration{3,1}],[burst_duration{1,2};burst_duration{2,2};burst_duration{3,2}],[burst_duration{1,3};burst_duration{2,3};burst_duration{3,3}],'Color',colorlist,'LineWidth',0.75,'Edges',5:10:305,'Normalization','pdf');
hold on;
for i=1:3
    lmd=1/(mean([burst_duration{1,i};burst_duration{2,i};burst_duration{3,i}])-10);
    x=0:305;
    y=lmd*exp(-lmd*x);
    plot(x,y*faclist(4-i)+3-i,'Color',colory,'LineWidth',1);
    str=sprintf("kon=%.4f/min",lmd);
    text(0.2,0.95-(i-1)/3,str,'FontSize',8,'Units','normalized','FontName','Arial')
end
xlabel("Burst duration (min)")
xlim([0,305])
set(gca,'linewidth',0.75);
set(gca,'FontSize',7,'Fontname','Arial')
box on;

subplot(133)
[faclist]=Indigo_stacking_his([burst_interval{1,1};burst_interval{2,1};burst_interval{3,1}],[burst_interval{1,2};burst_interval{2,2};burst_interval{3,2}],[burst_interval{1,3};burst_interval{2,3};burst_interval{3,3}],'Color',colorlist,'LineWidth',0.75,'Edges',5:10:505,'Normalization','pdf');
hold on;
for i=1:3
    lmd=1/(mean([burst_interval{1,i};burst_interval{2,i};burst_interval{3,i}])-10);
    x=0:505;
    y=lmd*exp(-lmd*x);
    plot(x,y*faclist(4-i)+3-i,'Color',colory,'LineWidth',1);
    str=sprintf("koff=%.4f/min",lmd);
    text(0.2,0.95-(i-1)/3,str,'FontSize',8,'Units','normalized','FontName','Arial')
end
xlabel("Burst interval (min)")
xlim([0,505])
set(gca,'linewidth',0.75);
set(gca,'FontSize',7,'Fontname','Arial')
box on;




%%
%=======================================================================
%Response time, burst duration and interval distribution of Gene A, B, C
%administrated by different DOX concentrations
%=======================================================================
flist=["A","B","C"];
for f=1:3
    figure('Position',[0,0,550,200])
    str=sprintf("Response time, burst duration and interval distribution\nof gene %s administrated with different DOX concentrations",flist(f));
    sgtitle(str,'FontSize',9,'FontName','Arial','FontWeight','bold')
    color=[colorb;colorp;colorr];

    subplot(131)
    [faclist]=Indigo_stacking_his(response_time{1,f},response_time{2,f},response_time{3,f},'Color',color,'LineWidth',0.75,'Edges',5:100:3005,'Normalization','pdf');
    hold on;
    for i=1:3
        lmd=1/(mean([response_time{i,f}])-10);
        x=0:3005;
        y=lmd*exp(-lmd*x);
        plot(x,y*faclist(4-i)+3-i,'Color',colory,'LineWidth',1);
        str=sprintf("kon=%.4f/min",lmd);
        text(0.2,0.95-(i-1)/3,str,'FontSize',8,'Units','normalized','FontName','Arial')
    end
    ks_p_list(3,1)=0;
    [~,ks_p_list(1)]=kstest2(response_time{1,f},response_time{2,f});
    [~,ks_p_list(2)]=kstest2(response_time{1,f},response_time{3,f});
    [~,ks_p_list(3)]=kstest2(response_time{2,f},response_time{3,f});
    str=sprintf("KS-tset\np value\n%s\n%s\n%s",Indigo_double2str(ks_p_list(1)),Indigo_double2str(ks_p_list(2)),Indigo_double2str(ks_p_list(3)));
    text(1,0.5,str,'FontSize',8,'Units','normalized','FontName','Arial')

    xlabel("Response time (min)")
    xlim([0,3005])
    set(gca,'linewidth',0.75);
    set(gca,'FontSize',7,'Fontname','Arial')
    box on;

    subplot(132)
    [faclist]=Indigo_stacking_his(burst_duration{1,f},burst_duration{2,f},burst_duration{3,f},'Color',color,'LineWidth',0.75,'Edges',5:10:305,'Normalization','pdf');
    hold on;
    for i=1:3
        lmd=1/(mean([burst_duration{i,f}])-10);
        x=0:305;
        y=lmd*exp(-lmd*x);
        plot(x,y*faclist(4-i)+3-i,'Color',colory,'LineWidth',1);
        str=sprintf("kon=%.4f/min",lmd);
        text(0.2,0.95-(i-1)/3,str,'FontSize',8,'Units','normalized','FontName','Arial')
    end
    ks_p_list(3,1)=0;
    [~,ks_p_list(1)]=kstest2(burst_duration{1,f},burst_duration{2,f});
    [~,ks_p_list(2)]=kstest2(burst_duration{1,f},burst_duration{3,f});
    [~,ks_p_list(3)]=kstest2(burst_duration{2,f},burst_duration{3,f});
    str=sprintf("KS-tset\np value\n%s\n%s\n%s",Indigo_double2str(ks_p_list(1)),Indigo_double2str(ks_p_list(2)),Indigo_double2str(ks_p_list(3)));
    text(1,0.5,str,'FontSize',8,'Units','normalized','FontName','Arial')

    xlabel("Burst duration (min)")
    xlim([0,305])
    set(gca,'linewidth',0.75);
    set(gca,'FontSize',7,'Fontname','Arial')
    box on;

    subplot(133)
    [faclist]=Indigo_stacking_his(burst_interval{1,f},burst_interval{2,f},burst_interval{3,f},'Color',color,'LineWidth',0.75,'Edges',5:10:505,'Normalization','pdf');
    hold on;
    for i=1:3
        lmd=1/(mean([burst_interval{i,f}])-10);
        x=0:505;
        y=lmd*exp(-lmd*x);
        plot(x,y*faclist(4-i)+3-i,'Color',colory,'LineWidth',1);
        str=sprintf("koff=%.4f/min",lmd);
        text(0.2,0.95-(i-1)/3,str,'FontSize',8,'Units','normalized','FontName','Arial')
    end
    ks_p_list(3,1)=0;
    [~,ks_p_list(1)]=kstest2(burst_interval{1,f},burst_interval{2,f});
    [~,ks_p_list(2)]=kstest2(burst_interval{1,f},burst_interval{3,f});
    [~,ks_p_list(3)]=kstest2(burst_interval{2,f},burst_interval{3,f});
    str=sprintf("KS-tset\np value\n%s\n%s\n%s",Indigo_double2str(ks_p_list(1)),Indigo_double2str(ks_p_list(2)),Indigo_double2str(ks_p_list(3)));
    text(1,0.5,str,'FontSize',8,'Units','normalized','FontName','Arial')

    xlabel("Burst interval (min)")
    xlim([0,505])
    set(gca,'linewidth',0.75);
    set(gca,'FontSize',7,'Fontname','Arial')
    box on;
end




%%
%=======================================================================
%Comparing of response time, burst duration and interval distribution of
%Gene A, B, C between cells with/without TSA 
%=======================================================================
flist=["A","B","C"];
for f=1:3
    figure('Position',[0,0,550,350])
    str=sprintf("Response time, burst duration and interval distribution\nof gene %s with/without TSA",flist(f));
    sgtitle(str,'FontSize',9,'FontName','Arial','FontWeight','bold')
    color=[colorr;colorb];
    indexlist=[1,3;4,5];
    for d=1:2
        subplot(2,3,d*3-2)
        [faclist]=Indigo_stacking_his(response_time{indexlist(1,d),f},response_time{indexlist(2,d),f},'Color',color,'LineWidth',0.75,'Edges',5:100:3005,'Normalization','pdf');
        hold on;
        for i=1:2
            lmd=1/(mean(response_time{indexlist(i,d),f})-10);
            x=0:3005;
            y=lmd*exp(-lmd*x);
            plot(x,y*faclist(3-i)+2-i,'Color',colory,'LineWidth',1);
            str=sprintf("kon1=%.4f/min",lmd);
            text(0.2,0.9-(i-1)/2,str,'FontSize',8,'Units','normalized','FontName','Arial')
        end
        [~,p]=kstest2(response_time{indexlist(1,d),f},response_time{indexlist(2,d),f});
        str=sprintf("KS-tset p: %s",Indigo_double2str(p));
        title(str,'FontSize',8,'FontName','Arial')

        xlabel("Response time (min)")
        xlim([0,3005])
        set(gca,'linewidth',0.75);
        set(gca,'FontSize',7,'Fontname','Arial')
        box on;

        subplot(2,3,d*3-1)
        [faclist]=Indigo_stacking_his(burst_duration{indexlist(1,d),f},burst_duration{indexlist(2,d),f},'Color',color,'LineWidth',0.75,'Edges',5:10:305,'Normalization','pdf');
        hold on;
        for i=1:2
            lmd=1/(mean(burst_duration{indexlist(i,d),f})-10);
            x=0:305;
            y=lmd*exp(-lmd*x);
            plot(x,y*faclist(3-i)+2-i,'Color',colory,'LineWidth',1);
            str=sprintf("kon=%.4f/min",lmd);
            text(0.2,0.9-(i-1)/2,str,'FontSize',8,'Units','normalized','FontName','Arial')
        end
        [~,p]=kstest2(burst_duration{indexlist(1,d),f},burst_duration{indexlist(2,d),f});
        str=sprintf("KS-tset p: %s",Indigo_double2str(p));
        title(str,'FontSize',8,'FontName','Arial')

        xlabel("Burst duration (min)")
        xlim([0,305])
        set(gca,'linewidth',0.75);
        set(gca,'FontSize',7,'Fontname','Arial')
        box on;

        subplot(2,3,d*3)
        [faclist]=Indigo_stacking_his(burst_interval{indexlist(1,d),f},burst_interval{indexlist(2,d),f},'Color',color,'LineWidth',0.75,'Edges',5:10:505,'Normalization','pdf');
        hold on;
        for i=1:2
            lmd=1/(mean(burst_interval{indexlist(i,d),f})-10);
            x=0:505;
            y=lmd*exp(-lmd*x);
            plot(x,y*faclist(3-i)+2-i,'Color',colory,'LineWidth',1);
            str=sprintf("koff=%.4f/min",lmd);
            text(0.2,0.9-(i-1)/2,str,'FontSize',8,'Units','normalized','FontName','Arial')
        end
        [~,p]=kstest2(burst_interval{indexlist(1,d),f},burst_interval{indexlist(2,d),f});
        str=sprintf("KS-tset p: %s",Indigo_double2str(p));
        title(str,'FontSize',8,'FontName','Arial')

        xlabel("Burst interval (min)")
        xlim([0,505])
        set(gca,'linewidth',0.75);
        set(gca,'FontSize',7,'Fontname','Arial')
        box on;
    end
end


%%
%=======================================================================
%
%=======================================================================
