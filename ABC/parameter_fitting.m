load('ABC.mat')
load('..\figure_scripts\color.mat')

ssnamelist=["Gene A";"Gene B";"Gene C";"Gene B fast";"Gene B slow"];
ssupstream=["0";"Gene A";"Gene B";"Gene A";"Gene A"];
ssupstreami=[0;1;2;1;1];
%"0": upstream strating from time 0
[lenssnamelist,~]=size(ssnamelist);

% manual_division={
%     [-inf,44,71,97,112,130,145,166,193,213,241,267,295,314,337,352,375,403,428,451,471,494];
%     [-inf,78,123,167,211,268,322,377,409,461,495,542,599,663,705,752,809,853,918,977];
%     [-inf,23,107,165,205,273,318,380,447,495,562,625,667,708,756,802,853,946]};
% manual_start=[1,1,1];

%These results are given by manually choosing peaks from the FFT density.
fft_result=[23.30,61.24,33.07];
fft_b=[155,291,245];
fft_ad=[6,5,8];


% global Ei Di Em

figure

no_tsa=vertcat(ABC.cell_sets(1:3).time_statistics);
for ss=1:5
    ABC.fitting(ss,1).name=ssnamelist(ss);
    burst_duration=vertcat(no_tsa(vertcat(no_tsa(:).name)==ssnamelist(ss)).burst_duration);
    burst_interval=vertcat(no_tsa(vertcat(no_tsa(:).name)==ssnamelist(ss)).burst_interval);
    %fitting koff with exponential distribution. Knowing the sample 
    %frequence is 10min/frame, the turning-off event happened before 10min
    %could not be detected. 
    ABC.fitting(ss).koff=1/(mean(burst_duration)-10);
    %fitting kon with exponential distribution. Knowing the sample 
    %frequence is 10min/frame, the turning-on event happened before 10min
    %could not be detected. 
    ABC.fitting(ss).kon=1/(mean(burst_interval)-10);
    
    if ss==1%calculate kon1 for gene A
        response_time=vertcat(no_tsa(vertcat(no_tsa(:).name)==ssnamelist(ss)).response_time);
        ABC.fitting(ss).kon1=1/(mean(response_time)-10);
    end
    
    if ss<4%Treat all gene B as the same in kT estimation
        intensity=[];
        for cs=1:3%cs: cell set
            temcs=ABC.cell_sets(cs);
            [lenc,~]=size(temcs.cells);
            for c=1:lenc%c: cell
                temc=temcs.cells(c);
                temss=temc.spot_sets(ss);
                [lens,~]=size(temss.spots);
                for s=1:lens
                    tems=temss.spots(s);
                    ti=tems.intensity;
                    %removal of intensity that are less than or equal to 0
                    ti=ti(ti>0);
                    intensity=[intensity;ti];
                end
            end
        end
        subplot(3,3,ss)
        xlimlist=[500,1000,1000];
        x=-xlimlist(ss)/10:xlimlist(ss)/100:xlimlist(ss);
        histogram(intensity,x,'Normalization','pdf')
        %histogram(intensity,x)
        k=20;
        hold on;
        
%         mixture_normal_Poisson(This did not work)
%         ================================================
%         [e_par,h_var]=intensity_fitting_mixture_normal_Poisson(intensity,k);
%         lambda=e_par{1}
%         p=e_par{2}
%         par=e_par{3};
%         u=par(1,:);
%         sd=par(2,:);
%         f=par(3,:);
%         x=1:1000;
%         for i=p+1:k
%             plot(f(i)*normpdf(x,u(i),sd(i)));
%         end
% %         bar(x,p,1,'facealpha',0.5);
        

%         mixture_normal(This did not work)
%         ================================================
%         T=fitgmdist(intensity,k);
%         u=(T.mu)';
%         sd=reshape(T.Sigma,1,k);
%         f=T.ComponentProportion;
%         for i=1:k
%             plot(x,f(i)*normpdf(x,u(i),sd(i)))
%         end


%         Possion, intensity=a*(mRNA count)+b (This did not work)
%         ================================================
%         lenb=101;
%         l=1:lenb;
%         lambdalist=1:lenb;
%         alist=1:lenb;
%         blist=1:lenb;
%         for i=1:lenb
%             i
%             b=i;
%             Ei=mean(intensity);%mean of intensity (i=a*m+b)
%             Di=var(intensity);%variance of intensity (i=a*m+b)
%             %Di=(Ei-b)^2*(1-exp(-lambda)-lambda*exp(-lambda))/lambda
%             syms symlambdaaaa;
%             eq=(Di-((Ei-b)^2)*(1-exp(-symlambdaaaa)-symlambdaaaa*exp(-symlambdaaaa))/symlambdaaaa)==0;
%             lambda=vpasolve(eq,symlambdaaaa);
%             lambisempty(lambda)
%              da=double(lambda);
%             if
%                 l(i)=-inf;
%             else
%                 a=(Ei-b)*(1-exp(-lambda))/lambda;
%                 m=(intensity-b)/a;
%                 [N,~]=histcounts(m,0.5:1:20.5);
%                 p=poisspdf(1:20,lambda);
%                 p=p/sum(p);
%                 l(i)=N*log(p');
%                 alist(i)=a;
%                 blist(i)=b;
%                 lambdalist(i)=lambda;
%             end
%         end
%         l(alist<0)=0;
%         [~,i]=max(l);
%         a=alist(i);
%         b=blist(i);
%         lambda=lambdalist(i);
%         p=poisspdf(1:20,lambda);
%         p=p/a;
%         x=(1:20)*a+b;
%         bar(x,p,1,'facealpha',0.5);
%         xlim([0,xlimlist(ss)])
%         str=sprintf('%s kT=%.2f I=%.2fm+%d',ssnamelist(ss),lambda,a,b);
%         title(str)
%         set(gca,'linewidth',0.75);
%         set(gca,'FontSize',7,'Fontname','Arial')
%         box on;


%         Possion, intensity=a*(mRNA count) (This did not work, solve it in a wrong way?)
%         ================================================
%         Ei=mean(intensity)%mean of intensity (i=a*m)
%         Di=var(intensity)%variance of intensity (i=a*m)
%         %Di=Ei^2*(1-exp(-lambda)-lambda*exp(-lambda))/lambda
%         lambda=fzero(@flambda,Ei);
%         lambda=double(lambda);
%         a=Ei*(1-exp(-lambda))/lambda;
%         m=intensity/a;
%         p=poisspdf(1:20,lambda);
%         p=p/a;
%         x=(1:20)*a;
%         bar(x,p,1,'facealpha',0.5);
%         xlim([0,xlimlist(ss)])
%         str=sprintf('%s kT=%.2f I=%.2fm',ssnamelist(ss),lambda,a);
%         title(str)
%         set(gca,'linewidth',0.75);
%         set(gca,'FontSize',7,'Fontname','Arial')
%         box on;


%         cdf Fourier transform intensity=a*(mRNA count)+b
%         ================================================
        x=(0:1001)';
        density=ksdensity(intensity,x,'Bandwidth',xlimlist(ss)/250);
        base=ksdensity(intensity,x,'Bandwidth',xlimlist(ss)/25);
        [L,~]=size(density);
        plot(x,density-base,'Color',colorr,'linewidth',2)
        a=fft_result(ss);
        b=fft_b(ss);
        br=b-a*fft_ad(ss);
        m=(intensity-br)/a;
        x=0.5:1:31.5;
        x=a*x+br;
        for i=1:16
            patch([x(i*2-1),x(i*2-1),x(i*2),x(i*2)],[-1/xlimlist(ss),4/xlimlist(ss),4/xlimlist(ss),-1/xlimlist(ss)],'red','FaceAlpha',0.3,'EdgeColor','none')
        end
        xlim([-xlimlist(ss)/10,xlimlist(ss)])
        str=sprintf('%s I=%.2fm+%.2f',ssnamelist(ss),a,br);
        title(str)
        xlabel('flourescent intensity (a.u.)')
        ylabel('Density')
        set(gca,'linewidth',0.75);
        set(gca,'FontSize',7,'Fontname','Arial')
        box on;
        
        Y = fft(density-base);
        P2 = abs(Y/L);%the two-sided spectrum
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);%the single-sided spectrum
        f=(0:(L/2))/L;%for plotting
        t=1./f;%for plotting
        subplot(3,3,ss+3)
        plot(t,P1,'Color',colorr,'linewidth',0.7)
        xlim([0,xlimlist(ss)/10])
        ylim([0,0.0001])
        xlabel('Perior (a.u.)')
        ylabel('FFT Density')
        set(gca,'linewidth',0.75);
        set(gca,'FontSize',7,'Fontname','Arial')
        box on;
        
        counts=histcounts(m,0.5:1:30.5);
        counts=counts/sum(counts);
        x=1:30;
%         Em=counts*x';
%         lambda=fzero(@Em_lambda,Em);
        lambdalist=1:0.1:15;
        LLlist=lambdalist;
        for i=1:141
            lambda=lambdalist(i);
            p_p=poisspdf(x,lambda);
            p_p=p_p/sum(p_p);
            LLlist(i)=log(p_p)*counts';
        end

        [~,I]=max(LLlist);
        lambda=lambdalist(I);
        p_p=poisspdf(x,lambda);
        p_p=p_p/sum(p_p);
        subplot(3,3,ss+6)
        hold on;
        bar(x,counts,1,'facecolor',colorb)
        bar(x,p_p,1,'facealpha',0.5,'facecolor',colorr)
        str=sprintf("kT=%.2f",lambda);
        title(str)
        xlabel('mRNA counts')
        ylabel('Density')
        ylim([0,0.2])
        set(gca,'linewidth',0.75);
        set(gca,'FontSize',7,'Fontname','Arial')
        box on;
        
%         cdf mutitaper (This worked but not good enough as fft)
%         ================================================
%         x=(0:1001)';
%         density=ksdensity(intensity,x,'Bandwidth',xlimlist(ss)/250);
%         base=ksdensity(intensity,x,'Bandwidth',xlimlist(ss)/25);
%         plot(x,density-base,'Color',colorr,'linewidth',2)
%         xlim([0,xlimlist(ss)])
%         title(ssnamelist(ss))
%         subplot(3,3,ss+3)
%         pmtm(density-base);
% %         [pxx,f]=pmtm(density-base);
% %         plot(f/pi,log10(pxx./(f/pi)));
% %         ylabel("Thomson¡¯s multitaper power spectral density")
% %         xlabel("frequency")
%         set(gca,'linewidth',0.75);
%         set(gca,'FontSize',7,'Fontname','Arial')
%         box on;
%        

    


%         manual_division(Don't use this method unless you got not solution)
%         ================================================
%         x=(0:1001)';
%         density=ksdensity(intensity,x,'Bandwidth',xlimlist(ss)/250);
%         base=ksdensity(intensity,x,'Bandwidth',xlimlist(ss)/25);
%         [L,~]=size(density);
%         plot(x,density-base,'Color',colorr,'linewidth',2)
%         xlim([0,xlimlist(ss)])
%         title(ssnamelist(ss))
%         xlabel('flourescent intensity (a.u.)')
%         ylabel('Density')
%         set(gca,'linewidth',0.75);
%         set(gca,'FontSize',7,'Fontname','Arial')
%         box on;
%         
%         bands=manual_division{ss};
%         bands=bands(2:end);
%         [~,nb]=size(bands);
%         nb=floor(nb/2);
%         for i=1:nb
%             patch([bands(i*2-1),bands(i*2-1),bands(i*2),bands(i*2)],[-1/xlimlist(ss),4/xlimlist(ss),4/xlimlist(ss),-1/xlimlist(ss)],'red','FaceAlpha',0.3,'EdgeColor','none')
%         end
%         counts=histcounts(intensity,manual_division{ss});
%         [~,nb]=size(counts);
%         counts=counts/sum(counts);
%         x=1:nb;
%         Em=counts*x';
% %         X=1:1000;
% %         plot(X,Em_lambda(X));
%         lambda=fzero(@Em_lambda,Em);
%         p_p=poisspdf(x,lambda);
%         p_p=p_p/sum(p_p);
%         subplot(3,3,ss+3)
%         hold on;
%         bar(x,counts)
%         bar(x,p_p,'facealpha',0.5)
%         str=sprintf("¦Ë=%.2f",lambda);
%         title(str)
%         xlabel('mRNA counts')
%         ylabel('Density')
%         set(gca,'linewidth',0.75);
%         set(gca,'FontSize',7,'Fontname','Arial')
%         box on;
        ABC.fitting(ss).mRNA_a=a;
        ABC.fitting(ss).mRNA_b=br;
        ABC.fitting(ss).kT=lambda/10;
        if ss==2% kT estimation of gene B are applyed to both the partition
            ABC.fitting(4).mRNA_a=a;
            ABC.fitting(4).mRNA_b=br;
            ABC.fitting(4).kT=lambda/10;
            ABC.fitting(5).mRNA_a=a;
            ABC.fitting(5).mRNA_b=br;
            ABC.fitting(5).kT=lambda/10;
        end
    end
    
    if ss>1%Gene A has no input function
        km_list=-5:0.05:-0.5;
        KA_list=1:0.1:10;
        [~,len_km]=size(km_list);
        [~,len_KA]=size(KA_list);
        LL=zeros(len_km,len_KA);
        u_ss=ssupstreami(ss);
        a=ABC.fitting(ss).mRNA_a;
        b=ABC.fitting(ss).mRNA_b;
        n=2;
        %fitting input function with Hill function, modeled by inhomogeous 
        %Poisson process. Making kon1=kon1_max*TF^n/(TF^n+KA^n). kon1_max is 
        %the upper limit of kon1, TF is the concentration of TF, n is the Hill 
        %coefficient, KA is the TF concentration that make kon1 reach the 1/2 
        %of the limit. n is given by literatures and TF is estimated by the 
        %transcription events. The last two parameters would given by fitting.  

        %n is given by literatures
        ABC.fitting(ss).n=2;

        %Presuming every mRNA produce protein at a same rate r, r is given by 
        %literature, equal to 10 protein/mRNA/min. data were collected
        %10min/frame, so counts of protein = (cumulative mRNA counts)*100. In 
        %this estimation, degradation of mRNA and Protein are ignored.
        for cs=1:3
            temcs=ABC.cell_sets(cs);
            [lenc,~]=size(temcs.cells);
            for c=1:lenc
                temc=temcs.cells(c);
                temss=temc.spot_sets(ss);
                [lens,~]=size(temss.spots);
                time=temc.sum_intensity.time;
                U_intensity=temc.sum_intensity.intensity;
                U_starting=temc.sum_intensity.starting;
                U_starting=(U_starting+10)/10;
                m=(U_intensity-b)/a;
                for s=1:lens
                    tems=temss.spots(s);
                    starting=tems.starting;
                    starting=(starting+10)/10;
                    tem_U_intensity=U_intensity(U_starting:starting);
                    c_m=cumsum(tem_U_intensity);
                    p=c_m*100;
                    for i=1:len_km
                        for j=1:len_KA
                            kon1_max=10^km_list(i);
                            KA=10^KA_list(j);
                            
                            kon1=10*kon1_max.*(p.^n)./(p.^n+KA.^n);
                            inte=sum(kon1);
                            if kon1>0
                                LL(i,j)=LL(i,j)+log(kon1(end))-inte;
                            end
                        end
                    end
                end
            end
        end
        LL(LL==0)=-inf;
        [M1,I1]=max(LL);
        [~,I2]=max(M1);
        Mi=I1(I2);
        Mj=I2;
        kon1_max=10^km_list(Mi);
        KA=10^KA_list(Mj);
        ABC.fitting(ss).kon1_max=kon1_max;
        ABC.fitting(ss).KA=KA;
    end
end       
    

save('ABC.mat','ABC')

function y=flambda(x)
    global Ei Di
    y=Di-(Ei^2).*(1-exp(-x)-x.*exp(-x))./x;
end

function y=Em_lambda(x)
    global Em
    y=Em-x./(1-exp(-x));
end


