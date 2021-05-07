function[A,B,C,Ac,Bc,Cc,Apc,Bpc,Cpc,btA,btB,btC,bsA,bsB,bsC]=TRSrandomcontrol_full(kon1A,konA,koffA,ktA,kpA,kdmA,kdpA,kon1B,konB,koffB,ktB,kpB,kdmB,kdpB,kon1C,konC,koffC,ktC,kpC,kdmC,kdpC,T,tau)

T=ceil(T/tau)+1;
kon1A=kon1A*tau;
konA=konA*tau;
koffA=koffA*tau;
ktA=ktA*tau;
kpA=kpA*tau;
kdmA=kdmA*tau;
kdpA=kdpA*tau;
kon1B=kon1B*tau;
konB=konB*tau;
koffB=koffB*tau;
ktB=ktB*tau;
kpB=kpB*tau;
kdmB=kdmB*tau;
kdpB=kdpB*tau;
kon1C=kon1C*tau;
konC=konC*tau;
koffC=koffC*tau;
ktC=ktC*tau;
kpC=kpC*tau;
kdmC=kdmC*tau;
kdpC=kdpC*tau;


A=zeros(T,1);%mRNA produce of gene A
B=zeros(T,1);%mRNA produce of gene B
C=zeros(T,1);%mRNA produce of gene C

Ac=zeros(T,1);%mRNA count of gene A
Bc=zeros(T,1);%mRNA count of gene B
Cc=zeros(T,1);%mRNA count of gene C

Apc=zeros(T,1);%Protein count of gene A
Bpc=zeros(T,1);%Protein count of gene B
Cpc=zeros(T,1);%Protein count of gene C

btA=[];%Burst start time of gene A
btB=[];%Burst start time of gene B
btC=[];%Burst start time of gene C
bsA=zeros(T,1);%Burst state of gene A
bsB=zeros(T,1);%Burst state of gene B
bsC=zeros(T,1);%Burst state of gene C



pon1A=1-exp(-kon1A);
ponA=1-exp(-konA);
pon1B=1-exp(-kon1B);
ponB=1-exp(-konB);
pon1C=1-exp(-kon1C);
ponC=1-exp(-konC);

poffA=1-exp(-koffA);
poffB=1-exp(-koffB);
poffC=1-exp(-koffC);
%state: 
%-1 Additional state
%1 on state
%0 off state
onA= -1;
onB= -1;
onC= -1;

for t=1:T
    %=========================================
    %=================Gene A==================
    %=========================================

    %==========%Additional State========
    if onA==(-1)
        if unifrnd(0,1)<pon1A % switch to on state
            onA=1;
            bsA(t)=1;
            btA=[btA;t];
            mApro=random('Poisson',ktA);%number of mRNA produce
            A(t)=mApro;
        else
            bsA(t)=0;
            mApro=0;
        end
    %==========%off State========
    elseif onA==0
        if unifrnd(0,1)<ponA % switch to on state
            onA=1;
            bsA(t)=1;
            btA=[btA;t];
            mApro=random('Poisson',ktA);%number of mRNA produce
            A(t)=mApro;
        else
            bsA(t)=0;
            mApro=0;
        end
    %==========%on State========
    else
        if unifrnd(0,1)<poffA % switch to on state
            onA=0;
            bsA(t)=0;
            mApro=0;
        else
            bsA(t)=1;
            mApro=random('Poisson',ktA);%number of mRNA produce
            A(t)=mApro;
        end
    end  
    
    %Totol mRNA and protein
    if t==1%the first frame
        Ac(t)=mApro;
    else
        Ac(t)=Ac(t-1)-random('Poisson',Ac(t-1)*kdmA);%degradation of mRNA
        if Ac(t)<0
            Ac(t)=0;
        end
        Ac(t)=Ac(t)+mApro;%production of mRNA

        Apc(t)=Apc(t-1)-random('Poisson',Apc(t-1)*kdpA);%degradation of protein
        if Apc(t)<0
            Apc(t)=0;
        end
        Apc(t)=Apc(t)+random('Poisson',Ac(t-1)*kpA);%production of protein
    end



    %=========================================
    %=================Gene B==================
    %=========================================

    %==========%Additional State========
    if onB==(-1)
        if unifrnd(0,1)<pon1B % switch to on state
            onB=1;
            bsB(t)=1;
            btB=[btB;t];
            mBpro=random('Poisson',ktB);%number of mRNA produce
            B(t)=mBpro;
        else
            bsB(t)=0;
            mBpro=0;
        end
    %==========%off State========
    elseif onB==0
        if unifrnd(0,1)<ponB % switch to on state
            onB=1;
            bsB(t)=1;
            btB=[btB;t];
            mBpro=random('Poisson',ktB);%number of mRNA produce
            B(t)=mBpro;
        else
            bsB(t)=0;
            mBpro=0;
        end
    %==========%on State========
    else
        if unifrnd(0,1)<poffB % switch to on state
            onB=0;
            bsB(t)=0;
            mBpro=0;
        else
            bsB(t)=1;
            mBpro=random('Poisson',ktB);%number of mRNA produce
            B(t)=mBpro;
        end
    end

    
    %Totol mRNA and protein
    if t==1%the first frame
        Bc(t)=mBpro;
    else
        Bc(t)=Bc(t-1)-random('Poisson',Bc(t-1)*kdmB);%degradation of mRNA
        if Bc(t)<0
            Bc(t)=0;
        end
        Bc(t)=Bc(t)+mBpro;%production of mRNA

        Bpc(t)=Bpc(t-1)-random('Poisson',Bpc(t-1)*kdpB);%degradation of protein
        if Bpc(t)<0
            Bpc(t)=0;
        end
        Bpc(t)=Bpc(t)+random('Poisson',Bc(t-1)*kpB);%production of protein
    end




    %=========================================
    %=================Gene C==================
    %=========================================

    %==========%Additional State========
    if onC==(-1)
        if unifrnd(0,1)<pon1C % switch to on state
            onC=1;
            bsC(t)=1;
            btC=[btC;t];
            mCpro=random('Poisson',ktC);%number of mRNA produce
            C(t)=mCpro;
        else
            bsC(t)=0;
            mCpro=0;
        end
    %==========%off State========
    elseif onC==0
        if unifrnd(0,1)<ponC % switch to on state
            onC=1;
            bsC(t)=1;
            btC=[btC;t];
            mCpro=random('Poisson',ktC);%number of mRNA produce
            C(t)=mCpro;
        else
            bsC(t)=0;
            mCpro=0;
        end
    %==========%on State========
    else
        if unifrnd(0,1)<poffC % switch to on state
            onC=0;
            bsC(t)=0;
            mCpro=0;
        else
            bsC(t)=1;
            mCpro=random('Poisson',ktC);%number of mRNA produce
            C(t)=mCpro;
        end
    end
    
    %Totol mRNA and protein
    if t==1%the first frame
        Cc(t)=mCpro;
    else
        Cc(t)=Cc(t-1)-random('Poisson',Cc(t-1)*kdmC);%degradation of mRNA
        if Cc(t)<0
            Cc(t)=0;
        end
        Cc(t)=Cc(t)+mCpro;%production of mRNA

        Cpc(t)=Cpc(t-1)-random('Poisson',Cpc(t-1)*kdpC);%degradation of protein
        if Cpc(t)<0
            Cpc(t)=0;
        end
        Cpc(t)=Cpc(t)+random('Poisson',Cc(t-1)*kpC);%production of protein
    end
end