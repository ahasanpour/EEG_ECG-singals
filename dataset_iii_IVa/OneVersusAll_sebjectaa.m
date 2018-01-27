
function [ perc ] = csp_Datatest_III_IVa()

clear all;
path='H:\My Computer\D\University Projects\artificial intellgence\DBtest\III_IVa\DB_true_label\aa_h.mat';


%----------------------------------separate train
load(path);
h.Classlabel=mrk.y';
h.TRIG=mrk.pos;
s=cnt;
clear cnt;
jj=1;
sig=[];
for kk=1:28
    if h.Classlabel(kk,1)==1 || h.Classlabel(kk,1)==2
        if length(sig)==0
           tri(jj,1)=1;
       else
            tri(jj,1)=length(sig)+1;
        end
       
            sig=[sig ; s(h.TRIG(kk,1):h.TRIG(kk,1)+499,:)];
            Labels1(jj,1)=h.Classlabel(kk,1);jj=jj+1;
        
    end
    
end

signals=sig;
trigs=tri;
Labels=Labels1; 

%--------------------------------------separate test


jj=1;
sigt=[];
for kk=29:280
    if h.Classlabel(kk,1)==1 || h.Classlabel(kk,1)==2
        if length(sigt)==0
           trit(jj,1)=1;
       else
            trit(jj,1)=length(sigt)+1;
        end
       
            sigt=[sigt ; s(h.TRIG(kk,1):h.TRIG(kk,1)+499,:)];
            Labels1t(jj,1)=h.Classlabel(kk,1);jj=jj+1;
        
    end
    
end

signalst=sigt;
trigst=trit;
Labelst=Labels1t; 
%---------------------------
startPoint=50;
endpoint=250;

[One Two ]=separateClass(signals ,Labels, trigs);
j=1;    
        for h=1:500:length(One)
           
            Cov1=(One(h+startPoint:h+endpoint,:)'*(One(h+startPoint:h+endpoint,:)))/trace((One(h+startPoint:h+endpoint,:)'*(One(h+startPoint:h+endpoint,:))));
            cz1(1,j)=struct('Cov1',{Cov1});
           
            j=j+1;
        end
        j=1; 
        for h=1:500:length(Two)
            
            
             Cov2=(Two(h+startPoint:h+endpoint,:)'*(Two(h+startPoint:h+endpoint,:)))/trace((Two(h+startPoint:h+endpoint,:)'*(Two(h+startPoint:h+endpoint,:))));
            cz2(1,j)=struct('Cov2',{Cov2});
            j=j+1;
        end
        
        co1=cz1(1,1).Cov1;
        co2=cz2(1,1).Cov2;
        
        for h=2:length(cz1)
            co1=co1+cz1(1,h).Cov1;
        end
        for h=2:length(cz2)
            co2=co2+cz2(1,h).Cov2; 
        end 
        co1=co1/length(cz1);
        co2=co2/length(cz2);
        [V D]=eig(co1,co1+co2);
        
        filters=V(:,[1:2 117:118]);
        j1=1; 
        j2=1;
        FVtrain=[];
        FVtest=[];
        for q=1:length(Labels)
                FVtrain(j1,:)=log(var(signals(trigs(q,1)+startPoint:trigs(q,1)+endpoint,:)*filters));
                j1=j1+1;           
        end
         for q=1:length(Labelst)
                FVtest(j2,:)=log(var(signalst(trigst(q,1)+startPoint:trigst(q,1)+endpoint,:)*filters));
                j2=j2+1;
         end
        ld=classify(FVtest,FVtrain,Labels);
        df=0;
        for ij=1:length(Labelst)
            if ld(ij,1)==Labelst(ij,1)
                df=df+1;
            end
        end
        percent=df/length(Labelst);
        perc=percent*100;

end
function [one two ]=separateClass(data,label,trig)
    one=[];
    two=[];
    for r=1:length(trig)
       if label(r,1)==1
           one=[one;data(trig(r,1):trig(r,1)+499,:)];
       else
           two=[two;data(trig(r,1):trig(r,1)+499,:)];
       end
    end
        
end
