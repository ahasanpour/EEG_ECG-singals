function [ perc,kappa ] = csp_Datatest_IV_IIa()

%common spatial pattern
pathTrain='/home/ahp/MyMatlabProjects/IV_2a/A2Th.mat';
pathTest='/home/ahp/MyMatlabProjects/IV_2a/A2Eh.mat';

%----------------------------------separate train
load(pathTrain);


jj=1;
sig=[];
for kk=1:288
    if h.Classlabel(kk,1)==1 || h.Classlabel(kk,1)==2 || h.Classlabel(kk,1)==3 || h.Classlabel(kk,1)==4
        if length(sig)==0
           tri(jj,1)=1;
       else
            tri(jj,1)=length(sig)+1;
        end
       
            sig=[sig ; s(h.TRIG(kk,1):h.TRIG(kk,1)+1749,1:22)];
            Labels1(jj,1)=h.Classlabel(kk,1);jj=jj+1;
        
    end
    
end

signals=sig;
trigs=tri;
Labels=Labels1; 

%--------------------------------------separate test
load(pathTest);

jj=1;
sigt=[];
for kk=1:288
    if h.Classlabel(kk,1)==1 || h.Classlabel(kk,1)==2 || h.Classlabel(kk,1)==3 || h.Classlabel(kk,1)==4
        if length(sigt)==0
           trit(jj,1)=1;
       else
            trit(jj,1)=length(sigt)+1;
        end
       
            sigt=[sigt ; s(h.TRIG(kk,1):h.TRIG(kk,1)+1749,1:22)];
            Labels1t(jj,1)=h.Classlabel(kk,1);jj=jj+1;
        
    end
    
end

signalst=sigt;
trigst=trit;
Labelst=Labels1t; 
%---------------------------
startPoint=750;
endpoint=1550;
m=2;
[One Two Three Four One1 Two1 Three1 Four1 ]=separateClass(signals ,Labels, trigs);
j=1;    
        for h=1:1750:length(One)
          
            Cov1=(One(h+startPoint:h+endpoint,:)'*(One(h+startPoint:h+endpoint,:)))/trace((One(h+startPoint:h+endpoint,:)'*(One(h+startPoint:h+endpoint,:))));
            cz1(1,j)=struct('Cov1',{Cov1});
           
            j=j+1;
        end
        j=1;
        for h=1:1750:length(Two)
            
            %Cov2=cov(Two(h+startPoint:h+endpoint,:));
             Cov2=(Two(h+startPoint:h+endpoint,:)'*(Two(h+startPoint:h+endpoint,:)))/trace((Two(h+startPoint:h+endpoint,:)'*(Two(h+startPoint:h+endpoint,:))));
            cz2(1,j)=struct('Cov2',{Cov2});
            j=j+1;
        end
        j=1;
        for h=1:1750:length(Three)
            
            %Cov2=cov(Two(h+startPoint:h+endpoint,:));
             Cov3=(Three(h+startPoint:h+endpoint,:)'*(Three(h+startPoint:h+endpoint,:)))/trace((Three(h+startPoint:h+endpoint,:)'*(Three(h+startPoint:h+endpoint,:))));
            cz3(1,j)=struct('Cov3',{Cov3});
            j=j+1;
        end
        j=1;
        for h=1:1750:length(Four)
            
            %Cov2=cov(Two(h+startPoint:h+endpoint,:));
             Cov4=(Four(h+startPoint:h+endpoint,:)'*(Four(h+startPoint:h+endpoint,:)))/trace((Four(h+startPoint:h+endpoint,:)'*(Four(h+startPoint:h+endpoint,:))));
            cz4(1,j)=struct('Cov4',{Cov4});
            j=j+1;
        end
        
        co1=cz1(1,1).Cov1;
        co2=cz2(1,1).Cov2;
        co3=cz3(1,1).Cov3;
        co4=cz4(1,1).Cov4;
        for h=2:length(cz1)
            co1=co1+cz1(1,h).Cov1;
        end
        for h=2:length(cz2)
            co2=co2+cz2(1,h).Cov2; 
        end 
        for h=2:length(cz3)
            co3=co3+cz3(1,h).Cov3; 
        end 
        for h=2:length(cz4)
            co4=co4+cz4(1,h).Cov4; 
        end 
        
        co1=co1/length(cz1);
        co2=co2/length(cz2);
        co3=co3/length(cz3);
        co4=co4/length(cz4);
        
        [V1 D1]=eig(co1,co2+co3+co4);
        [V2 D2]=eig(co2,co1+co3+co4);
        [V3 D3]=eig(co3,co1+co2+co4);
        [V4 D4]=eig(co4,co1+co2+co3);
        
        filters1=V1(:,[1:m end-m+1:end]);
        filters2=V2(:,[1:m end-m+1:end]);
        filters3=V3(:,[1:m end-m+1:end]);
        filters4=V4(:,[1:m end-m+1:end]);
        
        j1=1; 
        j2=1;
        FVtrain=[];
        FVtest=[];
        for q=1:length(Labels)
                a=log(var(signals(trigs(q,1)+startPoint:trigs(q,1)+endpoint,:)*filters1));
                b=log(var(signals(trigs(q,1)+startPoint:trigs(q,1)+endpoint,:)*filters2));
                c=log(var(signals(trigs(q,1)+startPoint:trigs(q,1)+endpoint,:)*filters3));
                d=log(var(signals(trigs(q,1)+startPoint:trigs(q,1)+endpoint,:)*filters4));
                FVtrain(j1,:)=[a b c d];
                j1=j1+1;           
        end
         for q=1:length(Labelst)
                a1=log(var(signalst(trigst(q,1)+startPoint:trigst(q,1)+endpoint,:)*filters1));
                a2=log(var(signalst(trigst(q,1)+startPoint:trigst(q,1)+endpoint,:)*filters2));
                a3=log(var(signalst(trigst(q,1)+startPoint:trigst(q,1)+endpoint,:)*filters3));
                a4=log(var(signalst(trigst(q,1)+startPoint:trigst(q,1)+endpoint,:)*filters4));
                FVtest(j2,:)=[a1 a2 a3 a4];
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
        kappa=((percent-0.25)/(1-0.25))*100;
        ld(:,2)=Labelst;
end
function [one two three four One1 Two1 Three1 Four1]=separateClass(data,label,trig)
    one=[];
    two=[];
    three=[];
    four=[];
    One1=[];
    oi=1;
    Two1=[];
    oi1=1;
    Three1=[]
    oi2=1;
    Four1=[];
    oi3=1;
    for r=1:length(trig)
       if label(r,1)==1
           one=[one;data(trig(r,1):trig(r,1)+1749,:)];
           One1(:,:,oi)=data(trig(r,1)+750:trig(r,1)+1500,:);
           oi=oi+1;
       elseif label(r,1)==2
           two=[two;data(trig(r,1):trig(r,1)+1749,:)];
            Two1(:,:,oi1)=data(trig(r,1)+750:trig(r,1)+1500,:);
           oi1=oi1+1;
       elseif label(r,1)==3
            three=[three;data(trig(r,1):trig(r,1)+1749,:)];   
            Three1(:,:,oi2)=data(trig(r,1)+750:trig(r,1)+1500,:);
            oi2=oi2+1;
       elseif label(r,1)==4
           four=[four;data(trig(r,1):trig(r,1)+1749,:)];
           Four1(:,:,oi3)=data(trig(r,1)+750:trig(r,1)+1500,:);
           oi3=oi3+1;
       end
    end
        
end


