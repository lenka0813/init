X=xlsread('20130329');


Bi=1;
Bj=0;
B1=zeros(50,1964)


% for i=1:24550
%   for j=1:4
%   Bj=Bj+1;
%   B1(Bi,Bj)=X(i,j);
%   end
%   if mod(i,491)==0
%   Bi=Bi+1
%   Bj=0
%   end
% end   

%按照变量展开
B=zscore(B1); 
%对展开后的矩阵进行标准化

B2=mat2cell(B,ones(1,1)*50,ones(491,1)*4);  %把矩阵B进行分块

B3=B2'; %将数据矩阵转置后得到按变量运行的结果

B4=cell2mat(B3);



[coeff,score,latent,T2]=princomp(B4);
explained=100*latent/sum(latent);
[m,n]=size(B4);


result1=cell(n+1,3);
result1(1,:)={'特征值','贡献率','累计贡献率'};
result1(2:end,1)=num2cell(latent);
result1(2:end,2)=num2cell(explained);
result1(2:end,3)=num2cell(cumsum(explained));
summ=cumsum(explained);

%the num for the data axis
for num=1:m 
  if summ(num)<85&&summ(num+1)>85 
  t=num+1;
  elseif summ(1)>85
  t=1;
  end
end %Attempted to access summ(5); index out of bounds because numel(summ)=4.

result=cell2mat(result1(2:5,1:2));

T=score(:,1:t);
T1=mat2cell(T,ones(491,1)*50,ones(1,1)*t);
T11=T1';
Cov=cell(491,1);
COV=cell(491,1);
for k=0:490
k=k+1
[S]=cov(T1{k,1});
Cov{k,1}=S;
COV(k,1)={inv(Cov{k,1})};
end  %计算出每个时刻的协方差矩阵及协方差矩阵的逆

Ti2=cell(50,491);
COVT=COV';
T12=cell2mat(T11);
T13=mat2cell(T12,ones(50,1),ones(491,1)*3);

for h=1:50
  for k=1:491
   Ti2(h,k)={[T13{h,k}]*[COVT{1,k}]*[T13{h,k}]'};
  end
end  %计算出每个时刻的T2值


SPE=cell(1,491);
P=coeff(:,1:t);
C=eye(size(B4,2))-P*P';


B5=mat2cell(B,ones(50,1),ones(491,1)*4);
for h=1:50
  for k=1:491
    SPE(h,k)={[B5{h,k}]*C*C'*[B5{h,k}]'};
  end
end

SPEa=cell2mat(SPE);
M=mean(SPEa,1);
N=std(SPEa,0,1);
N1=N.^2;
g=N1./[2*M];
h=[2*(M.^2)]./N1;


% SPE statistic limit 95% ,99%
pk1=chi2inv(0.95,h());
SPEk=g.*pk1;
pk2=chi2inv(0.99,h());
SPEk1=g.*pk2; 

I=50;
Fa1=finv(0.95,t,I-t);
Talpha=t*(I^2-1)*Fa1/(I*(I-t));
Fa2=finv(0.99,t,m-t);
Talpha1=t*(I^2-1)*Fa2/(I*(I-t));

ti2=cell2mat(Ti2);

t11=t1';
t12=cell2mat(t11);
plot(t12(:,1),t12(:,2));

figure(1);
plot(ti2(7,:)');
hold on;
TS1=Talpha*ones(1,size(ti2,2));
plot(TS1,'g--');
hold on;
TS2=Talpha1*ones(1,size(ti2,2));
plot(TS2,'r--');
title('T2 statistic of Batches ');
hold off;

figure(2);
plot(SPEa(7,:)');
hold on;
TS1=SPEk;
plot(TS1,'g--');
hold on;
TS2=SPEk1;
plot(TS2,'k--');
title('SPE statistic of Batches ');
hold off;
B6=mean(B1,1);
B7=std(B1,0,1);
XT=xlsread('gz26');
XT1=mat2cell(XT,ones(491,1),ones(1,1)*4);
XT2=XT1';
XT3=cell2mat(XT2);
XT4=[XT3-B6]./B7;
XT5=mat2cell(XT4,ones(1,1),ones(491,1)*4);
t1=cell(1,491);
for k=1:491
t1(1,k)={[XT5{1,k}]*P};
end

t2score=cell(1,491);
for k=1:491
t2score(1,k)={[t1{1,k}]*[COVT{1,k}]*[t1{1,k}]'};
end
t2score1=cell2mat(t2score);

spe=cell(1,491);
for k=1:491
spe(1,k)={[XT5{1,k}]*C*C'*[XT5{1,k}]'};
end
spe1=cell2mat(spe);

figure(3);
pareto(result(:,1));
hold off;



figure(4);
plot(t2score1,'b-*');
hold on;
TS1=Talpha*ones(size(XT,1),1);
plot(TS1,'g--');
hold on;
TS2=Talpha1*ones(size(XT,1),1);
plot(TS2,'r--');
M1=310;
M2=491;
axes('position',[0.4,0.65,0.2,0.2]);
hold on;
plot((M1:M2),t2score1(M1:M2),'b-+');
plot((M1:M2),TS1(M1:M2),'g--');
plot((M1:M2),TS2(M1:M2),'r--');
hold off;

figure(5);
plot(spe1,'b-*');
hold on;
TS1=SPEk;
plot(TS1,'g--');
hold on;
TS2=SPEk1;
plot(TS2,'r--');
hold off;


% the uncommon batich
D=xlsread('gz27');
D1=mat2cell(D,ones(491,1),ones(1,1)*4);
D2=D1';
D3=cell2mat(D2);
D4=[D3-B6]./B7;
D5=mat2cell(D4,ones(1,1),ones(491,1)*4);
tD=cell(1,491);
for k=1:491
tD(1,k)={[D5{1,k}]*P};
end


t2Dscore=cell(1,491);
for k=1:491
t2Dscore(1,k)={[tD{1,k}]*[COVT{1,k}]*[tD{1,k}]'};
end
t2Dscore1=cell2mat(t2Dscore);

speD=cell(1,491);
for k=1:491
speD(1,k)={[D5{1,k}]*C*C'*[D5{1,k}]'};
end
speD1=cell2mat(speD);

K=xlsread('gz28');
K1=mat2cell(K,ones(491,1),ones(1,1)*4);
K2=K1';
K3=cell2mat(K2);
K4=[K3-B6]./B7;
K5=mat2cell(K4,ones(1,1),ones(491,1)*4);
tK=cell(1,491);
for k=1:491
tK(1,k)={[K5{1,k}]*P};
end


t2Kscore=cell(1,491);
for k=1:491
t2Kscore(1,k)={[tK{1,k}]*[COVT{1,k}]*[tK{1,k}]'};
end
t2Kscore1=cell2mat(t2Kscore);

speK=cell(1,491);
for k=1:491
speK(1,k)={[K5{1,k}]*C*C'*[K5{1,k}]'};
end
speK1=cell2mat(speK);


L1=cell(1,491);
for k=1:491
L1(1,k)={P'*diag(XT5{1,k})};
end 

L11=cell(1,491);
for k=1:491
L11(1,k)={[L1{1,k}]'};
end

L2=cell(1,491);
for k=1:491
L2(1,k)={[L11{1,k}]*[COVT{1,k}]*[t1{1,k}]'};
end

L3=cell(1,491);
for k=1:491
L3(1,k)={[[XT5{1,k}]*C].^2};
end

figure(7);
bar(L2{1,79});
hold off;



figure(6);
plot(t2Kscore1,'b-*');
hold on;
TS1=Talpha*ones(size(K,1),1);
plot(TS1,'g--');
hold on;
TS2=Talpha1*ones(size(K,1),1);
plot(TS2,'r--');

hold on;
M1=50;
M2=491;
axes('position',[0.4,0.65,0.2,0.2]);
hold on;
plot((M1:M2),t2Dscore1(M1:M2),'b-+');
plot((M1:M2),TS1(M1:M2),'g--');
plot((M1:M2),TS2(M1:M2),'r--');
hold off;

figure(7);
plot(speK1,'b-*');
hold on;
TS1=SPEk;
plot(TS1,'g--');
hold on;
TS2=SPEk1;
plot(TS2,'r--');
hold off;

