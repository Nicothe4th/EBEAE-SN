clear all; 
addpath('EBEAE');
addpath('GraphL');
N=4;                % Number of End-members
load('Urban_F210.mat');
load(['end' num2str(N) '_groundTruth.mat']);

P0=M./sum(M,1);
A0=A./sum(A,1);
Yo=Y(SlectBands,:);
Z=Yo./sum(Yo,1);
[L,K]=size(Z);
Nsamples=nCol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEBEAE Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=0.1;                   % Weight on Regularization Term  >= 0
lambda=0.1;                                                                                                                                                            % Weight on Entropy of Abundances \in [0,1)
epsilon=1e-3;            % Threshold for convergence
maxiter=50;              % Maximum number of iteration in alternated least squares approach
downsampling=0.0;       % Downsampling factor to estimate end-members \in [0,1)
parallel=1;              % Parallelization in abundance estimation process
disp_iter=0;          % Display results of iterative optimization
initcond=6;
lm=0.01;

%%% AMLMPSO Parameters
alpha=0.01; 
beta=0.015; 
gamma=0.1; 

tic;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('NEBEAE-SN Analysis');
paramvec=[initcond,rho,lambda,lm,epsilon,maxiter,downsampling,parallel,disp_iter];
[P1,A1,D1,S1,Zh1,V1, J1]=NEBEAESN(Z,N,paramvec);
Tnebeaesn=toc;

ResultsYh(1)=norm(Zh1-Z,'fro')/norm(Z,'fro');
ResultsAh(1)=errorabundances(A0,A1);
ResultsPh(1)=errorendmembers(P0,P1);
ResultsPh2(1)=errorSAM(P0,P1);
ResultsTh(1)=Tnebeaesn;

tic;
disp('%%%%%%%%%%%%%%%%%%%%');
disp('NEBEAE Analysis');
paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
[P2,A2,D2,S2,Zh2, J2]=NEBEAE3(Z,N,paramvec);
Tnebeae=toc;
ResultsYh(2)=norm(Zh2-Z,'fro')/norm(Z,'fro');
ResultsAh(2)=errorabundances(A0,A2);
ResultsPh(2)=errorendmembers(P0,P2);
ResultsPh2(2)=errorSAM(P0,P2);
ResultsTh(2)=Tnebeae;

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Supervised MLM');
P_min=-100;
options=optimset('fmincon');
options = optimset(options,'Display','off','Algorithm','sqp','MaxFunEvals',1000,'TolFun',epsilon,'TolCon',1e-8,'GradObj','off');
Aa3=zeros(N+1,K); % The first p variables are the abundances, the p+1'th variable contains the P value
Zh3=zeros(L,K);
P3=VCA(Z,'Endmembers',N,'SNR',1,'verbose','no');
P3=P3./sum(P3,1);
Aini=pinv(P3)*Z;
Aini(Aini<0)=0;
Aini=Aini./repmat(sum(Aini),[N,1]);
for i=1:K
    a_ini=[Aini(:,i); 0.0]; % Initialize with linear unmixing results, P=0
    % Sum-to-one applies to abundances, not P. P is restricted to [P_min,1]
    Aa3(:,i) = fmincon(@(a) sum((Z(:,i)-model_MLM(a,P3)).^2), a_ini,[],[],[ones(1,N) 0],1,[zeros(N,1); P_min],ones(N+1,1),[],options);
    Zh3(:,i) = model_MLM(Aa3(:,i),P3);
end
Tsmlm=toc;
P3=P3./repmat(sum(P3),[L,1]);
A3=Aa3(1:N,:);
D3=Aa3(N+1,:);
ResultsYh(3)=norm(Zh3-Z,'fro')/norm(Z,'fro');
ResultsAh(3)=errorabundances(A0,A3);
ResultsPh(3)=errorendmembers(P0,P3);
ResultsPh2(3)=errorSAM(P0,P3);
ResultsTh(3)=Tsmlm;

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('Unsupervised MLM');
[P44,A44,D4,Zh4,~]=unmix(Z,N,maxiter);
P4=P44./repmat(sum(P44),[L,1]);
A4=A44./repmat(sum(A44),[N,1]);
Tumlm=toc;
ResultsYh(4)=norm(Zh4-Z,'fro')/norm(Z,'fro');
ResultsAh(4)=errorabundances(A0,A4);
ResultsPh(4)=errorendmembers(P0,P4);
ResultsPh2(4)=errorSAM(P0,P4);
ResultsTh(4)=Tumlm;

[P6, A6, D6, S6, d6, Zh6] = AMLMPSO(Z, N ,alpha, beta, gamma);
tam=toc;
ResultsYh(5)=norm(Zh6-Z,'fro')/norm(Z,'fro');
ResultsAh(5)=errorabundances(A0,A6);
ResultsPh(5)=errorendmembers(P0,P6);
ResultsPh2(5)=errorSAM(P0,P6);
ResultsTh(5)=tam;
%%
figure(1)
x=(1:20);
plot(J1);
grid on;
ylabel('Error')
xlabel('Iteration')
