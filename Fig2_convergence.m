clear all; clc;
addpath('EBEAE');

N=4;                % Number of End-members
Nsamples=64;
nCol=Nsamples;
nRow=Nsamples;

EndMembersSynth=1; 
sSNR=[30 35 40];    
pDensity=[0.01 0.0075 0.005];
Rep=50;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEBEAE Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=0.1;                   % Weight on Regularization Term  >= 0
lambda=0.1;                                                                                                                                                            % Weight on Entropy of Abundances \in [0,1)
epsilon=1e-3;            % Threshold for convergence
maxiter=51;              % Maximum number of iteration in alternated least squares approach
downsampling=0.0;       % Downsampling factor to estimate end-members \in [0,1)
parallel=1;              % Parallelization in abundance estimation process
disp_iter=0;          % Display results of iterative optimization
initcond=6;
lm=0.01;

ModelType=5; 
paramvec=[initcond,rho,lambda,lm,epsilon,maxiter,downsampling,parallel,disp_iter];
J1_=zeros(3,50,50);
for index=1:length(sSNR)

    SNR=sSNR(index);
    density=pDensity(index);
    disp('NEBEAE-SN Analysis');
    for j=1:Rep
        [Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
        [P1,A1,D1,S1,Zh1,V1, J1,c1]=NEBEAESN(Z,N,paramvec);
        J1_(index,j,:)=J1;
        c1_(index,j) = c1;
    end
end
%% EBEAE
initcond=6;             % Initial condition of end-members matrix: 6 (VCA) and 8 (SISAL).
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.1;             % Entropy weight for abundance estimation
epsilon=1e-2;
maxiter=50;
parallel=1;
downsampling=0.0;       % Downsampling in end-members estimation
display_iter=0;            % Display partial performance in BEAE
lm=0.05;

disp('EBEAE-SN');
paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,display_iter];
J2_=zeros(3,50,50);

for index=1:length(sSNR)
    SNR=sSNR(index);
    density=pDensity(index);
    disp('EBEAE-SN Analysis');
    for j=1:Rep
        [Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,0);
        [P2,A2,S2,Zh2,V2,J2,c2]=EBEAESN(Z,N,paramvec);
        c2_(index,j) = c2; 
        J2_(index,j,:)=J2;
    end
end
%%
J1_e1=squeeze(J1_(1,:,:));
J1_e2=squeeze(J1_(2,:,:));
J1_e3=squeeze(J1_(3,:,:));


c1_means=mean(c1_,2);
c1_std = std(c1_,0,2);
c2_means=mean(c2_,2);
c2_std = std(c2_,0,2);

J2_e1=squeeze(J2_(1,:,:));
J2_e2=squeeze(J2_(2,:,:));
J2_e3=squeeze(J2_(3,:,:));

J1_r1 = J1_e1(3,:);
J1_r2 = J1_e2(10,:);
J1_r3 = J1_e3(4,:);

J2_r1 = J2_e1(13,:);
J2_r2 = J2_e2(5,:);
J2_r3 = J2_e3(31,:);

figure(1)
clf()
hold on;
plot(J1_r1,'LineWidth',2,'Color','b'); axis tight; grid on
plot(J1_r2,'LineWidth',2,'Color','g');
plot(J1_r3,'LineWidth',2,'Color','k');
plot(J2_r1,'LineWidth',2); axis tight; grid on
plot(J2_r2,'LineWidth',2);
plot(J2_r3,'LineWidth',2);
plot(epsilon*ones(5),'LineWidth',2,'Color','r')
set(gca, 'YScale', 'log');
xlabel('Iteration','FontSize',14,'Interpreter','latex');
ylabel('Converge condition: $\frac{|J^{l}-J^{l+1}|}{J^l}$','FontSize',18,'Interpreter','latex');
legend('NEBEAE-SN with 30 dB Noise','NEBEAE-SN with 35 dB Noise','NEBEAE-SN with 40 dB Noise', ...
    'EBEAE-SN with 30 dB Noise','EBEAE-SN with 35 dB Noise','EBEAE-SN with 40 dB Noise', ...
    'epsilon')
%title('','FontSize',16)
