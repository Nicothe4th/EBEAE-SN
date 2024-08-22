%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Synthetic Evaluation of Unmixing with linear Model
%
% SLB7 End-members
%
% DUCD
% August/2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
addpath('./sub_func')
addpath('EBEAE');
addpath('NMF-QMV');
%addpath('GraphL');
addpath('gtvMBO');
addpath('gtvMBO/fractional')
%addpath('SeCoDe-main');
%addpath(genpath(pwd));


sSNR=[30 35 40];    
pDensity=[0.01 0.0075 0.005];

% sSNR=[30];     %% Test with 50, 52.5, 55 y 57.5 dB
% pDensity=[0.01];

analysisType=1;     % 0 --> Absorbance HSI, & 1 --> m-FLIM

Nsamples=128;
nCol=Nsamples;
nRow=Nsamples;
N=4;
ModelType=1;
Rep=50;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Paremeters of EBEAE-SN & EBEAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initcond=6;             % Initial condition of end-members matrix: 6 (VCA) and 8 (SISAL).
rho=0.1;               % Similarity weight in end-members estimation
lambda=0.1;             % Entropy weight for abundance estimation
epsilon=1e-3;
maxiter=50;
parallel=0;
downsampling=0.0;       % Downsampling in end-members estimation
display_iter=0;            % Display partial performance in BEAE
lm=0.01;


%%% PISINMF Parameters;

para.dimX= Nsamples;
para.dimY= Nsamples;
para.tven= 5 ;
para.tau=epsilon;
para.maxiter=maxiter;
para.delta=50;
para.mu= 0.005;
para.t=25;
para.alpha=lm;

%%% HU-JSTV

optsJSTV.m=Nsamples;
optsJSTV.n=Nsamples;
optsJSTV.lambda1=0.01; % for total variation
optsJSTV.lambda2=0.5; % for sparse term
optsJSTV.lambda3=0.01; % for joint sparsity i.e. L21 term
optsJSTV.mu1=0.01; % for TV regularization term
optsJSTV.mu2=0.01;  % for L21 regularization term
optsJSTV.iter=maxiter;

%%% GLNMF

para_nmf.tol = epsilon;
para_nmf.itermax = maxiter;
para_nmf.lambda = lm;
para_nmf.mu = 0.005;

%%% NMF-QMV

term ='totalVar';
beta_candidates =   10.^(-5:5);

%%% gtvMBO

para_mbo.method = 'gtvMBO';
para_mbo.tol = epsilon;
para_mbo.m = nCol;
para_mbo.n = nRow;
para_mbo.itermax = maxiter;
para_mbo.dt = 0.01;
para_mbo.lambda = 10^(-4.25);
para_mbo.rho = 10^(-3.75);
para_mbo.gamma = 10^(4);

%%% SeCoDe

lambdaSeCoDe = lm;
alfaSeCoDe=0.1;
betaSeCode=0.1249;
deltaSeCoDe=50;
optSeCoDe = [];
optSeCoDe.Verbose = 0;
optSeCoDe.MaxMainIter = maxiter;
optSeCoDe.rho = 100*lambda + 0.5;
optSeCoDe.sigma = 5;
optSeCoDe.AutoRho = 1;
optSeCoDe.AutoRhoPeriod = 10;
optSeCoDe.AutoSigma = 1;
optSeCoDe.AutoSigmaPeriod = 10;
optSeCoDe.XRelaxParam = 1.8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Performance Metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResultsYh=zeros(length(sSNR),6,Rep);
ResultsVh=zeros(length(sSNR),2,Rep);
ResultsAh=zeros(length(sSNR),6,Rep);
ResultsPh=zeros(length(sSNR),6,Rep);
ResultsTh=zeros(length(sSNR),6,Rep);
ResultsPh2=zeros(length(sSNR),6,Rep);

for index=1:length(sSNR)

    SNR=sSNR(index);
    density=pDensity(index);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Execute Blind Estimation Methodologies 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('Synthetic Datasets');
    if analysisType==0
        disp('Absorbance HSI');
    else
        disp('mFLIM image');
    end
    disp('Blind Unmixing Estimation');
    disp(['SNR =' num2str(SNR) ' dB']);
    disp(['density =' num2str(density) ]);
    disp(['Number of end-members=' num2str(N)]);
    for j=1:Rep 

        if analysisType==0
            [Z,P0,A0,V0]=SphericGaussian_Sparse_Synth(SNR,density,ModelType);
            disp('Absorbance HSI');
        else
            [Z,P0,A0,V0]=mFLIM_Sparse_Synth(N,Nsamples,0.25e-9,SNR,density);
            disp('mFLIM image');
        end
        [L,K]=size(Z);
        disp(['Iteration=' num2str(j)])   

        tic;
        disp('%%%%%%%%%%%%%%%%%%');
        disp('EBEAE-SN');
        paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,display_iter];
        [P1,A1,S1,Zh1,V1,J1]=EBEAESN(Z,N,paramvec);

        Tebeaesn=toc;
        ResultsYh(index,1,j)=norm(Zh1-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,1,j)=norm(V1-V0,'fro')/norm(V0,'fro');
        ResultsAh(index,1,j)=errorabundances(A0,A1);
        ResultsPh(index,1,j)=errorendmembers(P0,P1);
        ResultsPh2(index,1,j)=errorSAM(P0,P1);
        ResultsTh(index,1,j)=Tebeaesn;

        tic;
        disp('EBEAE');
        paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,display_iter];
        [P2,A2,S2,Zh2]=EBEAE(Z,N,paramvec);
        Tebeae=toc;
        ResultsYh(index,2,j)=norm(Zh2-Z,'fro')/norm(Z,'fro');
        ResultsAh(index,2,j)=errorabundances(A0,A2);
        ResultsPh(index,2,j)=errorendmembers(P0,P2);
        ResultsPh2(index,2,j)=errorSAM(P0,P2);
        ResultsTh(index,2,j)=Tebeae;
        
        tic;
        disp('PISINMF');
        tic
        [P33,A33] =  PISINMF(Z,N,para);
        Zh3=P33*A33;
        A33(A33<0)=0;
        P33(P33<0)=0;
        A3=A33./repmat(sum(A33,1),[N 1]);
        P3=P33./repmat(sum(P33,1),[L 1]);
        Tpisinmf=toc;
        ResultsYh(index,3,j)=norm(Zh3-Z,'fro')/norm(Z,'fro');
        ResultsAh(index,3,j)=errorabundances(A0,A3);
        ResultsPh(index,3,j)=errorendmembers(P0,P3);
        ResultsPh2(index,3,j)=errorSAM(P0,P3);
        ResultsTh(index,3,j)=Tpisinmf;
  

        disp('HU using Joint-Sparsity and Total-Variation');
        tic;
        Pp4=VCA(Z,'Endmembers',N,'SNR',1,'verbose','no');
        Pp4=Pp4./sum(Pp4,1);
        [Aa4,V4]=funJSTV(Pp4,Z,optsJSTV);
        A4=Aa4./repmat(sum(Aa4,1),[N 1]);
        P4=Pp4./repmat(sum(Pp4,1),[L 1]);
        Zh4=Pp4*Aa4+V4;
        Tjstv=toc;
        ResultsYh(index,4,j)=norm(Zh4-Z,'fro')/norm(Z,'fro');
        ResultsVh(index,2,j)=norm(V4-V0,'fro')/norm(V0,'fro');
        ResultsAh(index,4,j)=errorabundances(A0,A4);
        ResultsPh(index,4,j)=errorendmembers(P0,P4);
        ResultsPh2(index,4,j)=errorSAM(P0,P4);
        ResultsTh(index,4,j)=Tjstv

        
        disp('NMF-QMV with Total variation');        
        tic;
        Pini = VCA(Z,'Endmembers',N,'SNR',1,'verbose','no');
        Aini = FCLSU(Z, Pini)';
        [beta_best, Pp5, Aa5, results_save] = NMF_QMV(Z, N, beta_candidates, term, Aini);
        P5=Pp5./repmat(sum(Pp5,1),[L,1]);
        A5=Aa5./repmat(sum(Aa5,1),[N 1]);
        Zh5=Pp5*Aa5;
        Tqmv=toc;
        ResultsYh(index,5,j)=norm(Zh5-Z,'fro')/norm(Z,'fro');
        ResultsAh(index,5,j)=errorabundances(A0,A5);
        ResultsPh(index,5,j)=errorendmembers(P0,P5);
        ResultsPh2(index,5,j)=errorSAM(P0,P5);
        ResultsTh(index,5,j)=Tqmv;


        
        disp('Nystrom extension');
        seed = 1;
        rng(seed);
        sigma=5;
        tic
        [V, Sigma] = laplacian_nystrom(Z', 3, floor(0.001*N), sigma, seed);
        disp('gtvMBO');
        para_mbo.V = V;
        para_mbo.S = Sigma;
        [Pp6, Aa6, iter] = unmixing(Z, Pini, Aini, para_mbo);
        P6=Pp6./repmat(sum(Pp6,1),[L,1]);
        A6=Aa6./repmat(sum(Aa6,1),[N 1]);
        Zh6=Pp6*Aa6;
        TgtvMBO = toc;
        ResultsYh(index,6,j)=norm(Zh6-Z,'fro')/norm(Z,'fro');
        ResultsAh(index,6,j)=errorabundances(A0,A6);
        ResultsPh(index,6,j)=errorendmembers(P0,P6);
        ResultsPh2(index,6,j)=errorSAM(P0,P6);
        ResultsTh(index,6,j)=TgtvMBO;
       
        
    end
end


%%
NP=4;
AAnd=[]; PM=[]; EEnd=[];
for jj=1:length(sSNR)
        AAnd=[AAnd; ' & '];
        PM=[PM; ' $\pm$ '];
        EEnd=[EEnd; ' \\'];
end

%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Mean Responses in Performance Metrics')
disp('SNR/Density   EBEAE-SN       EBEAE         PISINMF      HU-JSTV    NMF-QMV    gtvMBO');
disp('%%%%%%%%%%%%%%%');
disp('Error in Output Estimation (%)');
disp([num2str(int8(sSNR')) AAnd num2str(mean(ResultsYh(:,1,:),3),NP) PM  num2str(std(ResultsYh(:,1,:),[],3),NP) ...
    AAnd num2str(mean(ResultsYh(:,2,:),3),NP) PM num2str(std(ResultsYh(:,2,:),[],3),NP) ...
    AAnd num2str(mean(ResultsYh(:,3,:),3),NP) PM num2str(std(ResultsYh(:,3,:),[],3),NP) ...
    AAnd num2str(mean(ResultsYh(:,4,:),3),NP) PM num2str(std(ResultsYh(:,4,:),[],3),NP) ...
    AAnd num2str(mean(ResultsYh(:,5,:),3),NP) PM num2str(std(ResultsYh(:,5,:),[],3),NP) ...
    AAnd num2str(mean(ResultsYh(:,6,:),3),NP) PM num2str(std(ResultsYh(:,6,:),[],3),NP) EEnd]);

disp('%%%%%%%%%%%%%%%');
disp('Error in Abundance Estimation (%)');
disp([num2str(int8(sSNR')) ...
    AAnd num2str(mean(ResultsAh(:,1,:),3),NP) PM  num2str(std(ResultsAh(:,1,:),[],3),NP) ...
    AAnd num2str(mean(ResultsAh(:,2,:),3),NP) PM num2str(std(ResultsAh(:,2,:),[],3),NP) ...
    AAnd num2str(mean(ResultsAh(:,3,:),3),NP) PM num2str(std(ResultsAh(:,3,:),[],3),NP) ...
    AAnd num2str(mean(ResultsAh(:,4,:),3,'omitnan'),NP) PM num2str(std(ResultsAh(:,4,:),[],3,'omitnan'),NP) ...
    AAnd num2str(mean(ResultsAh(:,5,:),3),NP) PM num2str(std(ResultsAh(:,5,:),[],3),NP) ...
    AAnd num2str(mean(ResultsAh(:,6,:),3),NP) PM num2str(std(ResultsAh(:,6,:),[],3),NP) EEnd]);

disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation');
disp([num2str(int8(sSNR')) AAnd num2str(mean(ResultsPh(:,1,:),3),NP) PM  num2str(std(ResultsPh(:,1,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh(:,2,:),3),NP) PM num2str(std(ResultsPh(:,2,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh(:,3,:),3),NP) PM num2str(std(ResultsPh(:,3,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh(:,4,:),3),NP) PM num2str(std(ResultsPh(:,4,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh(:,5,:),3),NP) PM num2str(std(ResultsPh(:,5,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh(:,6,:),3),NP) PM num2str(std(ResultsPh(:,6,:),[],3),NP) EEnd]);

disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation (SAM)');
disp([num2str(int8(sSNR')) AAnd num2str(mean(ResultsPh2(:,1,:),3),NP) PM  num2str(std(ResultsPh2(:,1,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh2(:,2,:),3),NP) PM num2str(std(ResultsPh2(:,2,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh2(:,3,:),3),NP) PM num2str(std(ResultsPh2(:,3,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh2(:,4,:),3),NP) PM num2str(std(ResultsPh2(:,4,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh2(:,5,:),3),NP) PM num2str(std(ResultsPh2(:,5,:),[],3),NP) ...
    AAnd num2str(mean(ResultsPh2(:,6,:),3),NP) PM num2str(std(ResultsPh2(:,6,:),[],3),NP) EEnd]);

disp('%%%%%%%%%%%%%%%');
disp('Computation Time');
disp([num2str(int8(sSNR')) AAnd num2str(mean(ResultsTh(:,1,:),3),NP) PM  num2str(std(ResultsTh(:,1,:),[],3),NP) ...
    AAnd num2str(mean(ResultsTh(:,2,:),3),NP) PM num2str(std(ResultsTh(:,2,:),[],3),NP) ...
    AAnd num2str(mean(ResultsTh(:,3,:),3),NP) PM num2str(std(ResultsTh(:,3,:),[],3),NP) ...
    AAnd num2str(mean(ResultsTh(:,4,:),3),NP) PM num2str(std(ResultsTh(:,4,:),[],3),NP) ...
    AAnd num2str(mean(ResultsTh(:,5,:),3),NP) PM num2str(std(ResultsTh(:,5,:),[],3),NP) ...
    AAnd num2str(mean(ResultsTh(:,6,:),3),NP) PM num2str(std(ResultsTh(:,6,:),[],3),NP) EEnd]);

disp('%%%%%%%%%%%%%%%');
disp('Error in Sparse Noise Estimation (%)');
disp([num2str(int8(sSNR')) AAnd num2str(mean(ResultsVh(:,1,:),3),NP) PM  num2str(std(ResultsVh(:,1,:),[],3),NP) ...
    AAnd num2str(mean(ResultsVh(:,2,:),3),NP) PM num2str(std(ResultsVh(:,2,:),[],3),NP) EEnd]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANOVA Analysis

%%
disp('%%%%%%%%%%%%%%%%')
disp('ANOVA OUTPUT STIMATION')
[n, algs, reps] = size(ResultsYh);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsYh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

% Initialize comparison names for only the first algorithm
comparisonNames = {' vs_EBEAE', '     vs_PISINMF', '     vs_JSTV', '     vs_NMF-QMV', '   vs_gtvMBO'};

% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data = squeeze(ResultsYh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data, [], 'off');
    
    % Perform post-hoc comparisons
    results = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova(index, k) = results(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable = array2table(TabAnova, 'VariableNames', comparisonNames);
pValuesTable.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable);
%%
disp('ANOVA Abundances')
[n, algs, reps] = size(ResultsAh);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova_A = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsAh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_A = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_A = squeeze(ResultsAh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_A, [], 'off');
    
    % Perform post-hoc comparisons
    results_A = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_A(index, k) = results_A(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_A = array2table(TabAnova_A, 'VariableNames', comparisonNames);
pValuesTable_A.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_A);
%%
disp('ANOVA EndMembers error')
[n, algs, reps] = size(ResultsPh);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova_Ph = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsPh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Ph = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Ph = squeeze(ResultsPh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Ph, [], 'off');
    
    % Perform post-hoc comparisons
    results_Ph = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Ph(index, k) = results_Ph(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Ph = array2table(TabAnova_Ph, 'VariableNames', comparisonNames);
pValuesTable_Ph.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Ph);
%%
disp('ANOVA EndMembers SAM rror')
[n, algs, reps] = size(ResultsPh2);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova_Ph2 = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsPh2);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Ph2 = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Ph2 = squeeze(ResultsPh2(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Ph2, [], 'off');
    
    % Perform post-hoc comparisons
    results_Ph2 = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Ph2(index, k) = results_Ph2(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Ph2 = array2table(TabAnova_Ph2, 'VariableNames', comparisonNames);
pValuesTable_Ph2.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Ph2);
%%
disp('ANOVA Computational Time')
[n, algs, reps] = size(ResultsTh);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova_Th = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsPh2);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Th = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Th = squeeze(ResultsTh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Th, [], 'off');
    
    % Perform post-hoc comparisons
    results_Th = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Th(index, k) = results_Th(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Th = array2table(TabAnova_Th, 'VariableNames', comparisonNames);
pValuesTable_Th.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Th);
%% 
disp('ANOVA Sparse Error Estimation')
[n, algs, reps] = size(ResultsVh);
TabAnova = zeros(length(sSNR), 1); % Initialize matrix to store p-values

% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the two algorithms in the current condition
    dataAlg1 = squeeze(ResultsVh(index, 1, :)); % reps x 1 for EBEAE-SN
    dataAlg4 = squeeze(ResultsVh(index, 2, :)); % reps x 1 for PISINMF
    
    % Perform two-sample t-test
    [~, p] = ttest2(dataAlg1, dataAlg4);
    
    % Store the p-value
    TabAnova(index) = p;
end

% Convert TabAnova to a table for better visualization
comparisonNames = {'vs_JSTV'};
pValuesTable = array2table(TabAnova, 'VariableNames', comparisonNames);
pValuesTable.Properties.RowNames = arrayfun(@(x) sprintf('Condition%d', x), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable);