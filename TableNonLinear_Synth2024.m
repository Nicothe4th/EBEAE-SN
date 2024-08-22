%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Synthetic Evaluation of Unmixing with Multi-linear Model
%
%
% DUCD/JNMC
% August/2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
addpath('EBEAE');
addpath('GraphL');
sSNR=[30 35 40];    
pDensity=[0.01 0.0075 0.005];

N=4;                % Number of End-members
Nsamples=64;
nCol=Nsamples;
nRow=Nsamples;
ModelType=5;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model
EndMembersSynth=1;  % 0--> USGS Spectral Library Version & 1 --> Craneotomy 
Rep=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEBEAE Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=0.1;                   % Weight on Regularization Term  >= 0
lambda=0.1;                                                                                                                                                            % Weight on Entropy of Abundances \in [0,1)
epsilon=1e-3;            % Threshold for convergence
maxiter=50;              % Maximum number of iteration in alternated least squares approach
downsampling=0.0;       % Downsampling factor to estimate end-members \in [0,1)
parallel=0;              % Parallelization in abundance estimation process
disp_iter=0;          % Display results of iterative optimization
initcond=6;
lm=0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Performance Metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResultsYh=zeros(length(sSNR),5,Rep);
ResultsAh=zeros(length(sSNR),5,Rep);
ResultsPh=zeros(length(sSNR),5,Rep);
ResultsTh=zeros(length(sSNR),5,Rep);
ResultsPh2=zeros(length(sSNR),5,Rep);
ResultsDh=zeros(length(sSNR),5,Rep);
ResultsVh=zeros(length(sSNR),Rep);

for index=1:length(sSNR)

    SNR=sSNR(index);
    density=pDensity(index);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Execute Blind Estimation Methodologies 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('Synthetic Datasets');
    if EndMembersSynth==0
        disp('Absorbance HSI');
    else
        disp('mFLIM image');
    end
    disp('Blind Unmixing Estimation');
    disp(['SNR =' num2str(SNR) ' dB']);
    disp(['density =' num2str(density) ]);
    disp(['Number of end-members=' num2str(N)]);
    for j=1:Rep 

        if EndMembersSynth==0
            [Z,P0,A0,V0,D0]=Legendre_Sparse_Synth(SNR,density,ModelType);
            titleL='Synthetic Dataset with USGS Spectral Library';
        elseif EndMembersSynth==1
            [Z,P0,A0,V0,D0]=MatternGaussian_Sparse_Synth(SNR,density,ModelType);
            titleL='Synthetic Dataset with Absorbance Spectral Response of Hb, HbO2, Fat and Water';
        end
        [L,K]=size(Z);
        disp(['Iteration=' num2str(j)])   

        tic;
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('NEBEAE-SN Analysis');
        paramvec=[initcond,rho,lambda,lm,epsilon,maxiter,downsampling,parallel,disp_iter];
        [P1,A1,D1,S1,Zh1,V1,J1]=NEBEAESN(Z,N,paramvec);
        Tnebeaesn=toc;
        
        ResultsYh(index,1,j)=norm(Zh1-Z,'fro')/norm(Z,'fro');
        ResultsDh(index,1,j)=norm(D1-D0,'fro')/norm(D0,'fro');
        ResultsVh(index,j)=norm(V1-V0,'fro')/norm(V0,'fro');
        ResultsAh(index,1,j)=errorabundances(A0,A1);
        ResultsPh(index,1,j)=errorendmembers(P0,P1);
        ResultsPh2(index,1,j)=errorSAM(P0,P1);
        ResultsTh(index,1,j)=Tnebeaesn;

        tic;
        disp('%%%%%%%%%%%%%%%%%%%%');
        disp('NEBEAE Analysis');
        paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,disp_iter];
        [P2,A2,D2,S2,Zh2,J2]=NEBEAE3(Z,N,paramvec);
        Tnebeae=toc;
        ResultsYh(index,2,j)=norm(Zh2-Z,'fro')/norm(Z,'fro');
        ResultsDh(index,2,j)=norm(D2-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,2,j)=errorabundances(A0,A2);
        ResultsPh(index,2,j)=errorendmembers(P0,P2);
        ResultsPh2(index,2,j)=errorSAM(P0,P2);
        ResultsTh(index,2,j)=Tnebeae;
        
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
        ResultsYh(index,3,j)=norm(Zh3-Z,'fro')/norm(Z,'fro');
        ResultsDh(index,3,j)=norm(D3'-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,3,j)=errorabundances(A0,A3);
        ResultsPh(index,3,j)=errorendmembers(P0,P3);
        ResultsPh2(index,3,j)=errorSAM(P0,P3);
        ResultsTh(index,3,j)=Tsmlm;
        
        
        tic;
        disp('%%%%%%%%%%%%%%%%%%');
        disp('Unsupervised MLM');
        [P44,A44,D4,Zh4,~]=unmix(Z,N,maxiter);
        P4=P44./repmat(sum(P44),[L,1]);
        A4=A44./repmat(sum(A44),[N,1]);
        Tumlm=toc;
        ResultsYh(index,4,j)=norm(Zh4-Z,'fro')/norm(Z,'fro');
        ResultsDh(index,4,j)=norm(D4-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,4,j)=errorabundances(A0,A4);
        ResultsPh(index,4,j)=errorendmembers(P0,P4);
        ResultsPh2(index,4,j)=errorSAM(P0,P4);
        ResultsTh(index,4,j)=Tumlm;
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%  Graph Multilinear Model Unmixing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        zo=reshape(Z',Nsamples,Nsamples,L);
        tic;
        disp('%%%%%%%%%%%%%%%%%%');
        disp('Graph-MLM');
        P5o=VCA(Z,'Endmembers',N,'SNR',1,'verbose','no');
        [Zh5, A55, P55, D5]=spxls_gmlm(Z,zo(:,:,1),N,32,maxiter,P5o);
        P5=P55./repmat(sum(P55),[L,1]);
        A5=A55./repmat(sum(A55),[N,1]);
        Tgmlm=toc;
        ResultsYh(index,5,j)=norm(Zh5-Z,'fro')/norm(Z,'fro');
        ResultsDh(index,5,j)=norm(D5'-D0,'fro')/norm(D0,'fro');
        ResultsAh(index,5,j)=errorabundances(A0,A5);
        ResultsPh(index,5,j)=errorendmembers(P0,P5);
        ResultsPh2(index,5,j)=errorSAM(P0,P5);
        ResultsTh(index,5,j)=Tgmlm;
        
    end
end

%%
NP=4;
AAnd=[]; PM=[]; EEnd=[]; slash=[]; mon=[];
for jj=1:length(sSNR)
        AAnd=[AAnd; ' & '];
        PM=[PM; ' \pm '];
        EEnd=[EEnd; ' \\'];
        slash=[slash; '/'];
        mon = [mon; '$'];
end

clc
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Mean Responses in Performance Metrics')
disp('SNR/Density      NEBEAE-SN       NEBEAE         Supervised MLM      Unsupervised MLM             GMLM');
disp('%%%%%%%%%%%%%%%');
disp('Error in Output Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) AAnd  mon num2str(mean(ResultsYh(:,1,:),3),NP) PM  num2str(std(ResultsYh(:,1,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsYh(:,2,:),3),NP) PM num2str(std(ResultsYh(:,2,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsYh(:,3,:),3),NP) PM num2str(std(ResultsYh(:,3,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsYh(:,4,:),3),NP) PM num2str(std(ResultsYh(:,4,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsYh(:,5,:),3),NP) PM num2str(std(ResultsYh(:,5,:),[],3),NP) mon...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Abundance Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) AAnd mon num2str(mean(ResultsAh(:,1,:),3),NP) PM  num2str(std(ResultsAh(:,1,:),[],3),NP)  mon...
    AAnd mon num2str(mean(ResultsAh(:,2,:),3),NP) PM num2str(std(ResultsAh(:,2,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsAh(:,3,:),3),NP) PM num2str(std(ResultsAh(:,3,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsAh(:,4,:),3),NP) PM num2str(std(ResultsAh(:,4,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsAh(:,5,:),3),NP) PM num2str(std(ResultsAh(:,5,:),[],3),NP) mon...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) AAnd mon num2str(mean(ResultsPh(:,1,:),3),NP) PM  num2str(std(ResultsPh(:,1,:),[],3),NP) mon...
    AAnd mon num2str(mean(ResultsPh(:,2,:),3),NP) PM num2str(std(ResultsPh(:,2,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsPh(:,3,:),3),NP) PM num2str(std(ResultsPh(:,3,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsPh(:,4,:),3),NP) PM num2str(std(ResultsPh(:,4,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsPh(:,5,:),3),NP) PM num2str(std(ResultsPh(:,5,:),[],3),NP) mon ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in End-member Estimation (SAM)');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) AAnd mon num2str(mean(ResultsPh2(:,1,:),3),NP) PM  num2str(std(ResultsPh2(:,1,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsPh2(:,2,:),3),NP) PM num2str(std(ResultsPh2(:,2,:),[],3),NP) mon  ...
    AAnd mon num2str(mean(ResultsPh2(:,3,:),3),NP) PM num2str(std(ResultsPh2(:,3,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsPh2(:,4,:),3),NP) PM num2str(std(ResultsPh2(:,4,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsPh2(:,5,:),3),NP) PM num2str(std(ResultsPh2(:,5,:),[],3),NP) mon ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Non-linear interaction levels Estimation');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) ...
    AAnd mon num2str(mean(ResultsDh(:,1,:),3),NP) PM  num2str(std(ResultsDh(:,1,:),[],3),NP) mon  ...
    AAnd mon num2str(mean(ResultsDh(:,2,:),3),NP) PM num2str(std(ResultsDh(:,2,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsDh(:,3,:),3),NP) PM num2str(std(ResultsDh(:,3,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsDh(:,4,:),3),NP) PM num2str(std(ResultsDh(:,4,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsDh(:,5,:),3),NP) PM num2str(std(ResultsDh(:,5,:),[],3),NP) mon ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Computation Time');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) ...
    AAnd mon num2str(mean(ResultsTh(:,1,:),3),NP) PM  num2str(std(ResultsTh(:,1,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsTh(:,2,:),3),NP) PM num2str(std(ResultsTh(:,2,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsTh(:,3,:),3),NP) PM num2str(std(ResultsTh(:,3,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsTh(:,4,:),3),NP) PM num2str(std(ResultsTh(:,4,:),[],3),NP) mon ...
    AAnd mon num2str(mean(ResultsTh(:,5,:),3),NP) PM num2str(std(ResultsTh(:,5,:),[],3),NP) mon ...
    EEnd]);
disp('%%%%%%%%%%%%%%%');
disp('Error in Sparse Noise Estimation (%)');
disp([num2str(int8(sSNR')) slash num2str((pDensity')) ...
    AAnd mon num2str(mean(ResultsVh(:,:),2),NP) PM  num2str(std(ResultsVh(:,:),[],2),NP)  mon EEnd]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
disp('ANOVA OUTPUT STIMATION')
[n, algs, reps] = size(ResultsYh);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsYh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

% Initialize comparison names for only the first algorithm
comparisonNames = {' vs_NEBEAE', '     vs_Supervised MLM', '     vs_Unsupervised MLM', '     vs_GMLM'};

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
disp('ANOVA NonLinear Index Error')
[n, algs, reps] = size(ResultsDh);
numComparisons = algs * (algs - 1) / 2; % Number of pairwise comparisons (6 choose 2)
TabAnova_Dh = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values

[n, algs, reps] = size(ResultsDh);
numComparisons = algs - 1; % Number of comparisons with the first algorithm
TabAnova_Dh = zeros(length(sSNR), numComparisons); % Initialize matrix to store p-values



% Loop over each condition
for index = 1:length(sSNR)
    % Extract data for the current condition
    data_Dh = squeeze(ResultsDh(index, :, :))'; % reps x algs matrix
    
    % Perform one-way ANOVA
    [~, ~, stats] = anova1(data_Dh, [], 'off');
    
    % Perform post-hoc comparisons
    results_Dh = multcompare(stats, 'Display', 'off');
    
    % Store the p-values for pairwise comparisons with the first algorithm
    for k = 1:numComparisons
        TabAnova_Dh(index, k) = results_Dh(k, 6); % 6th column contains p-values
    end
end

% Convert TabAnova to a table for better visualization
pValuesTable_Dh = array2table(TabAnova_Dh, 'VariableNames', comparisonNames);
pValuesTable_Dh.Properties.RowNames = arrayfun(@(x) sprintf('%d dB', sSNR(x)), 1:length(sSNR), 'UniformOutput', false);

% Display the table
disp(pValuesTable_Dh);

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

