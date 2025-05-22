clear all; clc; close all;
addpath('./sub_func')
addpath('EBEAE');
addpath('NMF-QMV');

N=4;
load('Urban_F210.mat');
load(['end' num2str(N) '_groundTruth.mat']);

P0=M./sum(M,1);
A0=A./sum(A,1);
Yo=Y(SlectBands,:);
Z=Yo./sum(Yo,1);
[L,K]=size(Z);
Nsamples=nCol;

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

%%% maxvoldual
lambda_mx = 0.015;

tic;
disp('%%%%%%%%%%%%%%%%%%');
disp('EBEAE-SN');
paramvec=[initcond,rho,lambda, lm, epsilon,maxiter,downsampling,parallel,display_iter];
[P1,A1,S1,Zh1,V1,J1]=EBEAESN(Z,N,paramvec);

Tebeaesn=toc;
ResultsYh(1)=norm(Zh1-Z,'fro')/norm(Z,'fro');
ResultsAh(1)=errorabundances(A0,A1);
ResultsPh(1)=errorendmembers(P0,P1);
ResultsPh2(1)=errorSAM(P0,P1);
ResultsTh(1)=Tebeaesn;
tic;
disp('EBEAE');
paramvec=[initcond,rho,lambda,epsilon,maxiter,downsampling,parallel,display_iter];
[P2,A2,S2,Zh2]=EBEAE(Z,N,paramvec);
Tebeae=toc;
ResultsYh(2)=norm(Zh2-Z,'fro')/norm(Z,'fro');
ResultsAh(2)=errorabundances(A0,A2);
ResultsPh(2)=errorendmembers(P0,P2);
ResultsPh2(2)=errorSAM(P0,P2);
ResultsTh(2)=Tebeae;


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
ResultsYh(3)=norm(Zh3-Z,'fro')/norm(Z,'fro');
ResultsAh(3)=errorabundances(A0,A3);
ResultsPh(3)=errorendmembers(P0,P3);
ResultsPh2(3)=errorSAM(P0,P3);
ResultsTh(3)=Tpisinmf;


disp('HU using Joint-Sparsity and Total-Variation');
tic;
Pp4=VCA(Z,'Endmembers',N,'SNR',1,'verbose','no');
Pp4=Pp4./sum(Pp4,1);
[Aa4,V4]=funJSTV(Pp4,Z,optsJSTV);
A4=Aa4./repmat(sum(Aa4,1),[N 1]);
P4=Pp4./repmat(sum(Pp4,1),[L 1]);
Zh4=Pp4*Aa4+V4;
Tjstv=toc;
ResultsYh(4)=norm(Zh4-Z,'fro')/norm(Z,'fro');
ResultsAh(4)=errorabundances(A0,A4);
ResultsPh(4)=errorendmembers(P0,P4);

ResultsPh2(4)=errorSAM(P0,P4);
ResultsTh(4)=Tjstv;


disp('maxvoldual')
tic
[P7, A7, Zh7] = maxvoldual(Z,N,lambda_mx);
Tmx=toc;
ResultsYh(5)=norm((Zh7./sum(Zh7,1) )-Z,'fro')/norm(Z,'fro');
ResultsAh(5)=errorabundances(A0,A7);
ResultsPh(5)=errorendmembers(P0,P7);
ResultsPh2(5)=errorSAM(P0,P7);
ResultsTh(5)=Tmx;

%%
