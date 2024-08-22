function [A,S]=funJSTV(M,Y,opts)

%This program assume following noise model: 
%Y=MA+S+G  and solve the following problem
% min_{A,S} ||Y-MA-S||_F^2 + lam1*||DhA'||_2,1 + lam1*||DvA'||_1
% +lam2||S||_1 + lam3||A||_{2,1}

% M: mixing matrix; 
% Y: Noisy Image ; 
% A: Abundance matrix
% S : sparse noise; Dh and Dv are Total variation operators.

% m: number of rows in one band.
% n: number of columns in one band.
% b: number of bands.
% e: total number of endmemebrs available.
% k: number of endmembers used to make the image.


% Refrence paper :
% Hyperspectral Unmixing in the Presence of Mixed Noise using Joint-Sparsity and Total-Variation 
% Hemant Kumar Aggarwal, Angshul Majumdar
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing (JSTARS), 2016
% Link : http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7414394
% paper PDF is available from author's homepage as well.


lambda1=opts.lambda1; 
lambda2=opts.lambda2; lambda3=opts.lambda3;
mu1=opts.mu1; mu2=opts.mu2;iter=opts.iter; m=opts.m; n=opts.n;
[~,e]=size(M); [~,p]=size(Y);

Dh=TVmatrix(m,n,'H'); %horizontal total variation
Dv=TVmatrix(m,n,'V'); %vertical total variation
D=Dh'*Dh+Dv'*Dv;  
MtM=M'*M;
A=zeros(e,p);
B1=zeros(p,e);B2=B1; Y1=Y; B3=zeros(e,p);

%% main iteration
for i =1:iter
     for j=1:5
         P=softTh(Dh*(A')+B1,lambda1/mu1);
         Q=softTh(Dv*(A')+B2,lambda1/mu1);
         R=softThL21(A+B3,lambda3/mu2);
         S=softTh(Y1-M*A,1/lambda2);
     
        RHS=(mu1*(P-B1)')*Dh+ (mu1*(Q-B2)')*Dv+ M'*(Y1-S) + mu2*(R-B3);
      
        [a,~]=pcg(@afun,RHS(:),1e-15,5);
      
        A=reshape(a,e,p);       
        
        B1=B1+Dh*(A')-P; 
        B2=B2+Dv*(A')-Q;
        B3=B3+A-R;
        
        A=max(0,A);
        %S=max(0,S);
    end
    Y1=Y1+Y-M*A-S;  
end

 function y = afun(x)      
       X=reshape(x,e,p);
       temp=MtM*X+mu1*X*D +mu2*X;
       %temp=mtimesx(MtM,X,'BLAS')+mu1*mtimesx(X,D,'BLAS');
       y=temp(:);
 end
end
%%  Soft-thresholding for L21-norm minimization
function X=softThL21(B,lambda)
       [m,~]=size(B);
       D= spdiags(sqrt(sum(B.^2,2))-lambda/2,0,m,m);
       X=max(0,D)*normrGood(B);     

end
%% Simple soft-thresholding
function X=softTh(B,lambda)
       X=sign(B).*max(0,abs(B)-(lambda/2));
end
%%  Total variation
function opD=TVmatrix(m,n,str)

if str=='H' % This will give matrix for Horizontal Gradient
    D = spdiags([-ones(n,1) ones(n,1)],[0 1],n,n);
    D(n,:) = 0;
    D = kron(D,speye(m));
    
    
elseif str=='V' %This will give matrix for Verticle Gradient
   D = spdiags([-ones(m,1) ones(m,1)],[0 1],m,m);
   D(m,:) = 0;
   D = kron(speye(n),D);
   
end
opD=D;

end
%% PSNR calculation
function [cpsnr,psnr]=myPSNR(org,recon,skip)
%skip : to skip some boundary lines
org=org(skip+1:end-skip,skip+1:end-skip,:);
recon=recon(skip+1:end-skip,skip+1:end-skip,:);
  [m, n,~]=size(org);

if strcmp(class(org),class(recon))
    sse=squeeze(sum(sum((org-recon).^2))); %square sum of error  
    mse=sse./(m*n);  %mean square error of each band.
    maxval=squeeze(max(max(org)));
    psnr= 10*log10( (maxval.^2) ./mse);
    cpsnr=mean(psnr);
end
end

%% histogram equilization of each band
function img=myhisteq(img)
%make each band in zero to one range
for i=1:size(img,3)
    img(:,:,i)=mat2gray(img(:,:,i));
end
end
%%  add Gaussian noise of given SNR
function [noisy,sigma]=addGaussianNoise(img,snr)
%This function will add Gaussian noise of given SNR.
%img is image in size 


noisy=zeros(size(img));
for i=1:size(img,3)
    band=img(:,:,i);
    varNoise= norm(band(:))^2/(length(band(:)) * (10^ (snr/10)));
    noisy(:,:,i)=band+sqrt(varNoise)*randn(size(band));
end
    sigma=sqrt(varNoise);

end
%% this is built-in matlab function "normr" with small change to avoid infinity
function y = normrGood(x)
%NORMR Normalize rows of matrices.
%
%  <a href="matlab:doc normr">normr</a>(X) takes a single matrix or cell array of matrices and returns
%  the matrices with rows normalized to a length of one.
%
%  Here the rows of a random matrix are randomized.
%
%    x = <a href="matlab:doc rands">rands</a>(4,8);
%    y = <a href="matlab:doc normr">normr</a>(x)
%
%  See also NORMC.

% Mark Beale, 1-31-92
% Copyright 1992-2010 The MathWorks, Inc.
% $Revision: 1.1.8.4 $  $Date: 2012/08/21 01:03:41 $

% Checks
if nargin < 1,error(message('nnet:Args:NotEnough')); end
wasMatrix = ~iscell(x);
x = nntype.data('format',x,'Data');

% Compute
y = cell(size(x));
for i=1:numel(x)
  xi = x{i};
  cols = size(xi,2);
  n = 1 ./ sqrt(sum(xi.*xi,2));
  yi = xi .* n(:,ones(1,cols));
  yi(~isfinite(yi)) = 0;  %%%% this is the change
  y{i} = yi;
end

% Format
if wasMatrix, y = y{1}; end

end