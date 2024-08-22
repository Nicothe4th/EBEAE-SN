
function [Y, A, Po, Pr]=spxls_gmlm(Z,z,n,nsp,maxiter,Po);
% Z-> Z hipercubo
% z-> 1 slide del hipercubo para obtener los superpixeles
% n-> #endmembers
%nsp-> # numero de superpixeles
%
%% superpixeles 
[eti,numeti]=superpixels(z,nsp);
[L,N]=size(Z);
%%
rho=1e-5;%0.005 base 1
l1=0.001;%.001 base 1%sparsity
l2=0.001; %smoothness of nonlinearity matrix
l3=0.004; %0.004 base 1smoothness of nonlinearity matrix
%Po=vca(Z,n);
etir=reshape(eti,N,1)';
Y=zeros(L,N);
At=zeros(n,N);
Prt=zeros(1,N);
%%
 for h=1:numeti
     n1=find(etir==h);
     [Y(:,n1), A(:,n1), P, Pr(n1)]= g_mlm(Z(:,n1),n,rho,l1,l2,l3,maxiter,Po);
 end

end

%% ADMM for the proposed graph regularided(G-MLM)
function [Y, At, Po, Prt]= g_mlm(X,n,rho,l1,l2,l3,maxiter,Po)
%%
% E endmember matriz-> VCA 
% Ao abundance matriz-> FCLS 
%%
[L,N]=size(X);
if nargin == 7
    Po=VCA(X./sum(X),n);
end
At=FCLS(X,Po);
elim=1e-4;
%% condiciones iniciales
Yh=(abs(Po)*At)./sum(abs(Po*At),1);
dmin=400*(norm(X-Yh,'Fro').^2)/(N*L);

W=zeros(N,N);%matriz de afinidad
for i=1:N
    for j=i:N
        value=norm(X(:,i)-X(:,j))^2;
        if (value<dmin)
            W(i,j)=1;
            W(j,i)=1;
        end
    end
end
 Prt =(ones(1,L)*((Yh-Yh.*X).*(Yh-X)))./((ones(1,L)*((Yh-Yh.*X).^2))+rho*ones(1,N));%ec26 H y M =0
 Gt=At;
 Go=Gt;
 Ht=Prt;
 Ho=Ht;
 M1=zeros(size(At));
 M2=zeros(size(Prt));
 %% 
 for i=1:maxiter
 At=Ec17(At,Po,X,rho,Prt,l3,W,Gt,M1);%update A Ec(17)
 Gt=Ec24(At,M1,l1,rho);%update  G Ec(24)
 Prt=Ec26(Po,At,X,rho,Ht,M2);%update Prt EC(26)
 Ht=Ec29(Prt,rho,M2,l3,W);%update Ht Ec(29)
 M1=M1+At-Gt;%update M1 Ec(30)
 M2=M2+Prt-Ht;%update M2 Ec(31)
 resp=norm([At-Gt;Prt-Ht],'fro');
 resd=norm([Go-Gt; Ho-Ht],'fro');
 %disp(['resp = ' num2str(resp)]);
 %disp(['resd = ' num2str(resd)]);
 Go=Gt;
 Ho=Ht;
 if(resp<=elim && resd<=elim)
     %disp('break');
     break
 end
 end
 ylm=(abs(Po)*At)./sum(abs(Po*At),1);
 for i=1:N
    Y(:,i)=((1-Prt(i)).*ylm(:,i))./(1-(Prt(i).*ylm(:,i)));
    end
end
 function Ht=Ec29(Prt,rho,M2,l3,W)
 [N,N]=size(W);
 sumw=sum(W);
 L=diag(sumw)-W;
 H=rho*(Prt+M2)*pinv(l3*L+rho*ones(N));
 Ht=min(1,H);
 end
 function Prt=Ec26(Po,At,X,rho,Ht,M2)
 [L,N]=size(X);
 [n,m]=size(At);
 Y=(Po*At)./sum(Po*At,1);
 Prt=((ones(1,L)*((Y-Y.*X).*(Y-X)))+rho*(Ht-M2))./((ones(1,L)*((Y-Y.*X).^2))+rho*ones(1,N));
 end
 function Gt1=Ec24(At,M1,l1,rho)
 b=l1/rho;
 Z=At+M1;
 [C1r,C1c]=find(Z>b);
 [C3r,C3c]=find(abs(Z<-b));
 S=zeros(size(At));
 if(prod(size(C1r))~=0)
    S(C1r,C1c)=Z(C1r,C1c)-b;
 end
 if(prod(size(C3r))~=0)
    S(C3r,C3c)=Z(C3r,C3c)+b;
 end
 Gt1=max(0,S);
end
 function At1=Ec17(At,Po,X,rho,Prt,l3,W,Gt,M1)
  [L,N]=size(X);
  [n,m]=size(At);
  for j=1:N
    Pj_hat=Po.*((1-Prt(j))*ones(L,n)+Prt(j)*X(:,j)*ones(1,n));
    B=double((Pj_hat'*Pj_hat)+(rho+l3*sum(W(:,j)))*eye(n));
    C=pinv(B)*ones(n,1)*inv(ones(1,n)*pinv(B)*ones(n,1));
    w=Pj_hat'*X(:,j)+(rho*(Gt(:,j)-M1(:,j)))+(l3*(sum(W(:,j).*At(:,j)'))');
    a=pinv(B)*w-C*(ones(1,n)*inv(B)*w-1);
    At1(:,j)=a;
    
  end
 end
%% vca
 function [M,indices,snrE]=VCA(R,p)
[L,N]=size(R);
SNRth=15+10*log10(p);
rm=mean(R,2);
rzm=R-repmat(rm,1,N);
[U,S,V]=svds(rzm*rzm'./N,p);
rd=U'*rzm;
pr=sum(R(:).^2)/N;
prp=sum(rd(:).^2)/N +rm'*rm;
snrE=abs(10*log10((prp-(p/L)*pr)/(pr-prp)));

if (snrE>SNRth)
    d=p;
    [Ud,S,V]=svds(R*R'./N,d);
    X=Ud'*R;
    u=mean(X,2);
    Y=X ./ repmat( sum( X .* repmat(u,[1 N]) ) ,[d 1]);
else
    d=p-1;
      r_line=mean(R')';
    Ud=pca((R-r_line)*(R-r_line)'./N,'NumComponents',d);
    R_r=R-repmat(r_line,1,N);
    X=Ud'*R_r;
    c=zeros(N,1);
    for j=1:N
        c(j)=norm(X(:,j));
    end
    c=repmat(max(c),1,N);
    Y=[X;c];
end
eu=zeros(p,1);
eu(p)=1;
A=zeros(p,p);
A(:,1)=eu;
I=eye(p);
k=zeros(N,1);
for i=1:p
    w=rand(p,1);
    tmpnum=(I-A*pinv(A))*w;
    f=tmpnum/norm(tmpnum);
    v=f'*Y;
    k=abs(v);
    [unused,k]=max(k);
    A(:,i)=Y(:,k);
    indices(i)=k;
end
if (snrE>SNRth)
    M=Ud*X(:,indices);
else
    M=Ud*X(:,indices)+repmat(r_line,1,p);
end
 end
function A=FCLS(X,P0)
    [L,N]=size(X);
    p=size(P0,2);
    A=zeros(p,N);
    delta = 1/(10*max(max(P0)));
    for i=1:N
        s = [delta.*X(:,i);1];
        M = [delta.*P0; ones(1,p)];
        A(:,i)= NCLS(s, M, -1e-6);
    end
    A=A./sum(A,1);
end



function [abundance]=NCLS(x, MatrixZ, tol)
    % input MatrixZ is the signatures of endmembers. It is of size [bands p].
    % input x is the signature whose abundance is to be estimated.
    % output abundance is the abundance of each material in r1. It is of size [p 1]. % This function is written according to Dr. Chang?s first book , P 47
    M=size(MatrixZ,2);
    R=zeros(M,1);
    P=ones(M,1);
    invMtM=(MatrixZ'*MatrixZ)^(-1);
    Alpha_ls=invMtM*MatrixZ'*x;
    Alpha_ncls=Alpha_ls;
    min_Alpha_ncls=min(Alpha_ncls);
    j=0;
    while(min_Alpha_ncls<-tol && j<500) 
        j = j+1;
        for II=1:M
            if((Alpha_ncls(II)<0)&&(P(II)==1)) R(II)=1;
                P(II)=0;
            end %%% end of if (Alpha_ncls(II)<0) 
        end % end of for II=1:M
        S = R;
    goto_step6=1; 
    counter = 0; 
    while(goto_step6==1)
        index_for_Lamda = find(R==1);
        Alpha_R = Alpha_ls(index_for_Lamda);
        Sai = invMtM(index_for_Lamda,index_for_Lamda);
        inv_Sai = (Sai)^(-1); % remember inversion of Sai 
        Lamda=inv_Sai*Alpha_R;
        [max_Lamda,index_Max_Lamda]=max(Lamda); 
        counter = counter+1;
        %disp(['max_Lamda= ',num2str(le(max_Lamda,0))])
        %disp(['counter= ',num2str(eq(counter, 200))])
        if (isempty(max_Lamda))
            break;
        end
        if ( le(max_Lamda,0) || eq(counter, 200) )
            break; 
        end
        temp_i = inv_Sai; % simplify the inversion of matrix 
        temp_i(1,:) = inv_Sai(index_Max_Lamda,:);
        if  (index_Max_Lamda>1)
            temp_i(2:index_Max_Lamda,:) = inv_Sai(1:index_Max_Lamda-1,:); 
        end
        inv_Sai_ex = temp_i;
        inv_Sai_ex(:,1) = temp_i(:,index_Max_Lamda);
        if index_Max_Lamda>1
            inv_Sai_ex(:,2:index_Max_Lamda) = temp_i(:,1:index_Max_Lamda-1); 
        end
        inv_Sai_next = inv_Sai_ex(2:end,2:end) - inv_Sai_ex(2:end,1)*inv_Sai_ex (1,2:end)/inv_Sai_ex(1,1);
        P(index_for_Lamda(index_Max_Lamda))=1; 
        R(index_for_Lamda(index_Max_Lamda))=0; 
        index_for_Lamda(index_Max_Lamda) = [];
        Alpha_R = Alpha_ls(index_for_Lamda); 
        Lamda=inv_Sai_next*Alpha_R;
        Phai_column = invMtM(:,index_for_Lamda);
        if (size(Phai_column,2)~=0) 
         Alpha_s=Alpha_ls-Phai_column*Lamda;
        else
            Alpha_s=Alpha_ls;
        end
        goto_step6=0;
        for II=1:M
            if ((S(II)==1)&&(Alpha_s(II)<0))
                P(II)=0;
                R(II)=1;
                goto_step6=1;
            end
        end
    end % end of while (gotostep6==1)
    index_for_Phai = find(R==1);
    Phai_column = invMtM(:,index_for_Phai);
    if (size(Phai_column,2)~=0) 
        Alpha_ncls=Alpha_ls-Phai_column*Lamda;
    else
            Alpha_ncls=Alpha_ls;
    end
    min_Alpha_ncls=min(Alpha_ncls); 
end % end of while
abundance=zeros(M,1); 
for II=1:M
    if (Alpha_ncls(II)>0) 
        abundance(II)=Alpha_ncls(II);
    end
end
return;
end