N=4;                % Number of End-members
Nsamples=120;       % Size of the Squared Image Nsamples x Nsamples 
SNR=30;             % Additive Gaussian noise
density=0.01;       % Density of sparse noise component
ModelType=0;        % 0 --> Linear Mixing Model and 5 --> Multilinear Model

[Z,P0,A0,V0]=SphericGaussian_Sparse_Synth(SNR,density,ModelType);
[L,K]=size(Z);
Zz=Z./repmat(sum(Z,1),[L,1]);
nRow=Nsamples;
nCol=Nsamples;
load('DatasetSynth/EndMembersHbHbO2FatWater.mat');

h1=figure;
h2=figure;
figure(h1);
subplot(1,2,1)
plot(wavelength,P0,'linewidth',2); grid on;
legend('Hb','HbO_2','Fat','Water','FontSize',12);
xlabel('wavelength (nm)'); axis([min(wavelength) max(wavelength) min(P0(:)) max(P0(:))]);
ylabel('Normalized absorbance','FontSize',12)
title('(a) HSI','FontSize',12)

figure(h2);
for i=1:4
    subplot(2,4,i);
    imagesc(reshape(A0(i,:),nRow,nCol),[0,1]);
    xlabel('pixels'); ylabel('pixels'); title(['End-member ' num2str(i)]);
end
subplot(2,4,2); title({'(a) HSI', 'End-member 2'});

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[Z,P0,A0,V0]=mFLIM_Sparse_Synth(N,Nsamples,0.25e-9,SNR,density);
[L,K]=size(Z);
Zz=Z./repmat(sum(Z,1),[L,1]);
nRow=Nsamples;
nCol=Nsamples;

figure(h1);
subplot(1,2,2)
plot((0.25)*(1:size(P0,1)),P0,'linewidth',2); grid on;
legend('End-member 1','End-member 2','End-member 3','End-member 4');
xlabel('time (ns)')
ylabel('Normalized intensity')
axis([0 0.25*size(P0,1) 0 max(P0(:))]);
title('(b) m-FLIM')
set(h1,'PaperPositionMode','auto')


figure(h2);
for i=1:4
    subplot(2,4,i+4);
    imagesc(reshape(A0(i,:),nRow,nCol),[0,1]);
    xlabel('pixels'); ylabel('pixels'); title(['End-member ' num2str(i)]);
end
subplot(2,4,6); title({'(b) m-FLIM', 'End-member 2'});
set(h2,'PaperPositionMode','auto')