
% """
% Created 2.10.2017
% ===============================================================
% Code for demonstration of analyzing ongoing EEG with spatial Fourier-ICA
% ===============================================================
% @author: Yongjie Zhu
% """
%%
%-----using STFT transfor data(channel*time) into channel*frequency*time---%
load('E:\EEG_data_HP_4Hz.mat')
load G
data=double(data(:,:,1)');
fs=256;
nfft=512;maxfreq=35;maxid=ceil(maxfreq*nfft/fs);
Nch=64;
dataChFrqT=[];
for i=1:Nch
    [S,F,T]=spectrogram(data(i,:),3*fs,2*fs,nfft,fs);
    dataChFrqT(i,:,:)=S(1:maxid,:);
end
%----liner inverse transform.from sensor space(channel*Frequency*Time) 
%----to cortical space(cortex*Frequency*time)
%----G:inverse oporation.obtained from brainstorm-----%
[~,Nt]=size(S); % 
Nfrq=maxid;
Nvoxel=size(G,1);
G=single(G);dataChFrqT=single(dataChFrqT); % save memmory

for i=1:Nt
    dataCrtFrqT(:,:,i)=G*dataChFrqT(:,:,i);
end
clear dataChFrqT G;
%---reshape matrix from cortex*Frequency*time to time*cortexFrequency---%
for i=1:Nt
    dataTVF(i,:)=reshape(dataCrtFrqT(:,:,i),[1,Nvoxel*Nfrq]);
end
% dataTCrtF=reshape(dataCrtFrqT,[NbT,NbCortex*NbFrq]);
clear dataCrtFrqT;
%---complex-valued fastICA---------%
%---Mix:time*Nbcomp;ICs:Nbcomp*cortex/Frequency---%
Ncomps=25;
options = struct('samplfreq',256,'pcadim',160,'minfreq',1,'maxfreq',30,...
    'complexmixing',true,'components',Ncomps);
[S,A,W]=cfastICA(dataTVF,options);

%---Extracting spectra, spatial map and timecourses---%

timecourse=abs(A);

for i=1:Ncomps
    specspa=reshape(S(i,:),[NbCortex,NbFrq]);
    spectra(i,:)=mean(abs(specspa),1);
    %spartial(i,:)=mean(abs(specspa),2)';
    maxV=max(spectra(i,:));
    topId=find(spectra(i,:)>0.95*maxV);
    spt(i,:)=mean(abs(specspa(:,topId)),2)';
end
