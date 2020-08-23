% """
% Created 5.27.2017
% ===============================================================
% Code for demonstrating group analysis of ongoing EEG with spatial Fourier-ICA
% ===============================================================
% @author: Yongjie Zhu
% """
%%
%---1.origdata(channel*time),spectrom(origdata),got dataChFrpF
%---2.mapping dataChFrpF to cortex space,dataTCrtFrq.
%---3.dimension redunction individual data
%---4.temporal catenation
%---5.group spatial fourier-ICA
%---6.back projction to individual space

%% dimension reduction

%-----using STFT transfor data(channel*time) to channel*frequency*time---%
% load('E:\EEG_data_HP_4Hz.mat')
% load G
[Nvoxel,Nch]=size(G);
Nsubs=14;
pcadim=25;
fs=256;
nfft=512;maxfreq=35;maxid=ceil(maxfreq*nfft/fs);
groupData=zeros(Nsubs*pcadim,Nvoxel*maxid);
Nfreq=maxid;
for subi=1:Nsubs
    disp(['PCA reduction:' num2str(subi)])
    X=double(data(:,:,subi)');
    dataChFrqT=[];
    for i=1:Nch
        [S,F,T]=spectrogram(X(i,:),3*fs,2*fs,nfft,fs);
        dataChFrqT(i,:,:)=S(1:maxid,:);
    end
    clear X;
    %----liner inverse transform.from sensor space(channel*Frequency*Time)
    %----to cortical space(cortex*Frequency*time)
    %----G:inverse oporation.obtained from brainstorm-----%
    [~,Nt]=size(S); % [NbCortex,NbFrq,NbT]=size(dataCrtFrqT);
    G=single(G);dataChFrqT=single(dataChFrqT); % save memmory
    %dataCrtFrqT=zeros(NbCortex,NbFrq,NbT);
    for i=1:Nt
        dataCrtFrqT(:,:,i)=G*dataChFrqT(:,:,i);
    end
    clear dataChFrqT;
    %---reshape matrix from cortex*Frequency*time to time*cortexFrequency---%
    for i=1:Nt
        dataTCrtF(i,:)=reshape(dataCrtFrqT(:,:,i),[1,Nvoxel*Nfreq]);
    end
    clear dataCrtFrqT;
    % dimension reduction
    [Z,whiteningmatrix,dewhiteningmatrix]=f_reduceDim(dataTCrtF,pcadim,true);
    groupData((subi-1)*pcadim+1:subi*pcadim,:)=Z; clear Z ;
    %filePath='C:\Users\yozhu\matlabdata\groupSpatialFourierICA\';
    fileName=[num2str(subi) '_' 'writenMatrix.mat'];
    save(fileName,'whiteningmatrix','dewhiteningmatrix');
end
clear dataTCrtF;
%%

Ncomps=25;
options = struct('samplfreq',256,'pcadim',Nsubs*pcadim,'minfreq',1,'maxfreq',30,...
    'complexmixing',true,'components',Ncomps);
[S,A,W]=cfastICA(groupData,options);

for i=1:Ncomps
    specspa=reshape(S(i,:),[Nvoxel,Nfreq]);
    spectra(i,:)=mean(abs(specspa),1);
    %spartial(i,:)=mean(abs(specspa),2)';
    maxV=max(spectra(i,:));
    topId=find(spectra(i,:)>0.95*maxV);
    spt(i,:)=mean(abs(specspa(:,topId)),2);
end
save('AWS.mat','A','W','S');
spectrum=spectra;spt=spt';
%% back project to recover time courses
groupTime=reshape(A,pcadim,Nsubs,Ncomps);
groupTime=permute(groupTime,[1 3 2]);
for i=1:Nsubs
    fileName=[num2str(i) '_' 'writenMatrix.mat'];
    load(fileName)
    bpgroupTime(:,:,i)=dewhiteningmatrix*groupTime(:,:,i);
end
Timecourses=abs(bpgroupTime);
%%