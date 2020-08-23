%---This is the PCA demension reduction part of group spatial fourier-ICA---
%---dimension reduction individual EEG data and temporal model catenation.
%---Then do group-level spatial fourier ICA.
% by Yongjie 5.22.2017

%%
function [Z,whiteningmatrix,dewhiteningmatrix]=f_reduceDim(data,pcadim,complexmixing)
%-Input
%---data:voxel positions x time sampals
%---pcadim:reduce dimension to number pcadim
%---complexmixing: complex value?
%-Output
%---Z: reduced data
%---whiteningmatrix
%---dewhiteningmatrix
%-------------%
disp('---Do individual PCA and writening')
zerotolerance=1e-7;% Value of the smallest eigenvalue relative to the maximum eigenvalue that is considered
[~,samplesN]=size(data);
X=data';
clear data;
% Substract mean value from voxels
Xmean=mean(X);
X=X-Xmean(ones(samplesN,1),:);
clear Xmean;
% PCA
Xcov=cov(X);
if ~complexmixing, Xcov=real(Xcov); else Xcov=conj(Xcov); end
[Ex,Dx]=eig(Xcov);
[lambda,order]=sort(diag(Dx),'descend');
clear Xcov Dx;
% checks if negative eigenvalues
if any(lambda(1:pcadim)<0), warning('Negative eigenvalues! Reducting PCA dimension...'),end
%Chek for eigenvalues near zero(relative to the maximum eigenvalue)
zeroEig=sum((lambda(1:pcadim)/lambda(1))<zerotolerance); % 
% Adjust dimensions if necessary 
pcadim=pcadim-zeroEig;

% Construct whitening and dewritening matrices
lamsqrt=sqrt(lambda(1:pcadim));
lamsqrtinv=1./lamsqrt;
Ex=Ex(:,order(1:pcadim));
whiteningmatrix=diag(lamsqrtinv)*Ex';
dewhiteningmatrix=Ex*diag(lamsqrt);
clear lambda order zeroEig lamsqrt lamsqrtinv Ex;

% Reduce dimensions and whiten data.
Z=whiteningmatrix*transpose(X);
disp('---PCA and wrtening Done.')
end
