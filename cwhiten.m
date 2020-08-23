function [whiteningMatrix, dewhiteningMatrix] = cwhiten(origdata,options)
%WHITENV - Whitenv vectors.
%
% [newVectors, whiteningMatrix, dewhiteningMatrix] = ...
%                               whitenv(vectors, E, D, verbose);
%
% Whitens the data (row vectors) and reduces dimension. Returns
% the whitened vectors (row vectors), whitening and dewhitening matrices.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% E             Eigenvector matrix from function 'pcamat'
% D             Diagonal eigenvalue matrix from function 'pcamat'
% verbose       Optional. Default is 'on'
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also PCAMAT

% @(#)$Id: whitenv.m,v 1.3 2003/10/12 09:04:43 jarmo Exp $

% ========================================================
%% cfastICA modifyed by yongjie 
disp('Spatially whitening data with PCA dimension reduction');
pcadim=options.pcadim;
components=options.components;
complexmixing=options.complexmixing;
[~,N]=size(origdata);
Xmat_c=origdata';
clear origdata;
% Substract mean value from channels
Xmat_c_mv=mean(Xmat_c);
Xmat_c=Xmat_c-Xmat_c_mv(ones(N,1),:);
clear Xmat_c_mv

%Do PCA (on matrix Xmat_c)
covmat=cov(Xmat_c);
if ~complexmixing, covmat=real(covmat); else covmat=conj(covmat); end
[Ec, Dc] = eig(covmat);
[d,order] = sort(diag(Dc),'descend');
clear covmat Dc

% Checks for negative eigenvalues
if any(d(1:pcadim)<0), warning('Negative eigenvalues! Reducing PCA and ICA dimension...'), end
% Check for eigenvalues near zero (relative to the maximum eigenvalue)
zeroeigval=sum((d(1:pcadim)/d(1))<1e-7);
% Adjust dimensions if necessary (because zero eigenvalues were found)
pcadim=pcadim-zeroeigval;
if pcadim<components, components=pcadim; end
if zeroeigval, fprintf('PCA dimension is %d and ICA dimension is %d\n',pcadim,components); end


% Construct whitening and dewhitening matrices
dsqrt = sqrt(d(1:pcadim));
dsqrtinv = 1./dsqrt;
Ec = Ec(:,order(1:pcadim));
whiteningMatrix=diag(dsqrtinv)*Ec';
dewhiteningMatrix=Ec*diag(dsqrt);
clear d order zeroeigval dsqrt dsqrtinv Ec


