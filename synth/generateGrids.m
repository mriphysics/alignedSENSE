function [kGrid,kkGrid,rGrid,rkGrid]=generateGrids(N)

%GENERATEGRIDS   Generates spectral and spatial coordinate grids
%   [KGRID,KKGRID,RGRID,RKGRID] = GENERATEGRIDS(N)
%   generates coordinate grids using the image size (no spectral
%   subsampling is assumed so that spatial and spectral grids match in size)
%   N specifies the dimensions of the image to be reconstructed
%   It returns:
%   KGRID, a 1x3 cell array with the spectral coordinates along each axis
%   KKGRID, a 1x6 cell array with the spectral-spectral coordinates along 
%   each possible axes pair (indexed by [1 2 3 1 1 2;1 2 3 2 3 3])
%   RGRID, a 1x3 cell array with the spatial coordinates along each axis
%   RKGRID, a 2x3 cell array where the first dimension of the cell indexes
%   the shearing orientation and the second indexes the spatio-spectral
%   axes pair used for that particular shearing (given by [1 3 2;2 1 3]).
%   See ec. 4 of the paper for further details
%

%Spectral grid
kGrid=cell(1,3);
for m=1:3
    kGrid{m}=-floor(N(m)/2):ceil(N(m)/2)-1;
end
kGrid{1}=single(permute(kGrid{1},[2 1]));
kGrid{3}=single(permute(kGrid{3},[3 2 1]));
for m=1:3
    kGrid{m}=2*pi*kGrid{m}/N(m);
end
fact=[1 2 3 1 1 2;1 2 3 2 3 3];
kkGrid=cell(1,6);
for m=1:6
    kkGrid{m}=bsxfun(@times,kGrid{fact(1,m)},kGrid{fact(2,m)});
end

%Spatial grid
cent=ceil(N/2);%Center of rotation
per=[1 3 2;2 1 3];
rGrid=cell(1,3);
for m=1:3
    rGrid{m}=(1:N(m))-cent(m);
end    
rGrid{1}=single(permute(rGrid{1},[2 1]));
rGrid{3}=single(permute(rGrid{3},[3 2 1]));
rkGrid=cell(2,3);
for n=1:2
    for m=1:3
        rkGrid{n}{m}=bsxfun(@times,rGrid{per(3-n,m)},kGrid{per(n,m)});%First index denotes that the k takes the first dimension of xkgrid; second index denotes the rotation        
    end
end