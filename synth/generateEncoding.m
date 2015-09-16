function [A,W]=generateEncoding(N,kGrid,NS,EncMeth)

%GENERATENCODING   Generates spectral encoding
%   [A,W] = GENERATEENCODING(N,KGRID,NS,ENCMETH)
%   generates encoding patterns included in Fig. 4 of the paper
%   N specifies the dimensions of the image to be reconstructed
%   KGRID is a 1x3 cell array with the spectral coordinates along each axis
%   NS is a vector with each element specifying the number of shots (or of 
%   transform states)
%   ENCMETH is a cell array with each element specifying the encoding 
%   pattern (see Fig. 4 of the paper), one in 'LinSeq', 'RadSeq', 'LinPar',
%   'LinSqu', 'Random'
%   It returns:
%   A, cell array of LENGTH(NS) x LENGTH(ENCMETH) containing the sampling
%   scheme corresponding to the encoding method
%   W, k-space mask for converting to isotropic resolution
%

kk(:,:,1)=repmat(kGrid{1},[1 N(2)]);%k1 coordinates
kk(:,:,2)=repmat(kGrid{2},[N(1) 1]);%k2 coordinates
kkang=atan(kk(:,:,2)./kk(:,:,1));%Angular k-coordinates
kkrad=sqrt(sum(kk.^2,3));%Radial k-coordinates
[Y,indRsort]=sort(kkang(:));%Angular sorting of k-coordinates
[I,J]=ind2sub(N,indRsort);

W=single(kkrad<pi);%1 if |k|<pi, for isotropic resolution
W=ifftshift(ifftshift(W,1),2);%DC at first element

A=cell(length(NS),length(EncMeth));
for S=1:length(NS)
    for E=1:length(EncMeth)
        A{S}{E}=single(zeros([N(1:2) 1 1 NS(S)]));
        if strcmp(EncMeth{E},'LinSeq')%Cartesian sequential encoding
            for s=1:NS(S)
                Block=N(2)/NS(S);
                A{S}{E}(:,1+(s-1)*Block:s*Block,1,1,s)=1;
            end
        elseif strcmp(EncMeth{E},'RadSeq')%Radial sequential encoding
            for s=1:NS(S)
                Block=N(1)*N(2)/NS(S);
                IBlock=I(1+(s-1)*Block:s*Block);
                JBlock=J(1+(s-1)*Block:s*Block);
                for k=1:length(IBlock)
                    A{S}{E}(IBlock(k),JBlock(k),1,1,s)=1;
                end
            end
        elseif strcmp(EncMeth{E},'LinPar')%Cartesian parallel 1D encoding
            for s=1:NS(S)                  
                A{S}{E}(:,s:NS(S):end,1,1,s)=1;
            end
        elseif strcmp(EncMeth{E},'LinSqu')%Cartesian parallel 2D encoding                
            for m=1:N(1)
                for n=1:N(2)
                    A{S}{E}(m,n,1,1,mod(m+n-1,NS(S))+1)=1;
                end
            end                
        elseif strcmp(EncMeth{E},'Random')%Random encoding
            indRand=randperm(N(1)*N(2));
            [IR,JR]=ind2sub(N,indRand);               
            for s=1:NS(S)
                Block=N(1)*N(2)/NS(S);
                IRBlock=IR(1+(s-1)*Block:s*Block);
                JRBlock=JR(1+(s-1)*Block:s*Block);
                for k=1:length(IRBlock)
                    A{S}{E}(IRBlock(k),JRBlock(k),1,1,s)=1;
                end
            end                
        else
            error('Undefined Encoding method');
        end
        A{S}{E}=ifftshift(ifftshift(A{S}{E},1),2);%DC at first element
    end
end
