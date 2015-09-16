function [et,etg,eth]=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,T,di,cg)

%PRECOMPUTATIONSSINCRIGIDTRANSFORM   Precomputes diagonal operators for 3D
%sinc-based rigid transform
%   [ET ETG ETH]=PRECOMPUTATIONSSINCRIGIDTRANSFORM(KGRID,KKGRID,RKGRID,T,DI,CG)
%   precomputes the k-space phase multiplicative factors used for 
%   transforming the images with sinc-based interpolation.
%   KGRID is a grid of points in the spectral domain as provided by the
%   function generateGrids
%   KKGRID is a grid of points in the spectral-spectral domain as provided 
%   by the function generateGrids
%   RKGRID is a grid of points in the spatial-spectral domain as provided
%   by the function generateGrids
%   T are the parameters of the transform
%   DI is a flag to indicate whether to perform direct or inverse transform
%   CG is a flag to indicate whether to calculate the derivative terms
%   It returns (see the paper for further details):
%   ET, the transform factors
%   ETG, the transform gradient factors
%   ETH, the transform Hessian factors
%

theta=T(:,:,:,:,:,4:6);%Rotation parameters
t=((-1)^di)*1i*T(:,:,:,:,:,1:3);%Translation parameters

tantheta2=theta/2;
tantheta2=tan(tantheta2);
tantheta2j=((-1)^(di-1))*1i*tantheta2;
sintheta=sin(theta);
sintheta=((-1)^di)*1i*sintheta;    
clear T 
if cg>0
    tantheta=tan(theta);
    tanthetacuad=tantheta2.*tantheta2;
    tanthetacuad=(1+tanthetacuad)/2;
    costheta=cos(theta);    
end       
clear theta

per=[1 3 2;2 1 3];
for m=1:3
    et{2}{m}=exp(bsxfun(@times,tantheta2j(:,:,:,:,:,m),rkGrid{1}{m}));%Tan exponential
    et{3}{m}=exp(bsxfun(@times,sintheta(:,:,:,:,:,m),rkGrid{2}{m}));%Sin exponential
end

clear tantheta2j sintheta

if cg>0%Gradient terms
    for m=1:3     
        etg{2}{m}=bsxfun(@times,tanthetacuad(:,:,:,:,:,m),1i*rkGrid{1}{m});%Tan derivative
        etg{3}{m}=bsxfun(@times,costheta(:,:,:,:,:,m),-1i*rkGrid{2}{m});%Sin derivative
        if cg==2%Hessian terms
            eth{2}{m}=bsxfun(@plus,tantheta2(:,:,:,:,:,m),etg{2}{m});
            eth{3}{m}=bsxfun(@plus,-tantheta(:,:,:,:,:,m),etg{3}{m});        
        end
    end
    clear tanthetacuad costheta tantheta2 tantheta

    for m=2:3        
        for n=1:3        
            etg{m}{n}=etg{m}{n}.*et{m}{n};
            if cg==2
                eth{m}{n}=eth{m}{n}.*etg{m}{n};
            end
            etg{m}{n}=ifftshift(etg{m}{n},per(m-1,n));
            if cg==2
                eth{m}{n}=ifftshift(eth{m}{n},per(m-1,n));
            end
        end
    end
end
for m=1:3
    for n=2:3
        et{n}{m}=ifftshift(et{n}{m},per(n-1,m));
    end
end

et{1}=bsxfun(@plus,bsxfun(@plus,bsxfun(@times,t(:,:,:,:,:,1),kGrid{1}),bsxfun(@times,t(:,:,:,:,:,2),kGrid{2})),bsxfun(@times,t(:,:,:,:,:,3),kGrid{3}));
et{1}=exp(et{1});
clear t
if cg>0   
    for m=1:3
        etg{1}{m}=bsxfun(@times,-1i*kGrid{m},et{1});        
        for n=1:3
            etg{1}{m}=ifftshift(etg{1}{m},n);
        end
    end    
    if cg==2
        for m=1:6
            eth{1}{m}=bsxfun(@times,-kkGrid{m},et{1});
            for n=1:3
                eth{1}{m}=ifftshift(eth{1}{m},n);        
            end
        end
    end
end

for m=1:3
    et{1}=ifftshift(et{1},m);    
end
