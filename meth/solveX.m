function x=solveX(x,y,M,T,S,SH,Precond,A,kGrid,rkGrid,nX,gpu,BlockSize)

%SOLVEX   Reconstructs an image using CG SENSE with multishot alignment
%   X=SOLVEX(X,Y,M,T,S,SCONJ,PRECOND,A,KGRID,RKGRID,NX,GPU,BLOCKSIZE) 
%   computes the best X for a given T
%   X is the image to be reconstructed
%   Y is the measured data
%   M is a spatial mask to constrain the solution
%   T are the transform parameters
%   S is the coil-array sensitivity map
%   SH is the conjugate transpose of the coil-array sensitivity map
%   PRECOND is the preconditioner
%   A is a sampling mask
%   KGRID is a grid of points in the spectral domain
%   RKGRID is a grid of points in the spatial-spectral domain
%   NX is the number of iterations of the CG algorithm
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   BLOCKSIZE indicates the number of shots to be processed in a chunk
%   It returns:
%   X, the reconstructed image
%

NX=size(x);NX(end+1:3)=1;
NT=size(T);
NY=size(y);
NS=size(S);

NRun=ceil(NT(5)/BlockSize);
NRunRem=mod(NT(5),BlockSize);
vS=cell(NRun,1);
for s=1:NRun
    if s~=NRun || NRunRem==0
        vS{s}=(s-1)*BlockSize+1:s*BlockSize;
    else
        vS{s}=(s-1)*BlockSize+1:(s-1)*BlockSize+NRunRem;
    end
    etDir{s}=precomputationsSincRigidTransform(kGrid,[],rkGrid,T(:,:,:,:,vS{s},:),1,0);
    etInv{s}=precomputationsSincRigidTransform(kGrid,[],rkGrid,T(:,:,:,:,vS{s},:),0,0);      
end

%SENSE undersampling parameters; just works for SENSE factors lower than 3
FOV=single(ones(2,2));
iFOV=single(ones(2,2));
for n=1:2
    if NS(n)~=NY(n)
        disc=3*NY(n)-NS(n);
        iFOV(n,1)=floor(disc/2);
        iFOV(n,2)=ceil(disc/2);
        over=NS(n)-NY(n);
        FOV(n,1)=floor(over/2);
        FOV(n,2)=ceil(over/2);
    end
end

if gpu
    x=gpuArray(x);A=gpuArray(A);y=gpuArray(y);S=gpuArray(S);SH=gpuArray(SH);M=gpuArray(M);Precond=gpuArray(Precond);
end

%To balance high-frequency armonics in the reconstruction for isotropic
%resolution
ATot=1-sum(A,5);
ATot=ATot/NT(5);
A=bsxfun(@plus,A,ATot);

clear ATot

if gpu
    yEnd=gpuArray(single(zeros(NX))); 
else
    yEnd=single(zeros(NX)); 
end

for s=1:NRun    
    if gpu
        etInvS{1}=gpuArray(etInv{s}{1});
        for m=1:3
            for n=2:3
                etInvS{n}{m}=gpuArray(etInv{s}{n}{m});                
            end
        end
    else
        etInvS{1}=etInv{s}{1};
        for m=1:3
            for n=2:3
                etInvS{n}{m}=etInv{s}{n}{m};                
            end
        end
    end
    %A'y
    yS=bsxfun(@times,y,A(:,:,:,:,vS{s}));
    %F'A'y
    for m=1:2
        yS=ifftGPU(yS,m,gpu);
    end    
    for m=1:2
        yS=isense(yS,m,NS(m),NY(m),iFOV(m,:));
    end 
    %S'F'A'y
    yS=sum(bsxfun(@times,yS,SH),4);  
    %T'S'F'A'y
    yEnd=yEnd+sincRigidTransform(yS,etInvS,0,gpu);        
end

y=yEnd;
y=M.*y;

Ap=applyCG(x);
r=y-Ap;
clear y yEnd yS

z=Precond.*r;
p=z; 
rsold=sum(sum(sum(conj(z).*r)));

%Iterations
for n=1:nX
    Ap=applyCG(p);
    al=conj(rsold)/sum(sum(sum(conj(p).*Ap)));
    
    x=x+al*p;
    r=r-al*Ap;
    z=Precond.*r;
    rsnew=sum(sum(sum(conj(z).*r)));
    be=rsnew/rsold;
    p=z+be*p;
    rsold=rsnew;
    if sqrt(abs(rsnew))<1e-10
        break;
    end
end

if gpu
    x=gather(x);
end

%%
function x=applyCG(x)
    if gpu
        xEnd=gpuArray(single(zeros(NX(1:3))));
    else
        xEnd=single(zeros(NX(1:3)));
    end   
    for s=1:NRun
        if gpu    
            etDirS{1}=gpuArray(etDir{s}{1});
            etInvS{1}=gpuArray(etInv{s}{1});
            for m=1:3
                for l=2:3
                    etDirS{l}{m}=gpuArray(etDir{s}{l}{m});
                    etInvS{l}{m}=gpuArray(etInv{s}{l}{m});                
                end
            end
        else
            etDirS{1}=etDir{s}{1};
            etInvS{1}=etInv{s}{1};
            for m=1:3
                for l=2:3
                    etDirS{l}{m}=etDir{s}{l}{m};
                    etInvS{l}{m}=etInv{s}{l}{m};     
                end
            end
        end
        %Tx
        xS=sincRigidTransform(x,etDirS,1,gpu);
        %STx
        xS=bsxfun(@times,xS,S);
        %FSTx
        for m=1:2
            xS=sense(xS,m,NS(m),NY(m),FOV(m,:));
        end        
        for m=1:2
            xS=fftGPU(xS,m,gpu);
        end
        %A'AFSTx
        xS=bsxfun(@times,xS,A(:,:,:,:,vS{s}));
        %F'A'AFSTx
        for m=1:2
            xS=ifftGPU(xS,m,gpu);
        end               
        for m=1:2
            xS=isense(xS,m,NS(m),NY(m),iFOV(m,:));
        end       
        %S'F'A'AFSTx
        xS=sum(bsxfun(@times,xS,SH),4);
        %T'S'F'A'AFSTx
        xEnd=xEnd+sincRigidTransform(xS,etInvS,0,gpu);
    end
    x=xEnd;
    x=x.*M;
end

end
