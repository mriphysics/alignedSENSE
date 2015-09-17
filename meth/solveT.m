function [T,x,w,flagw]=solveT(x,y,M,T,S,A,kGrid,kkGrid,rGrid,rkGrid,nT,w,flagw,gpu,BlockSize,meanT)

%SOLVET   Estimates the rigid transforms each shot has been subject to
%using Newton's method
%   [T X W FLAGW]=SOLVET(X,Y,T,S,A,KGRID,KKGRID,RKGRID,NT,W,FLAGW,GPU,BLOCKSIZE)
%   computes the best T for a given X
%   X is the reconstructed data
%   Y is the measured data
%   M is a spatial mask to constrain the solution
%   T are the transform parameters
%   S is the coil-array sensitivity map
%   A is a sampling mask
%   KGRID is a grid of points in the spectral domain
%   KKGRID is a grid of points in the spectral-spectral domain
%   RKGRID is a grid of points in the spatial-spectral domain
%   NT is the number of iterations of the Newton's method (not carefully
%   checked for NT>1)
%   W is the w in ec. (27) of the paper
%   FLAGW keeps track of previous change in the energy of the solution to 
%   update w
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   BLOCKSIZE indicates the number of shots to be processed in a chunk
%   MEANT is a flag that determines whether to use ec. (28) of the paper,
%   defaults to 1
%   It returns:
%   T, the estimated set of rigid transforms
%   X, the reconstructed image after applying correction in ec. (28) of the
%   paper
%   W, the updated w in ec. (27) of the paper
%   FLAGW, the updated track of previous change in the energy of the
%   solution
%   

if ~exist('meanT','var');meanT=1;end

multA=1.2;%Factor to divide the weight (that regularizes the Hessian matrix) when H(end)<H(end-1)
multB=2;%Factor to multiplicate the weight (that regularizes the Hessian matrix) when H(end)>=H(end-1)
winic=w;
flagwinic=flagw;
Tup=T;

NX=size(x);NX(end+1:3)=1;
NY=size(y);
NS=size(S);
NT=size(T);

%Second order derivatives
a=[1 2 3 1 1 2 1 2 3 1 2 3 1 2 3 4 4 5 4 5 6;1 2 3 2 3 3 4 4 4 5 5 5 6 6 6 5 6 6 4 5 6];
NHe=size(a,2);

%SENSE undersampling parameters; just works for SENSE factors lower than 3
FOV=single(ones(2,2));
for n=1:2
    if NS(n)~=NY(n)
        over=NS(n)-NY(n);
        FOV(n,1)=floor(over/2);
        FOV(n,2)=ceil(over/2);
    end
end

dHe=single(zeros([NHe NT(5)]));
dH=single(zeros([NT(6) NT(5)]));
Eprev=single(zeros(NT(1:5)));
E=single(zeros(NT(1:5)));

if gpu
    Eprev=gpuArray(Eprev);E=gpuArray(E);x=gpuArray(x);y=gpuArray(y);M=gpuArray(M);A=gpuArray(A);S=gpuArray(S);
end

NRun=ceil(NT(5)/BlockSize);
NRunRem=mod(NT(5),BlockSize);
vS=cell(NRun,1);
for s=1:NRun
    if s~=NRun || NRunRem==0
        vS{s}=(s-1)*BlockSize+1:s*BlockSize;
    else
        vS{s}=(s-1)*BlockSize+1:(s-1)*BlockSize+NRunRem;
    end     
end

%Iterations
for n=1:nT    
    xPrev=x;    
    %Update the weights                   
    w(flagw==2)=w(flagw==2)/multA;
    w(flagw==1)=w(flagw==1)*multB;

    for s=1:NRun                                  
        [et,etg,eth]=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,T(:,:,:,:,vS{s},:),1,2);
        if gpu
            et{1}=gpuArray(et{1});
            for m=1:3
                etg{1}{m}=gpuArray(etg{1}{m});
                for r=2:3
                    et{r}{m}=gpuArray(et{r}{m});etg{r}{m}=gpuArray(etg{r}{m});eth{r}{m}=gpuArray(eth{r}{m});
                end
            end
            for m=1:6
                eth{1}{m}=gpuArray(eth{1}{m});
            end
        end
        [xT,xB]=sincRigidTransform(x,et,1,gpu);    
        xT=forwardApplication(xT);        
        xT=bsxfun(@minus,xT,y);
        xT=bsxfun(@times,xT,A(:,:,:,:,vS{s}));
        xTconj=conj(xT);
        Eprev(1,1,1,1,vS{s})=sum(sum(sum(sum(real(xT.*xTconj),1),2),3),4);

        [G,GB,GC]=sincRigidTransformGradient(xB,et,etg,gpu);
        for m=1:NT(6)
            G{m}=forwardApplication(G{m});
            G{m}=bsxfun(@times,G{m},A(:,:,:,:,vS{s}));
            Gconj{m}=conj(G{m});        
        end
        
        for m=1:NHe               
            GG=sincRigidTransformHessian(xB,GB,GC,et,etg,eth,m,gpu);     
            GG=forwardApplication(GG);
            GG=bsxfun(@times,GG,A(:,:,:,:,vS{s}));
            GG=real(GG.*xTconj);
            GG=GG+real(G{a(1,m)}.*Gconj{a(2,m)});
            if ~gpu
                dHe(m,vS{s})=permute(sum(sum(sum(sum(GG,1),2),3),4),[1 5 3 2 4]);
            else
                dHe(m,vS{s})=gather(permute(sum(sum(sum(sum(GG,1),2),3),4),[1 5 3 2 4]));
            end
        end
        
        for m=1:NT(6)             
            G{m}=real(bsxfun(@times,G{m},xTconj));                                            
            if ~gpu
                dH(m,vS{s})=permute(sum(sum(sum(sum(G{m},1),2),3),4),[1 5 3 2 4]);                       
            else
                dH(m,vS{s})=gather(permute(sum(sum(sum(sum(G{m},1),2),3),4),[1 5 3 2 4]));
            end
        end
    end

    MHe=single(1000*eye(NT(6)));      
    for s=1:NT(5)
        for k=1:NHe
            if a(1,k)==a(2,k)
                MHe(a(1,k),a(2,k))=dHe(k,s)+w(1,1,1,1,s);                            
            else
               MHe(a(1,k),a(2,k))=dHe(k,s);
               MHe(a(2,k),a(1,k))=dHe(k,s);
            end              
        end    
        dH(:,s)=single(double(MHe)\double(dH(:,s)));
    end    
    Tup=permute(dH,[3 4 5 6 2 1]); 
    Tup=T-Tup;    
    %Restrict to allowed ranges to prevent overshooting of Newton's method
    Tang=Tup(:,:,:,:,:,4:6);
    while ~isempty(Tang(Tang>pi))
       Tang(Tang>pi)=Tang(Tang>pi)-2*pi;
    end
    while ~isempty(Tang(Tang<-pi))
       Tang(Tang<-pi)=Tang(Tang<-pi)+2*pi;
    end
    Tup(:,:,:,:,:,4:6)=Tang;
    for m=1:3
        Ttra=Tup(:,:,:,:,:,m);
        if NX(m)>1
            while ~isempty(Ttra(Ttra>rGrid{m}(end)))
               Ttra(Ttra>rGrid{m}(end))=Ttra(Ttra>rGrid{m}(end))-NX(m);
            end
            while ~isempty(Ttra(Ttra<rGrid{m}(1)))
               Ttra(Ttra<rGrid{m}(1))=Ttra(Ttra<rGrid{m}(1))+NX(m);
            end
        end
        Tup(:,:,:,:,:,m)=Ttra;
    end   

    for s=1:NRun
        et=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tup(:,:,:,:,vS{s},:),1,0);
        if gpu
            et{1}=gpuArray(et{1});
            for m=1:3
                for r=2:3
                    et{r}{m}=gpuArray(et{r}{m});
                end
            end
        end
        xT=sincRigidTransform(x,et,1,gpu);          
        xT=forwardApplication(xT);                
        xT=bsxfun(@minus,xT,y);        
        xT=bsxfun(@times,xT,A(:,:,:,:,vS{s}));
        xTconj=conj(xT);
        E(1,1,1,1,vS{s})=sum(sum(sum(sum(real(xT.*xTconj),1),2),3),4);
    end
    
    flagw(E<Eprev)=2;
    flagw(E>=Eprev)=1;
    for s=1:size(Tup,6)
        TauxA=T(:,:,:,:,:,s);
        TauxB=Tup(:,:,:,:,:,s);
        TauxA(flagw==2)=TauxB(flagw==2);
        T(:,:,:,:,:,s)=TauxA;
    end    
    %permute(T,[5 6 1 2 3 4])    
       
    %This would diminish drifting, as suggested in ec. (28) in the paper
    if meanT
        Tmed=mean(T,5);
        etD=precomputationsSincRigidTransform(kGrid,kkGrid,rkGrid,Tmed,1,0);
        if gpu
            etD{1}=gpuArray(etD{1});
            for m=1:3
                for n=2:3
                    etD{n}{m}=gpuArray(etD{n}{m});
                end
            end
        end    
        x=sincRigidTransform(xPrev,etD,1,gpu);   
        x=M.*x;
        T=bsxfun(@minus,T,Tmed);    
    end    
end
if gpu
    x=gather(x);
end
    
function x=forwardApplication(x)
    x=bsxfun(@times,x,S);
    for r=1:2
        x=sense(x,r,NS(r),NY(r),FOV(r,:));
    end
    for r=1:2
        x=fftGPU(x,r,gpu);
    end
end

end
