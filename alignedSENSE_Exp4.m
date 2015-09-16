%ALIGNEDSENSE_EXP4 script performs the experiment included in Section III.D
%of the manuscript ''Sensitivity encoding for aligned multishot magnetic 
%resonance reconstruction'', L. Cordero-Grande, R. P. A. G. Teixeira, E. J.
%Hughes, J. Hutter, A. N. Price, and J. V. Hajnal

pathOu='./data';%Data path
gpu=1;%0->Use CPU / 1->Use GPU
debug=1;%0/1/2 increasing the amount of debug info provided 

%Load synthetic data:
% - Ground truth image xGT of size 128x128
% - Sensitivity maps S of size 128x128x1x32
nameX=sprintf('%s/xGT',pathOu);load(nameX);


%%%%%%%%%%%%%%%%%%%%
%%%DATA SYNTHESIS%%%
%%%%%%%%%%%%%%%%%%%%

%Transform generation
theta=[2 5 10];%Displacement
NS=[2 4 8 16 32 64 128];%Number of shots
randomMotion=1;%0->deterministic motion / 1->random motion
TGT=synthesizeT(NS,theta,randomMotion);

%Grid generation
N=size(xGT);N(end+1:3)=1;%Image size
[kGrid,kkGrid,rGrid,rkGrid]=generateGrids(N);

%Encoding generation
EncMeth{1}='LinPar';
[A W]=generateEncoding(N,kGrid,NS,EncMeth);

%Isotropic resolution
xGT=ifft(ifft(W.*fft(fft(xGT,[],1),[],2),[],2),[],1);

%Coil array compression to accelerate
perc=0.99;%Percentage of energy retained
S=coilArrayCompression(S,[],perc,gpu);

%Data generation
y=synthesizeY(xGT,TGT,S,A,kGrid,rkGrid);


%%%%%%%%%%%%%%%%%%%%
%%%RECONSTRUCTION%%%
%%%%%%%%%%%%%%%%%%%%

%Solver parameters
testing=0;%0->Full convergence / 1->Partial convergence (mainly to debug)
if testing
    toler=1e-6;%Threshold for convergence
else
    toler=1e-12;
end
nExtern=10000;%Maximum number of iterations of the joint method
nX=3;%Number of iterations of CG
nT=1;%Number of iterations of Newton's method
winic=1;%initial w in ec. (22) of the paper

%Preconditioner for CG
SH=conj(S);
reg=0.001;
Precond=sum(real(SH.*S),4);
Precond=(Precond+reg).^(-1);

%Spatial mask (which, if available, could be used to constrain solution).
%Disabled here. To enable, the voxels where the solution is constrained to
%be zero should be zero in the mask
M=ones(size(xGT));

typeRec{1}='Multishot SENSE';
typeRec{2}='Aligned Multishot SENSE';


for eT=1:2%1->Multishot SENSE reconstruction / 2->Aligned multishot SENSE reconstruction
    estT=eT-1;%Flag to indicate whether to estimate T
    for s=1:length(NS)%Shots
        BlockSize=NS(s);%# shots to be processed in a chunk; reduce for low-memory systems
        for v=1:length(theta)%Motion levels
            for e=1:length(EncMeth)%Encoding methods 
                %Initialization
                x=zeros(size(xGT));
                T=zeros(size(TGT{s}{v}));
                NT=size(T);
                w=winic*ones(NT(1:5));
                flagw=zeros(NT(1:5));flagwprev=zeros(NT(1:5)); 

                %Precomputations
                yX=ifft(ifft(y{s}{v}{e},[],1),[],2);%F'y
                maximNormalize=max(max(max(max(abs(yX)))));
                yIn=y{s}{v}{e}/maximNormalize;  
                
                if debug>0;fprintf('Solving for %s with %d shots, theta=%.0f and encoding %s\n',typeRec{eT},NS(s),theta(v),EncMeth{e});end
                for n=1:nExtern
                    %Solve for x
                    if debug==2;tstart=tic;end
                    xant=x;%To check convergence    
                    x=solveX(x,yIn,M,T,S,SH,Precond,A{s}{e},kGrid,rkGrid,nX,gpu,BlockSize);            
                    if debug==2;telapsed=toc(tstart);fprintf('Solving x time: %f\n',telapsed);end
                    %Solve for T    
                    if estT
                        if debug==2;tstart=tic;end    
                        flagwprev=flagw;%To check convergence
                        [T,x,w,flagw]=solveT(x,yIn,M,T,S,A{s}{e},kGrid,kkGrid,rGrid,rkGrid,nT,w,flagw,gpu,BlockSize);
                        w(w<1e-4)=2*w(w<1e-4);w(w>1e16)=w(w>1e16)/2;%To avoid numeric instabilities
                        if debug==2;telapsed=toc(tstart);fprintf('Solving T time: %f\n',telapsed);end
                    end               
                    %Measure errors                    
                    xaux=x*maximNormalize-xGT;
                    if estT==1
                        errX{v}{e}{s}(n)=sum(sum((xaux).*conj(xaux)))/(N(1)*N(2));
                        errT{v}{e}{s}(n)=max(max(abs(T-TGT{s}{v})));
                        errF{v}{e}{s}(n)=errorFit(x,T,S,A{s}{e},yIn,kGrid,rkGrid);
                    end
                    xant=x-xant;
                    xant=real(xant.*conj(xant));
                    xant=max(max(max(xant)));
                    if debug==2;fprintf('Iteration %04d - Error %0.2g \n',n,xant);end
                    if xant<toler && sum(flagwprev~=1)>0
                        if debug>0;fprintf('Convergence reached at iteration %04d \n',n);end
                        break
                    elseif n==nExtern
                        if debug>0;fprintf('Maximum number of iterations reached without convergence\n');end                    
                    end
                end
                xEst{eT}{v}{e}{s}=x*maximNormalize;
                if estT==1
                    TEst{v}{e}{s}=T;
                end 
            end
        end
    end
end
nameExp=sprintf('%s/Exp4',pathOu);
save(nameExp,'xEst','xGT','TEst','TGT','A','errX','errT','errF','theta','NS','EncMeth');
