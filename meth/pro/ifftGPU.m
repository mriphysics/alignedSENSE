function x=ifftGPU(x,m,gpu)

%IFFTGPU   Optimal arrangement of GPU-based IFFT computation
%   X=IFFTGPU(X,M,GPU) applies IFFT's over rearranged arrays so that the
%   IFFT is applied across the first dimension. This has accelerated GPU 
%   computations in our setting
%   X is the array on which to apply the IFFT
%   M is the direction across which to apply the IFFT
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   It returns:
%   X, the IFFT-transformed array
%

perm=1:ndims(x);
perm(1)=m;
perm(m)=1;

if ~gpu
    x=ifft(x,[],m);
elseif m==1
    if size(x,1)~=1
        x=ifft(x);
    end
else
    x=permute(x,perm);
    if size(x,1)~=1
        x=ifft(x);
    end
    x=permute(x,perm);
end
