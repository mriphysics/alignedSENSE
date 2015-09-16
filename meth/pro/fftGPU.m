function x=fftGPU(x,m,gpu)

%FFTGPU   Optimal arrangement of GPU-based FFT computation
%   X=FFTGPU(X,M,GPU) applies FFT's over rearranged arrays so that the FFT 
%   is applied across the first dimension. This has accelerated GPU 
%   computations in our setting
%   X is the array on which to apply the FFT
%   M is the direction across which to apply the FFT
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   It returns:
%   X, the FFT-transformed array
%

perm=1:ndims(x);
perm(1)=m;
perm(m)=1;

if ~gpu
    x=fft(x,[],m);
elseif m==1
    if size(x,1)~=1
        x=fft(x);
    end
else
    x=permute(x,perm);
    if size(x,1)~=1
        x=fft(x);
    end
    x=permute(x,perm);
end
