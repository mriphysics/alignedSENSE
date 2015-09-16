function [G,GB,GC]=sincRigidTransformGradient(xB,et,etg,gpu)

%SINCRIGIDTRANSFORMGRADIENT   Computes the gradient of sinc-interpolated 
%3D rigid transforms
%   [G,GB,GC]=SINCRIGIDTRANSFORMGRADIENT(xB,ET,ETG,GPU) obtains the 
%   gradient of the transform of the images
%   XB contains images before the first, second and third rotations and 
%   before the translation
%   ET are the diagonal operators of the transform as given by
%   precomputationSincRigidTransform
%   ETG are the diagonal operators of the transform gradient as given by
%   precomputationSincRigidTransform
%   GPU is a flag that determines whether to use gpu (1) or cpu (0) 
%   computation
%   It returns:
%   G, the gradient of the transformed image
%   GB, the gradient of the first, second and third rotations before 
%   applying the translation
%   GC, the gradient of the first, first and second rotations before
%   applying the second, third and third rotations respectively
%

GB=[];
GC=[];

%Translation parameters
for m=1:3    
    x{1}=bsxfun(@times,xB{4},etg{1}{m});
    for n=1:3       
        x{1}=ifftGPU(x{1},n,gpu);
    end    
    G{m}=x{1};
end

%First rotation
x{1}=bsxfun(@times,xB{1},et{2}{1});
x{2}=bsxfun(@times,xB{1},etg{2}{1});
for m=1:2
    x{m}=ifftGPU(x{m},1,gpu);
    x{m}=fftGPU(x{m},2,gpu);
end
x{2}=bsxfun(@times,etg{3}{1},x{1})+bsxfun(@times,et{3}{1},x{2});
x{1}=bsxfun(@times,et{3}{1},x{1});    
for m=1:2
    x{m}=ifftGPU(x{m},2,gpu);       
    x{m}=fftGPU(x{m},1,gpu);
end
x{1}=bsxfun(@times,etg{2}{1},x{1})+bsxfun(@times,et{2}{1},x{2});
x{1}=ifftGPU(x{1},1,gpu);

x{1}=fftGPU(x{1},3,gpu);
GC{1}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftGPU(x{1},3,gpu);
x{1}=fftGPU(x{1},1,gpu);
x{1}=bsxfun(@times,x{1},et{3}{2});
x{1}=ifftGPU(x{1},1,gpu);
x{1}=fftGPU(x{1},3,gpu);
x{1}=bsxfun(@times,x{1},et{2}{2});
x{1}=ifftGPU(x{1},3,gpu);

x{1}=fftGPU(x{1},2,gpu);
GC{2}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftGPU(x{1},2,gpu);
x{1}=fftGPU(x{1},3,gpu);
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftGPU(x{1},3,gpu);
x{1}=fftGPU(x{1},2,gpu);
x{1}=bsxfun(@times,x{1},et{2}{3});

for m=1:2:3
    x{1}=fftGPU(x{1},m,gpu);
end
GB{1}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1    
    x{1}=ifftGPU(x{1},m,gpu);   
end
G{4}=x{1};


%Second rotation
x{1}=bsxfun(@times,xB{2},et{2}{2});
x{2}=bsxfun(@times,xB{2},etg{2}{2});
for m=1:2
    x{m}=ifftGPU(x{m},3,gpu);
    x{m}=fftGPU(x{m},1,gpu);
end
x{2}=bsxfun(@times,etg{3}{2},x{1})+bsxfun(@times,et{3}{2},x{2});
x{1}=bsxfun(@times,et{3}{2},x{1});
for m=1:2
    x{m}=ifftGPU(x{m},1,gpu);
    x{m}=fftGPU(x{m},3,gpu);
end
x{1}=bsxfun(@times,etg{2}{2},x{1})+bsxfun(@times,et{2}{2},x{2});
x{1}=ifftGPU(x{1},3,gpu);

x{1}=fftGPU(x{1},2,gpu);
GC{3}=x{1};
x{1}=bsxfun(@times,x{1},et{2}{3});
x{1}=ifftGPU(x{1},2,gpu);
x{1}=fftGPU(x{1},3,gpu);
x{1}=bsxfun(@times,x{1},et{3}{3});
x{1}=ifftGPU(x{1},3,gpu);
x{1}=fftGPU(x{1},2,gpu);
x{1}=bsxfun(@times,x{1},et{2}{3});

for m=1:2:3
    x{1}=fftGPU(x{1},m,gpu);
end
GB{2}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1    
    x{1}=ifftGPU(x{1},m,gpu);    
end
G{5}=x{1};


%Third rotation
x{1}=bsxfun(@times,xB{3},et{2}{3});
x{2}=bsxfun(@times,xB{3},etg{2}{3});
for m=1:2
    x{m}=ifftGPU(x{m},2,gpu);
    x{m}=fftGPU(x{m},3,gpu);
end
x{2}=bsxfun(@times,etg{3}{3},x{1})+bsxfun(@times,et{3}{3},x{2});
x{1}=bsxfun(@times,et{3}{3},x{1});
for m=1:2
    x{m}=ifftGPU(x{m},3,gpu);
    x{m}=fftGPU(x{m},2,gpu);
end
x{1}=bsxfun(@times,etg{2}{3},x{1})+bsxfun(@times,et{2}{3},x{2});

for m=1:2:3
    x{1}=fftGPU(x{1},m,gpu);        
end
GB{3}=x{1};
x{1}=bsxfun(@times,x{1},et{1});
for m=3:-1:1       
    x{1}=ifftGPU(x{1},m,gpu);        
end
G{6}=x{1};