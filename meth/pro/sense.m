function x=sense(x,m,NS,NY,FOV)

%ISENSE   SENSE folding operation
%   X=SENSE(X,M,NS,NY,FOV) applies the SENSE folding operation along a 
%   given dimension of the input array (up to 8D arrays)
%   X is the array on which to operate
%   M is the direction across which to operate
%   NS is the spatial size of the array
%   NY is the spectral size of the array
%   FOV indicates the size of the slab to be folded
%   It returns:
%   X, the inverse SENSE folded array
%

perm=1:ndims(x);
perm(1)=m;
perm(m)=1;

if NS~=NY%SENSE
    if m~=1
        x=permute(x,perm);
    end    
    xb=x(FOV(2)+1:end-FOV(1),:,:,:,:,:,:,:);       
    xb(1:FOV(1),:,:,:,:,:,:,:)=xb(1:FOV(1),:,:,:,:,:,:,:)+x(end-FOV(1)+1:end,:,:,:,:,:,:,:);
    xb(end-FOV(2)+1:end,:,:,:,:)=xb(end-FOV(2)+1:end,:,:,:,:)+x(1:FOV(2),:,:,:,:);
    x=xb;
    if m~=1
        x=permute(x,perm);
    end
end
