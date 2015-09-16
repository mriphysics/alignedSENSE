function x=isense(x,m,NS,NY,FOV)

%ISENSE   Inverse SENSE folding operation
%   X=ISENSE(X,M,NS,NY,FOV) applies the inverse SENSE folding operation
%   along a given dimension of the input array (up to 8D arrays)
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
    if NS<=3*NY%For a SENSE factor less than 3             
        x=[x(1+FOV(1):end,:,:,:,:,:,:,:);x;x(1:end-FOV(2),:,:,:,:,:,:,:)];            
    else
        error('The reconstruction code does not work for SENSE factors larger than 3');
    end   
    if m~=1
        x=permute(x,perm);
    end
end
