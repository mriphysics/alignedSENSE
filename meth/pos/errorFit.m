function E=errorFit(x,T,S,A,y,kGrid,rkGrid)

%ERRORFIT   Measures the error for the least-squares reconstruction
%formulation
%   E=ERRORFIT(X,T,S,A,Y,KGRID,RKGRID) returns the value of the objective 
%   function of the reconstruction (assuming no SENSE is used)
%   X is the reconstructed image
%   T are the estimated transforms for each shot
%   S is the coil-array sensitivity map
%   A is the sampling scheme
%   KGRID is a grid of points in the spectral domain
%   RKGRID is a grid of points in the spatial-spectral domain
%   It returns:
%   E, the energy value
%

%Tx
et=precomputationsSincRigidTransform(kGrid,[],rkGrid,T,1,0);
x=sincRigidTransform(x,et,1,0);
%STx
x=bsxfun(@times,x,S);
%FSTx
for m=1:2
    x=fftGPU(x,m,0);
end
%AFSTx
x=bsxfun(@times,x,A);
y=bsxfun(@times,y,A);
%AFSTx-y
x=x-y;
E=sum(sum(sum(sum(sum(real(x.*conj(x)))))));
