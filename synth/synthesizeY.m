function y=synthesizeY(xGT,TGT,S,A,kGrid,rkGrid)

%SYNTHESIZEY   Synthesizes motion corrupted data
%   Y = SYNTHESIZEY(XGT,TGT,S,A,NS,THETA,ENCMETH,XKGRID,KGRID) synthesizes
%   motion corrupted data according to the forward model in Fig. 1 of the
%   paper,
%   XGT is the ground truth image
%   TGT is a cell array of cell array of #TestedShots x #MotionLevels 
%   ground truth transforms
%   S is the coil-array sensitivity map
%   A is a cell array of #TestedShots x #EncodingMethods containing the 
%   sampling scheme
%   KGRID is a 1x3 cell array with the spectral coordinates along each axis
%   RKGRID is a 2x3 cell array where the first dimension of the cell
%   indexes the shearing orientation and the second indexes the 
%   spatio-spectral axes pair used for that particular shearing (given by 
%   [1 3 2;2 1 3]). See ec. 4 of the paper for further details.
%

for s=1:length(TGT)
    for v=1:length(TGT{s})  
        %Tx        
        et=precomputationsSincRigidTransform(kGrid,[],rkGrid,TGT{s}{v},1,0);
        xS=sincRigidTransform(xGT,et,1,0);
        %STx
        xS=bsxfun(@times,xS,S);
        %FSTx
        for m=1:2
            xS=fftGPU(xS,m,0);
        end
        %AFSTx
        for e=1:length(A{s})
            y{s}{v}{e}=sum(bsxfun(@times,xS,A{s}{e}),5);
        end        
    end
end
