function TGT=synthesizeT(NS,theta,randomMotion)

%SYNTHESIZET   Synthesizes a set of rigid transforms affecting the 
%   different shots.
%   TGT = SYNTHESIZET(NS,THETA,RANDOMMOTION)
%   synthesizes a set of rigid transforms according to a rotation parameter
%   NS is a vector with each element indicating the number of transforms 
%   (corresponding to the number of shots)
%   THETA is a vector with each element specifying the rotation parameter 
%   RANDOMMOTION specifies whether motion has to be generated in a random
%   manner by sampling from a uniform distribution in [-theta/2,theta/2]
%   (default) or, provided NS=2 and RANDOMMOTION=0, it should be generated
%   deterministically as two motion states at -theta/2 and theta/2.
%   It returns a cell array of LENGTH(NS) x LENGTH(THETA) transforms
%

if ~exist('randomMotion','var');randomMotion=1;end

TGT=cell(length(NS),length(theta));
for S=1:length(NS)%Shots
    for V=1:length(theta)%Rotation levels
        %Cell array of transforms
        TGT{S}{V}=single(zeros([1 1 1 1 NS(S) 6]));%6th dimension for transform parameters: 1-3 translations / 4-6 rotations
        if length(NS)==1 && NS(S)==2 && ~randomMotion
            TGT{S}{V}(1,1,1,1,1,4)=theta(V)*pi/180;
        else
            TGT{S}{V}(1,1,1,1,:,4)=pi*theta(V)*(rand([1 1 1 1 NS(S)])-0.5)/180;
        end                
        TGTm=mean(TGT{S}{V},5);
        TGT{S}{V}=bsxfun(@minus,TGT{S}{V},TGTm);%0-mean due to considerations above ec (23) in the paper
    end
end
