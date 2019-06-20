function[Out]=MINDy_Extract_BIGmat(ooP,Ind1,Ind2)
%% Extract subblock of the weight matrix from MINDy output ooP along indices
%% [Ind1] x [Ind2] (both vectors
ParcMask=ooP.ParcMask;jScale=ooP.jScale;
if ~isfield(ooP,'kScale')
    %% For MINDy BIG
Out=jScale{1}(Ind1).*(ParcMask(:,Ind1)'*(ooP.ParcCount.*ooP.Param{5})*ParcMask(:,Ind2)).*(jScale{2}(Ind2)');
else
    %% For MINDy BIGGER
%% Use j-Scaling for Sparse part
kScale=ooP.kScale;
jPart=jScale{1}(Ind1).*(ParcMask(:,Ind1)'*(ooP.ParcCount.*ooP.Param{1})*ParcMask(:,Ind2)).*jScale{2}(Ind2);
WK=ooP.Param{5}-ooP.Param{1};
%% Use k-Scaling for Low-Rank part
Out=jPart+kScale{1}(Ind1).*(ParcMask(:,Ind1)'*(ooP.ParcCount.*WK)*ParcMask(:,Ind2)).*kScale{2}(Ind2);
end
end