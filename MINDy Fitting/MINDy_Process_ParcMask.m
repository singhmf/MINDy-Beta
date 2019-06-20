function[Out]=MINDy_Process_ParcMask(ParcMask)
%% Change dlabel format to a logical vector
nParc=max(ParcMask);
Out=false(nParc,numel(ParcMask));
for ii=1:nParc
    Out(ii,:)=ParcMask==ii;
end
end