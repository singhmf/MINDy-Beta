function[Out]=MINDy_MakeStarts(X,ParStr)
nBack=ParStr.nBack;
nStep=ParStr.nStep;
CellLength=Uncellfun(@(xx)(size(xx,2)),X);
CellLength=[0 [CellLength{:}]];
nC=numel(X);
Out=cell(1,nC);
for i=1:nC
Out{i}=CellLength(i)+(nBack+1):(CellLength(i+1)-nStep);
end
Out=[Out{:}];
end
