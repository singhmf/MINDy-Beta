function[Out]=SplitMINDy(ooP,NetInd)
UnNet=unique(NetInd);
Out=cell(1,max(UnNet));

%% Automatically generates rest functions
%% Assumes that only weight matrices are square (so nParc>2)

for i=UnNet
    Out{i}=ooP;
    for j=1:numel(ooP.Param)
        if ~isempty(ooP.Param{j})
            if size(ooP.Param{j},2)~=size(ooP.Param{j},1)
                Out{i}.Param{j}=ooP.Param{j}(NetInd==i,:);
            else
                Out{i}.Param{j}=ooP.Param{j}(NetInd==i,NetInd==i);
            end
        end
    end
    Out{i}=MakeRestFunctions(Out{i},'z','n');
end
end
    