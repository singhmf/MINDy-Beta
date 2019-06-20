function [zz,offDiagonals]=FitDiagnostics(X1,X2,varargin)
%% varargin indicates whether to plot
if ~isempty(varargin)
    isPlot=varargin;
else
    isPlot='y';
end


%% Creates Plots, gives diagnostics for repeated fits.
zz=zeros(1,9);
if strcmpi(isPlot(1),'y')
figure;
end

zz(7)=corr(X1.Param{1}(:),X2.Param{5}(:));
zz(8)=corr(X1.Param{5}(:),X2.Param{1}(:));
if (~isempty(X1.Param{3}))&&(~isempty(X2.Param{3}))
zz(9)=corr(X1.Param{3}(:,1),X2.Param{3}(:,1));
end
%% Correction from D^2 back to D
if isfield(X1,'In')
        X1.Param{6}=sqrt(X1.Param{6}-X1.In{1});
        X2.Param{6}=sqrt(X2.Param{6}-X2.In{1});
else
disp('Not using Dmin')    
        X1.Param{6}=sqrt(X1.Param{6});
        X2.Param{6}=sqrt(X2.Param{6});
end
        
Names={'bW','A','Thresh','C','wW','D'};
for i=1:6
    if i==3
    zz(i)=corr(X1.Param{i}(:,2),X2.Param{i}(:,2));
    else        
    zz(i)=corr(X1.Param{i}(:),X2.Param{i}(:));
    end
    if strcmpi(isPlot(1),'y')&&~isnan(zz(i))
    subplot(2,3,i)
if i==3
    scatter(X1.Param{i}(:,2),X2.Param{i}(:,2));
elseif or(i==1,i==5)        
    scatter(X1.Param{i}(:),X2.Param{i}(:),.5,'filled');
else
    scatter(X1.Param{i}(:),X2.Param{i}(:))
end

if i==1
    offDiagonals(1)=corr(OffDiag(X1.Param{1}),OffDiag(X2.Param{1}));
elseif i==5
    offDiagonals(2)=corr(OffDiag(X1.Param{5}),OffDiag(X2.Param{5}));
end

title(strcat(Names{i},' r= ',num2str(zz(i))))
EqAx;
lsline
hold on
plot(xlim,xlim,'-k')
    end
end

cc1=numel(X1.Param);cc2=numel(X2.Param);

for j=1:min([cc1 cc2])
    if (isnumeric(X1.Param{j}))&&(isnumeric(X2.Param{j}))&&(~isempty(X1.Param{j}))&&(~isempty(X2.Param{j}))
    zz(9+j)=corr(X1.Param{j}(:),X2.Param{j}(:));
    end
end

end