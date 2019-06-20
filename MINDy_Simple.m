function[Out]=MINDy_Simple(Dat,varargin)
%% Just Input Data as region x time (a single matrix or cell of matrices)
%% Output.Param={Wsparse,A,b,c,Wfull,D}

%% Optional input argument is whether to do preprocessing (filtering and deconvolution)
if isempty(varargin)
    doPreProc='y';
else
    doPreProc=varargin{1};
end

ChosenPARSTR
ParStr.BatchSz=300;ParStr.NBatch=5000;
doRobust='n';

ParStr.H1min=5;ParStr.H1max=7;ParStr.H2min=.7;ParStr.H2max=1.3;ParStr.H1Rate=.1;ParStr.H2Rate=.1;

if ~iscell(Dat)
Dat={Dat};
end
for i=1:numel(Dat)
Dat{i}=Dat{i}(:,~isnan(sum(Dat{i},1)));
end


Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);

%% Preprocessing (if desired)
if strcmpi(doPreProc,'y')
Dat=MINDy_RestingPreProcInterp(Dat,Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR);
Dat=cellfun(@(xx)(zscore(xx(:,20:(end-20))')'),Dat,'UniformOutput',0);
end

%% Compute Finite-Difference Derivatives
dDat=cellfun(@(xx)(convn(xx,[1 -1],'valid')),Dat,'UniformOutput',0);
Dat=cellfun(@(xx)(xx(:,1:end-1)),Dat,'UniformOutput',0);

%% Fit MINDY models
Out=MINDy_Base(Dat,dDat,Pre,ParStr);

%% Censor filtered data-points for subsequent analyses/adjustment
[X1,dX1]=MINDy_Censor(Out,Dat,dDat);

%% Global regression on W and D to get rid of global regularization bias
Out=MINDy_Inflate(Out,X1,dX1,doRobust);
Out=MakeMINDyFunction(Out);

%% Recalculate goodness-of-fit
Out.Corr=DiagCorr(Out.FastFun([X1{:}])',[dX1{:}]');
Out=rmfield(Out,{'FastFun','Tran'});
Out.Pre=Pre;Out.ParStr=ParStr;
end