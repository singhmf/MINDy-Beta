function[Out]=MINDy_Simple_BIG(Dat,ParcMask,varargin)
%% Just Input Data as region x time (a single matrix or cell of matrices)
%% Output.Param={Wsparse,A,b,c,Wfull,D}
doCell='y';

ChosenPARSTR

ParStr.wPC=50;
ParStr.BatchSz=200;ParStr.NBatch=2000;
ParStr.H1min=5;ParStr.H1max=7;ParStr.H2min=.7;ParStr.H2max=1.3;ParStr.H1Rate=.1;ParStr.H2Rate=.1;

if isempty(varargin)
    doPreProc='y';
else
    doPreProc=varargin{1};
end


RegNames={'SpScale','SpDiag','SpPls1','SpPls2','L2SpPlsEN'};
RegScale=3;%numel(ParcMask)/419;
for ii=1:numel(RegNames)
    ParStr.(RegNames{ii})=ParStr.(RegNames{ii})*RegScale;
end

if ~iscell(Dat)
Dat={Dat};
end
for i=1:numel(Dat)
Dat{i}=Dat{i}(:,~isnan(sum(Dat{i},1)));
end
Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);

if strcmpi(doPreProc,'y')
Dat=MINDy_RestingPreProcInterp(Dat,Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR);
Dat=cellfun(@(xx)(zscore(xx(:,20:(end-20))')'),Dat,'UniformOutput',0);
Dat=cellfun(@(xx)(convn(xx,[1 1],'valid')),Dat,'UniformOutput',0);
Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);
end


dDat=cellfun(@(xx)(convn(xx,[1 0 -1]/2,'valid')),Dat,'UniformOutput',0);
Dat=cellfun(@(xx)(xx(:,1:end-2)),Dat,'UniformOutput',0);

Out=MINDy_BIG(Dat,dDat,ParcMask,Pre,ParStr,doCell);

Out.Pre=Pre;Out.ParStr=ParStr;
end