function[Out]=MINDy_FULL(Dat,Pre,ParStr,doSplit,doPreProc,doInflate)
%% THIS IS THE MASTER FUNCTION FOR MINDy in the first paper!!!

%% Dat is input as a cell of timeseries (region x time) or a single time-series
%% Pre is a structure of parameters for doing preprocessing (includes TR, deconvolution filtering etc.)
%% ParStr is a structure of hyperparameters for model fitting
%% doSplit denotes whether to use   dx=x(t+1)-x(t) (doSplit='n') or 
%%                                  dx=(x(t+2)-x(t))/2 (doSplit='y')
%% doPreProc=whether to do deconvolution/spike removal
%% doInflate=whether to do global readjustment of weights vs. decay post-fitting




if ~iscell(Dat)
    Dat={Dat};
end



%% Preprocessing (if desired)
if strcmpi(doPreProc,'y')
for i=1:numel(Dat)
Dat{i}=Dat{i}(:,~isnan(sum(Dat{i},1)));
end
Dat=cellfun(@(xx)(zscore(xx')'),Dat,'UniformOutput',0);
Dat=MINDy_RestingPreProcInterp(Dat,Pre.FiltAmp,Pre.ConvLevel,Pre.DownSamp,Pre.TR);
Dat=cellfun(@(xx)(zscore(xx(:,50:(end-50))')'),Dat,'UniformOutput',0);

end

%% Calculating derivatives
if ischar(doSplit)
    if strcmpi(doSplit(1),'y')
        DerivDat=Uncellfun(@(xx)(convn(xx,[1 0 -1]/2,'valid')),Dat);
        Dat=Uncellfun(@(xx)(xx(:,1:end-2)),Dat);
    elseif strcmpi(doSplit(1),'n')        
        DerivDat=Uncellfun(@(xx)(convn(xx,[1 -1],'valid')),Dat);
        Dat=Uncellfun(@(xx)(xx(:,1:end-1)),Dat);
    else
        error('doSplit should be y/n or numeric')
    end
else
        DerivDat=Uncellfun(@(xx)(convn(xx,[1 -1]*doSplit,'valid')),Dat);
        Dat=Uncellfun(@(xx)(xx(:,1:end-1)),Dat);
end 
Out=MINDy_Base(Dat,DerivDat,Pre,ParStr);
[X1,dX1]=ModularCensor(Out,Dat,DerivDat);
if strcmpi(doInflate(1),'y')
    if ~isfield(ParStr,'doRobust')        
    Out=MINDy_Inflate(Out,X1,dX1,'y');
    else
    Out=MINDy_Inflate(Out,X1,dX1,ParStr.doRobust);
    end
end
Out=MakeMINDyFunction(Out);
Out.Corr=DiagCorr(Out.FastFun([X1{:}])',[dX1{:}]');
end