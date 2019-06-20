function[OutStruct]=MINDy_MakeHRFpoly(X0,Pre,ParStr,NormDeriv)

ParamRes=10;

StartDrop=39;
EndDrop=40;


if ~iscell(X0)
    X0={X0};
end

FiltAmp=Pre.FiltAmp;

PolyC=cell(1,numel(X0));
dPolyC=cell(1,numel(X0));
OrigTime=cell(1,numel(X0));



for iXcell=1:numel(X0)
X=X0{iXcell};
X=zscore(X')';
%find(isnan(X))    
for iParc=1:size(X,1)
        gg=find(abs(X(iParc,:))>FiltAmp);
%        BadFrames{i,iParc}=gg;
%% Remove bad end points (can't interpolate)
gg(ismember(gg,[1 size(X,2)]))=[];
        goodFrames=setdiff(1:size(X,2),gg);
        X(iParc,gg)=interp1(goodFrames,X(iParc,goodFrames),gg);
end
%find(isnan(X))

%% ConvRes=TR
a1=6;
a2=16;
b1=1;
b2=1;
c=1/6;
nX=size(X,1);

ConvRes=Pre.TR;

h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));
ConvVec=[0 h(ConvRes*(1:(size(X,2)-1)))];

ChangeVec=linspace(ParStr.H1min,ParStr.H1max,ParamRes);
ChangeVec2=linspace(ParStr.H2min,ParStr.H2max,ParamRes);

Outbase=deconvwnr(X(1,:),ConvVec,Pre.ConvLevel);


Out=nan(numel(ChangeVec),numel(ChangeVec),nX,numel(Outbase));
nC=[numel(ChangeVec2) numel(ChangeVec2)];
Coeff1=nan(nC);Coeff2=nan(nC);
%ConvAll=cell(nC);
for i=1:numel(ChangeVec)
    for j=1:numel(ChangeVec)
%a1=6;
a2=16;
%b1=1;
b2=1;
c=1/6;
b1=ChangeVec2(j);
a1=ChangeVec(i);
h=@(t)((t.^(a1-1).*exp(-b1*t)*(b1^a1))/gamma(a1)-c*((t.^(a2-1).*exp(-b2*t)*b2^a2)/gamma(a2)));
ConvVec=[0 h(ConvRes*(1:(size(X,2)-1)))];

    
%ConvAll{i,j}=(ConvVec(1:EndDrop)/std(ConvVec(1:EndDrop)))';

%Out(i,j,:,:)=(deconvwnr(X,ConvVec,Pre.ConvLevel));
Out(i,j,:,:)=zscore(deconvwnr(X,ConvVec,Pre.ConvLevel)')';
disp([i j numel(ChangeVec)])
Coeff1(i,j)=a1;
Coeff2(i,j)=b1;
    end
end
Out=(Out-mean(Out,4))./std(Out,[],4);

PolyTerms={[0 0],[1 0],[0 1],[2 0],[1 1],[0 2],[3 0],[2 1],[1 2],[0 3]};
for i=1:numel(PolyTerms)
    PolyReg(i,:)=((Coeff1(:)).^PolyTerms{i}(1)).*((Coeff2(:)).^PolyTerms{i}(2)); %#ok<AGROW>
end

Out=Out(:,:,:,(1+StartDrop):end-EndDrop);

Poly=nan(nX,size(Out,4)-1,numel(PolyTerms));
dPoly=nan(nX,size(Out,4)-1,numel(PolyTerms));
for iP=1:nX
    disp([iP nX])
%% Do Original
Out2=permute(squeeze(Out(:,:,iP,1:end-1)),[3 1 2]);
Poly(iP,:,:)=(PolyReg'\(Out2(:,:))')';
%% Do Derivatives
Out3=permute(squeeze(Out(:,:,iP,2:end)-Out(:,:,iP,1:end-1)),[3 1 2]);
%Out3=RegressOut(Out3,Out2,1);
if strcmpi(NormDeriv,'y')
Out3=Out3./std(Out3,[],1);
end
dPoly(iP,:,:)=(PolyReg'\(Out3(:,:))')';
end
dPolyC{iXcell}=dPoly;
PolyC{iXcell}=Poly;
OrigTime{iXcell}(1,:)=(1+StartDrop):(size(X,2)-(1+EndDrop));
OrigTime{iXcell}(2,:)=iXcell;
end
Poly=[PolyC{:}];
dPoly=[dPolyC{:}];
PolyTerms=Uncellfun(@(xx)(xx'),PolyTerms);
PolyMat=[PolyTerms{:}];

OutStruct.OrigTime=OrigTime;
OutStruct.StartDrop=StartDrop;
OutStruct.EndDrop=EndDrop;

OutStruct.Poly=Poly;
OutStruct.dPoly=dPoly;
OutStruct.PolyMat=PolyMat;
end
