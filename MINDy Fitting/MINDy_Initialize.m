%% Initialize the MINDy Software
%% For non-HRF clean data

%% Set Defaults
%ParStr=MINDy_SetDefaults(ParStr);

if strcmpi(MINDy_type,'HRF')
SpecificRates={'ARate','DRate','wk1Rate','wk2Rate','H1Rate','H2Rate'};
else
SpecificRates={'ARate','DRate','wk1Rate','wk2Rate'};
end


for iSR=1:numel(SpecificRates)
    ParStr.(SpecificRates{iSR})=ParStr.(SpecificRates{iSR})*ParStr.Rate;
end
ParStr.W1Rate=ParStr.Rate;
mu=ParStr.Dec;v=ParStr.Dec2;


%% Extract
MINDy_ExtractParStr

[X,dX,Out]=MINDy_Filter(X,dX,ParStr);
nX=size(X,1);
nT=size(X,2);

jRand=RandScale;
if (~strcmpi(MINDy_type,'BIG'))&&(~strcmpi(MINDy_type,'BIGGER'))
nParc=nX;
elseif strcmpi(MINDy_type,'BIG')
nParc=max(ParcMask);

%% Initialize vertex scaling
J1=abs(1+randn(nX,1)/jRand);J2=abs(1+randn(nX,1)/jRand);
%% NADAM for vertex scaling
mJ1=zeros(size(J1));mJ2=mJ1;nJ1=mJ1;nJ2=mJ1;
else
nParc=max(ParcMask);
%% Initialize vertex scaling
J1=abs(1+randn(nX,1)/jRand);J2=abs(1+randn(nX,1)/jRand);
K1=abs(1+randn(nX,1)/kRand);K2=abs(1+randn(nX,1)/kRand);
%% NADAM for vertex scaling
mJ1=zeros(size(J1));mJ2=mJ1;nJ1=mJ1;nJ2=mJ1;
mK1=zeros(size(K1));mK2=mK1;nK1=mK1;nK2=mK1;
end

%% Initialize Parameter Values
B1(:,1)=Bmin;
A=.5+(1+abs(randn(nX,1))).^(.5)*RandScale+AdiffMin+B1(:,1).^2/4;
B1(:,2)=0;
W1=RandScale*randn(nParc);
D=Dstart+5*sqrt(RandScale*abs(randn(nX,1)));
Wk1=RandScale*randn(nParc,wPC);Wk2=RandScale*randn(wPC,nParc);


%% Make NADAM Variables
mA=zeros(size(A));         mW1=zeros(size(W1));    mD=zeros(size(D));     
mWk1=zeros(size(Wk1));  mWk2=mWk1';
%%%%
nA=zeros(size(A));      nW1=zeros(size(W1));    nD=zeros(size(D));
nWk1=zeros(size(Wk1));  nWk2=nWk1';
%%%%


if strcmpi(MINDy_type,'HRF')
    H1=repmat(H1start,nX,1);
    H2=repmat(H2start,nX,1);
    mH1=zeros(size(H1));nH1=mH1;
    mH2=zeros(size(H2));nH2=mH2;
end