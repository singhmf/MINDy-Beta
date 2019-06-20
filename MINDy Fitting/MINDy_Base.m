function[Out]=MINDy_Base(X0,dX0,Pre,ParStr)  %#ok<INUSL>
%% Axi denotes the XI transformation for A
%% Aep denotes the EPSILON transformation for A (EP=XI.^-2)

MINDy_type='Base'; %#ok<NASGU>
X=X0;
dX=dX0;
%% To overrule functional definition of mu
mu=1;
%% Preallocation, Extract variables from Pre, ParStr
MINDy_Initialize

Aep=A; %#ok<NODEF>
Axi=Aep.^-.5;
iRep=1; %#ok<NASGU>
doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0;

nX=size(X0,1);

C=zeros(nX,1);
if iscell(X)
    X0=[X{:}];
else
    X0=X;
end
if iscell(dX)
    dX0=[dX{:}];
else
    dX0=dX;
end

nT=size(X0,2);
BatchInd=randi(nT,NBatch,BatchSz);


Decay=Dmin+D.^2; %#ok<NODEF,NASGU>
WK=Wk1*Wk2; %#ok<NODEF>
W0=W1+WK; %#ok<NODEF,NASGU>
Aep=Axi.^-2;
P1=sqrt(Aep+(X.*(X+B1(:,1)))); %#ok<NODEF>
P2=sqrt(Aep+(X.*(X-B1(:,1))));
PolyZ=(B1(:,1).^-1).*(P1-P2);     %#ok<NASGU>


for iBatch=1:NBatch
    if mod(iBatch,50)==1
        disp([iBatch NBatch])
        if doRecord
            ModularRecord
        end
    end
    A=Aep; %#ok<NASGU>
bInd=BatchInd(iBatch,:);

X=X0(:,bInd);
dX=dX0(:,bInd);

%% Calculate useful variables
Decay=Dmin+D.^2;
WK=Wk1*Wk2;
W0=W1+WK;
%% Transfer Function
Aep=Axi.^-2;
P1=sqrt(Aep+(X.*(X+B1(:,1)))); 
P2=sqrt(Aep+(X.*(X-B1(:,1))));
PolyZ=(B1(:,1).^-1).*(P1-P2);    

%% Model Fit Error
E=dX-(W0*PolyZ-(Decay.*X));

%% Calculate Gradients

W1E=W0'*E;
gradD=-2*D.*mean(E.*X,2);
gradC=mean(E,2); %#ok<NASGU>
gradW1=(E*PolyZ')/BatchSz;

if L2SpPlsEN~=0
WcTemp=gradW1-L2SpPlsEN*WK;
gradWk1=(WcTemp)*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*WcTemp-SpPls2*sign(Wk2);
else
gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2);
end

%% Took away the normalization part
%mP=mean(abs(PolyZ),2)'; 
%wCost=((SpScale)*sign(W1)+(2*ENL2)*W1).*(mP/mean(mP));
if ENL2~=0
wCost=((SpScale)*sign(W1)+(2*ENL2)*W1);
else
    wCost=((SpScale)*sign(W1));
end
gradW1=gradW1-wCost-diag(SpDiag*diag(gradW1));
DpTemp=(B1(:,1).^-1).*W1E.*((1./P1)-(1./P2));
%% Really fitting grad (A.^-.5)
gradA=-Aep.^(1.5).*mean(DpTemp,2);

%% Record Error
    if mod(iBatch,50)==1
        Out.E(:,ceil(iBatch/50))=mean(E.^2,2);
    end

%% Moving Gradient
mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;

%% Moving Squared Gradient
nW1=v*nW1+(1-v)*gradW1.^2;
nD=v*nD+(1-v)*gradD.^2;
nA=v*nA+(1-v)*gradA.^2;
nWk1=v*nWk1+(1-v)*gradWk1.^2;
nWk2=v*nWk2+(1-v)*gradWk2.^2;

%% Moving NADAM parameters to remove initialization bias
gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);
Rtv0=sqrt(v0);

%% Update Values
Wk1=Wk1+(((wk1Rate*mHat/Rtv0)*mWk1)+(wk1Rate*gHat/Rtv0)*gradWk1)./(sqrt(nWk1)+(Reg/Rtv0));
Wk2=Wk2+(((wk2Rate*mHat/Rtv0)*mWk2)+(wk2Rate*gHat/Rtv0)*gradWk2)./(sqrt(nWk2)+(Reg/Rtv0));
W1=W1+(((W1Rate*mHat/Rtv0)*mW1)+(W1Rate*gHat/Rtv0)*gradW1)./(sqrt(nW1)+(Reg/Rtv0));

Axi=Axi+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0));

if DRegScale>0 
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end
Aep=Axi.^-2;
Aep=min(max(Aep,.25*(B1(:,1).^2)+AdiffMin),.25*(B1(:,1).^2)+AdiffMax);
Axi=Aep.^-.5;
end

X=X0;dX=dX0;

%% Evaluate Autocorrelation
Out.AutoCorr=DiagCorr(X',dX');

Decay=Dmin+D.^2;
WK=Wk1*Wk2;
W0=W1+WK;

Aep=Axi.^-2;

P1=sqrt(Aep+(X.*(X+B1(:,1))));
P2=sqrt(Aep+(X.*(X-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    
%% Evaluate Uncorrected goodness-of-fit
Out.Corr=DiagCorr(dX',(W0*PolyZ-(Decay.*X))');

%% Convert Variables
A0=sqrt((Aep./(B1(:,1).^2)-.25));
B10(:,1)=B1(:,1).^-1;   B10(:,2)=B1(:,2)./B1(:,1);
A=A0;   B1=B10;

%% Store Estimate Timeseries (if requested)
if doRecord
    gField=fields(Rec);
    for kk=1:numel(gField)
        Out.(gField{kk})=Rec.(gField{kk});
    end
end
%% Store Variables
Out.Param={W1,A,B1,C,W1+Wk1*Wk2,Decay};
Out.Mat={W1,Wk1,Wk2};
end



