function[Out]=MINDy_NeuralMass(X0,dX0,Pre,ParStr)  %#ok<INUSL>

%% Fits neural mass models of the form  dx/dt=-Dx+A(1-x)phi(Wx+c)
%%                                     phi(y)=y/(1-exp(-y))

if ~isfield(ParStr,'CRate')
    ParStr.CRate=ParStr.Rate;
end

MINDy_type='Base'; %#ok<NASGU>
X=X0;
dX=dX0;
%% To overrule functional definition of mu
mu=1;
MINDy_Initialize


iRep=1; %#ok<NASGU>
doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0;

nX=size(X0,1);
%nX=size(xHRFpoly,1);
C=zeros(nX,1);

nC=zeros(size(C));mC=zeros(size(C));

%H2=repmat(1.1,nX,1);
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
%nT=size(xHRFpoly,2);
BatchInd=randi(nT,NBatch,BatchSz);

Decay=Dmin+D.^2; %#ok<NODEF,NASGU>

WK=Wk1*Wk2; %#ok<NODEF>
W0=W1+WK; %#ok<NASGU,NODEF>



for iBatch=1:NBatch
    
    if mod(iBatch,50)==1
        disp([iBatch NBatch])
        if doRecord
            ModularRecord
        end
    end
bInd=BatchInd(iBatch,:);

X=X0(:,bInd);
dX=dX0(:,bInd);

Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;

Y=W0*X+C;
Zeta=1./(1-exp(-Y));
Phi=Y.*Zeta;

NonLin=(1-X).*Phi;

E=dX-(A.*NonLin-(Decay.*X)); 


gradD=-2*D.*mean(E.*X,2);

gradA=mean(E.*NonLin,2);

dPhidY=Zeta+Phi.*(1-Zeta);

tmpGrad=E.*(1-X).*dPhidY;

gradC=A.*mean(tmpGrad,2);

gradW1=(A/BatchSz).*(tmpGrad*X');

if L2SpPlsEN~=0
WcTemp=gradW1-L2SpPlsEN*WK;
gradWk1=(WcTemp)*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*WcTemp-SpPls2*sign(Wk2);
else
gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2);
end

wCost=((SpScale)*sign(W1)+(2*ENL2)*W1);
gradW1=gradW1-wCost-diag(SpDiag*sign(diag(W1)));

    if mod(iBatch,50)==1
        Out.E(:,ceil(iBatch/50))=mean(E.^2,2);
    end


mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;
mC=mu*mC+(1-mu)*gradC;

nW1=v*nW1+(1-v)*gradW1.^2;
nD=v*nD+(1-v)*gradD.^2;
nA=v*nA+(1-v)*gradA.^2;
nWk1=v*nWk1+(1-v)*gradWk1.^2;
nWk2=v*nWk2+(1-v)*gradWk2.^2;
nC=v*nC+(1-v)*gradC.^2;

gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);


Rtv0=sqrt(v0);

C=C+(((CRate*mHat/Rtv0)*mC)+(CRate*gHat/Rtv0)*gradC)./(sqrt(nC)+(Reg/Rtv0));

Wk1=Wk1+(((wk1Rate*mHat/Rtv0)*mWk1)+(wk1Rate*gHat/Rtv0)*gradWk1)./(sqrt(nWk1)+(Reg/Rtv0));
Wk2=Wk2+(((wk2Rate*mHat/Rtv0)*mWk2)+(wk2Rate*gHat/Rtv0)*gradWk2)./(sqrt(nWk2)+(Reg/Rtv0));
W1=W1+(((W1Rate*mHat/Rtv0)*mW1)+(W1Rate*gHat/Rtv0)*gradW1)./(sqrt(nW1)+(Reg/Rtv0));

A=A+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0));
%A=((A.^-.5)+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0))).^-2;

if DRegScale>0 
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end
end


Out.AutoCorr=DiagCorr(X',dX');

Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;



Y=W0*X+C;
Zeta=1./(1-exp(-Y));
Phi=Y.*Zeta;

NonLin=(1-X).*Phi;

Out.Corr=DiagCorr(dX',(A.*NonLin-(Decay.*X))');

if doRecord
    gField=fields(Rec);
    for kk=1:numel(gField)
        Out.(gField{kk})=Rec.(gField{kk});
    end
end
MINDy_StoreVal;
Out.Param={Val.W1,Val.A,Val.B1,Val.C,Val.W1+Val.Wk1*Val.Wk2,Val.D};
end



