function[Out]=MINDy_BIG(X0,dX0,ParcMask,Pre,ParStr,doCell)  %#ok<INUSL>

MINDy_type='BIG'; %#ok<NASGU>
X=X0;
dX=dX0;
%% To overrule functional definition of mu
mu=1;
MINDy_Initialize

jRate=W1Rate;

ParcMask=sparse(MINDy_Process_ParcMask(ParcMask));
ParcMaskT=sparse(ParcMask');

%% count in parcel's place
ParcCount=sum(ParcMask,2);
ParcCM=ParcMaskT*ParcCount;


if mean(ParcCount==0)~=0
    error(strcat('No vertices assigned to parcel(s):',num2str(find(ParcCount==0))))
end

iRep=1; %#ok<NASGU>
doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0; %#ok<NASGU>

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

Aep=A; %#ok<NODEF>
Axi=Aep.^-.5;
iRep=1; %#ok<NASGU>
doRecord=(strcmpi(ParStr.RecCorr,'y')+strcmpi(ParStr.RecW,'y')+strcmpi(ParStr.RecDdiff,'y')+strcmpi(ParStr.RecA,'y'))~=0;

nX=size(X0,1);

C=zeros(nX,1);
nT=size(X0,2);
BatchInd=randi(nT,NBatch,BatchSz);

Decay=Dmin+D.^2; %#ok<NODEF,NASGU>

WK=Wk1*Wk2; %#ok<NODEF>
W0=W1+WK; %#ok<NODEF,NASGU>


Aep=Axi.^-2;

if strcmpi(doCell(1),'y')
%% Faster Indexing
X0=mat2cell(X0,nX,ones(1,nT));
dX0=mat2cell(dX0,nX,ones(1,nT));
end


    J1=abs(J1);J2=abs(J2); %#ok<NODEF>
    J1=J1./(ParcMaskT*(ParcMask*J1));
    J2=J2./(ParcMaskT*(ParcMask*J2));

    
for iBatch=1:NBatch
    
    if mod(iBatch,50)==1
        disp([iBatch NBatch])
        if doRecord
 %           ModularRecord_00
        end
    end
    A=Aep; %#ok<NASGU>
bInd=BatchInd(iBatch,:);

if strcmpi(doCell(1),'y')
X=[X0{bInd}];
dX=[dX0{bInd}];
else
    X=X0(:,bInd);dX=dX0(:,bInd);
end
Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;

Aep=Axi.^-2;

P1=sqrt(Aep+(X.*(X+B1(:,1))));  %#ok<NODEF>
P2=sqrt(Aep+(X.*(X-B1(:,1))));

PolyZhat=(J2./B1(:,1)).*(P1-P2);
%PolyZfull=(1./B1(:,1)).*(P1-P2);
PolyZ=ParcMask*PolyZhat;

parcW0=ParcCount.*W0;

PW=ParcMaskT*((parcW0)*PolyZ);
%PW=ParcMaskT*(W0)*PolyZ;


E=dX-(J1.*PW-Decay.*X);

gradD=-2*D.*mean(E.*X,2);

%E=E./ParcCM;

PJE=(ParcMask*(J1.*E));

%PreParc=ParcMaskT*(W0')*(ParcMask*(J1.*E));
PreParc=ParcMaskT*((parcW0')*PJE);

%PreParc=ParcMaskT*(W0'.*ParcCount')*((ParcMask*(J1.*E))./ParcCount);

EPW=mean(E.*PW,2);
%% Switched Polyfull to PolyZhat./J2
PPZ=mean(PreParc.*PolyZhat,2);

gradJ1=EPW-(ParcMaskT*((ParcMask*(J1.*EPW))));
gradJ2=(PPZ./J2)-(ParcMaskT*(ParcMask*(PPZ)));

%% Really fitting grad (A.^-.5)
gradA=((-Aep.^(1.5).*J2)./B1(:,1)).*mean(PreParc.*(1./P1-1./P2),2);
%gradW1=ParcCount.*((ParcMask*(J1.*E))./ParcCount)*(PolyZ)'/BatchSz;
%gradW1=(ParcMask*(J1.*E))*(PolyZ)'/BatchSz;
gradW1=(ParcCount.*(PJE)*(PolyZ)')/BatchSz;


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
gradW1=gradW1-wCost-diag(SpDiag*sign(diag(W1)));

%    if mod(iBatch,50)==1
%       Out.E(:,ceil(iBatch/50))=mean(E.^2,2);
%    end

mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;
mJ1=mu*mJ1+(1-mu)*gradJ1;
mJ2=mu*mJ2+(1-mu)*gradJ2;


nW1=v*nW1+(1-v)*gradW1.^2;
nD=v*nD+(1-v)*gradD.^2;
nA=v*nA+(1-v)*gradA.^2;
nWk1=v*nWk1+(1-v)*gradWk1.^2;
nWk2=v*nWk2+(1-v)*gradWk2.^2;
nJ1=v*nJ1+(1-v)*gradJ1.^2;
nJ2=v*nJ2+(1-v)*gradJ2.^2;

gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);


Rtv0=sqrt(v0);


J1=J1+(((jRate*mHat/Rtv0)*mJ1)+(jRate*gHat/Rtv0)*gradJ1)./(sqrt(nJ1)+(Reg/Rtv0));
J2=J2+(((jRate*mHat/Rtv0)*mJ2)+(jRate*gHat/Rtv0)*gradJ2)./(sqrt(nJ2)+(Reg/Rtv0));

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


    J1=abs(J1);J2=abs(J2);
    J1=J1./(ParcMaskT*(ParcMask*J1));
    J2=J2./(ParcMaskT*(ParcMask*J2));

end

if strcmpi(doCell(1),'y')
X=[X0{:}];dX=[dX0{:}];
else
    X=X0;dX=dX0;
end

Out.AutoCorr=VecDiagCorrT(X,dX);

Decay=Dmin+D.^2;
WK=Wk1*Wk2;
W0=W1+WK;

Aep=Axi.^-2;

P1=sqrt(Aep+(X.*(X+B1(:,1))));
P2=sqrt(Aep+(X.*(X-B1(:,1))));

PolyZ=ParcMask*((J2./B1(:,1)).*(P1-P2));    
%Out.Corr=VecDiagCorrT(dX,J1.*(ParcMaskT*(ParcCount.*W0)*PolyZ)-(Decay.*X));

%PW=ParcMaskT*(ParcCount.*W0)*PolyZ;

%Out.Corr=VecDiagCorrT(dX,J1.*(ParcMaskT*(W0)*PolyZ)-(Decay.*X));



A0=sqrt((Aep./(B1(:,1).^2)-.25));
B10(:,1)=B1(:,1).^-1;   B10(:,2)=B1(:,2)./B1(:,1);
A=A0;   B1=B10;

if doRecord
    gField=fields(Rec);
    for kk=1:numel(gField)
        Out.(gField{kk})=Rec.(gField{kk});
    end
end
Out.Param={W1,A,B1,C,W1+Wk1*Wk2,Decay};
Out.Mat={W1,Wk1,Wk2};
Out.ParcMask=ParcMask;
Out.jScale={J1,J2};
Out.ParcCount=ParcCount;
end



