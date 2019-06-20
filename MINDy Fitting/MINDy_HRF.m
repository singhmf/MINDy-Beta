function[Out]=MINDy_HRF(X,Pre,ParStr)

%function[Out]=MINDy_HRF(xHRFpoly,DxHRFpoly,PolyVec,ParStr)
%% Fit MINDy HRF with polynomials for both derivative and original time series
%% {'p00'} {'p10'} {'p01'} {'p20'} {'p11'} {'p02'} {'p30'} {'p21'} {'p12'} {'p03'}
MINDy_type='HRF'; %#ok<NASGU>

mu=1;
if ~isstruct(X)
HRFout=MINDy_MakeHRFpoly(X,Pre,ParStr,'n');
else
    HRFout=X;
end
xHRFpoly=HRFout.Poly;
DxHRFpoly=HRFout.dPoly;
PolyVec=HRFout.PolyMat;

MINDy_Initialize

nX=size(xHRFpoly,1);
C=zeros(nX,1);

nT=size(xHRFpoly,2);
BatchInd=randi(nT,NBatch,BatchSz);

%H2=repmat(1.1,nX,1);


for iBatch=1:NBatch
    
    if mod(iBatch,50)==1
        disp([iBatch NBatch])
        Out.RecH1(:,ceil(iBatch/50))=H1;
        Out.RecH2(:,ceil(iBatch/50))=H2;
    end
    
HRFcomb=permute((H1.^PolyVec(1,:)).*(H2.^PolyVec(2,:)),[1 3 2]);

bInd=BatchInd(iBatch,:);

tX=xHRFpoly(:,bInd,:);
tDX=DxHRFpoly(:,bInd,:);

X=sum(HRFcomb.*tX,3);
dX=sum(HRFcomb.*tDX,3);
%dX=dX./std(dX,[],2);


Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;
bbTemp=X+B1(:,2);

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    

E=dX-(W0*PolyZ-(Decay.*X));

W1E=W0'*E;

gradD=-2*D.*mean(E.*X,2);
gradC=mean(E,2);
gradW1=(E*PolyZ')/BatchSz;

if L2SpPlsEN~=0
%gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1)-L2SpPlsEN*WK*Wk2'/sqrt(sum(WK(:).^2));
%gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2)-L2SpPlsEN*Wk1'*WK/sqrt(sum(WK(:).^2));
WcTemp=gradW1-L2SpPlsEN*WK;
gradWk1=(WcTemp)*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*WcTemp-SpPls2*sign(Wk2);
else
gradWk1=gradW1*Wk2'-SpPls1*sign(Wk1);
gradWk2=Wk1'*gradW1-SpPls2*sign(Wk2);
end

mP=mean(abs(PolyZ),2)'; %J=mean(abs(W1),1);%./(mP.^2);
wCost=((SpScale)*sign(W1)+(2*ENL2)*W1).*(mP/mean(mP));
gradW1=gradW1-wCost-diag(SpDiag*diag(gradW1));


DpTemp=(B1(:,1).^-1).*W1E.*((1./P1)-(1./P2));
dZdX=(bbTemp+(B1(:,1))/2).*DpTemp+(W1E./P2);

%% Really fitting grad (A.^-.5)
gradA=-A.^(1.5).*mean(DpTemp,2)/2;

dH1=permute(PolyVec(1,:).*(H1.^(PolyVec(1,:)-1)).*(H2.^PolyVec(2,:)),[1 3 2]);
dH2=permute(PolyVec(2,:).*(H1.^(PolyVec(1,:))).*(H2.^(PolyVec(2,:)-1)),[1 3 2]);

H1X=sum(dH1.*tX,3);
H1DX=sum(dH1.*tDX,3);
H2X=sum(dH2.*tX,3);
H2DX=sum(dH2.*tDX,3);

gradH1=mean(-E.*(H1DX+Decay.*H1X),2)+(mean(dZdX.*H1X,2));
gradH2=mean(-E.*(H2DX+Decay.*H2X),2)+(mean(dZdX.*H2X,2));

    if mod(iBatch,50)==1
        Out.E(:,ceil(iBatch/50))=mean(E.^2,2);
        Out.gradH1(:,ceil(iBatch/50))=gradH1;
        Out.gradH2(:,ceil(iBatch/50))=gradH2;
        
       delta=.01;
       ddH1=H1;ddH1(1)=ddH1(1)+delta;
       ddH2=H2;
       
HRFcomb=permute((ddH1.^PolyVec(1,:)).*(ddH2.^PolyVec(2,:)),[1 3 2]);
X=sum(HRFcomb.*tX,3);
dX=sum(HRFcomb.*tDX,3);
bbTemp=X+B1(:,2);
P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));
PolyZ=(B1(:,1).^-1).*(P1-P2);    
dE=dX-(W0*PolyZ-(Decay.*X));
Out.Eh1(ceil(iBatch/50))=(sum(sum(dE.^2))-sum(sum(E.^2)))/delta;

       ddH1=H1;
       ddH2=H2;ddH2(1)=ddH2(1)+delta;
       
HRFcomb=permute((ddH1.^PolyVec(1,:)).*(ddH2.^PolyVec(2,:)),[1 3 2]);
X=sum(HRFcomb.*tX,3);
dX=sum(HRFcomb.*tDX,3);
bbTemp=X+B1(:,2);
P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));
PolyZ=(B1(:,1).^-1).*(P1-P2);    
dE=dX-(W0*PolyZ-(Decay.*X));
Out.Eh2(ceil(iBatch/50))=(sum(sum(dE.^2))-sum(sum(E.^2)))/delta;

    end



%gradH2=gradH2-.05*(H2-1);
%gradH1=gradH1-.01*(H1-6);

%gradH1=mean(-E.*(H1DX),2)+(.5*(B1(:,1).^-1).*mean(DpTemp.*H1X,2));
%gradH2=mean(-E.*(H2DX),2)+(.5*(B1(:,1).^-1).*mean(DpTemp.*H2X,2));

mH1=mu*mH1+(1-mu)*gradH1;
mH2=mu*mH2+(1-mu)*gradH2;
mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;

nH1=v*nH1+(1-v)*gradH1.^2;
nH2=v*nH2+(1-v)*gradH2.^2;
nW1=v*nW1+(1-v)*gradW1.^2;
nD=v*nD+(1-v)*gradD.^2;
nA=v*nA+(1-v)*gradA.^2;
nWk1=v*nWk1+(1-v)*gradWk1.^2;
nWk2=v*nWk2+(1-v)*gradWk2.^2;


gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);


Rtv0=sqrt(v0);

%if iBatch<1000
Wk1=Wk1+(((wk1Rate*mHat/Rtv0)*mWk1)+(wk1Rate*gHat/Rtv0)*gradWk1)./(sqrt(nWk1)+(Reg/Rtv0));
Wk2=Wk2+(((wk2Rate*mHat/Rtv0)*mWk2)+(wk2Rate*gHat/Rtv0)*gradWk2)./(sqrt(nWk2)+(Reg/Rtv0));
W1=W1+(((W1Rate*mHat/Rtv0)*mW1)+(W1Rate*gHat/Rtv0)*gradW1)./(sqrt(nW1)+(Reg/Rtv0));

A=((A.^-.5)+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0))).^-2;

if DRegScale>0 
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end
%end

%if iBatch>1000
H1=H1+(((H1Rate*mHat/Rtv0)*mH1)+(H1Rate*gHat/Rtv0)*gradH1)./(sqrt(nH1)+(Reg/Rtv0));
H2=H2+(((H2Rate*mHat/Rtv0)*mH2)+(H2Rate*gHat/Rtv0)*gradH2)./(sqrt(nH2)+(Reg/Rtv0));
%end

%H1=(H1+H1Rate*((mHat*mH1)+gHat*gradH1));
%H2=(H2+H2Rate*((mHat*mH2)+gHat*gradH2));

W1(~Mask)=0;

A=min(max(A,.25*(B1(:,1).^2)+AdiffMin),.25*(B1(:,1).^2)+AdiffMax);

H1=max(H1min,min(H1max,H1));
H2=max(H2min,min(H2max,H2));
end

HRFcomb=permute((H1.^PolyVec(1,:)).*(H2.^PolyVec(2,:)),[1 3 2]);
X=sum(HRFcomb.*xHRFpoly,3);
dX=sum(HRFcomb.*DxHRFpoly,3);

%D=D.*std(dX,[],2);
%W1=W1.*std(dX,[],2);
%Wk1=Wk1.*std(dX,[],2);


Out.AutoCorr=DiagCorr(X',dX');
Out.HRF={H1,H2};

Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;
bbTemp=X+B1(:,2);

P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));

PolyZ=(B1(:,1).^-1).*(P1-P2);    
Out.Corr=DiagCorr(dX',(W0*PolyZ-(Decay.*X))');


A0=sqrt((A./(B1(:,1).^2)-.25));
B10(:,1)=B1(:,1).^-1;   B10(:,2)=B1(:,2)./B1(:,1);
A=A0;   B1=B10;

MINDy_StoreVal;
Out.Param={Val.W1,Val.A,Val.B1,Val.C,Val.W1+Val.Wk1*Val.Wk2,Val.D};
end



