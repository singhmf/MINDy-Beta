function[Out]=MINDy_Base_Rec(X,Pre,ParStr)
%% Also Requires dt and nStep


%function[Out]=MINDy_HRF(xHRFpoly,DxHRFpoly,PolyVec,ParStr)
%% Fit MINDy HRF with polynomials for both derivative and original time series
%% {'p00'} {'p10'} {'p01'} {'p20'} {'p11'} {'p02'} {'p30'} {'p21'} {'p12'} {'p03'}
MINDy_type='Rec'; %#ok<NASGU>


MINDy_Initialize

nTx=size(X,2)-1;
nT=nStep/dt;

nX=size(X,1);
C=zeros(nX,1);

BatchInd=randi(nTx-(1+nStep),NBatch,BatchSz);

for iBatch=1:NBatch
    
    if mod(iBatch,50)==1
        disp([iBatch NBatch])        
    end
    
%%%%%%%%%%%%%%%%
%% Forward Pass
%%%%%%%%%%%%%%%%

bInd=BatchInd(iBatch,:);

%% Note: Could compute unique elements so don't need to recompute when there's overlap

%size(X)
%max(bInd)
for iRec=1:(nStep+1)
XC{iRec}=X(:,bInd+iRec-1);
end

Decay=Dmin+D.^2;

WK=Wk1*Wk2;
W0=W1+WK;

X0=XC{1};    

%% Definitions of dtW0 and dtDecay
dtW0=dt*W0;
dtDecay=dt*Decay;

for iT=1:nT
bbTemp=X0+B1(:,2);
P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));
dZdA{iT}=(B1(:,1).^-1).*((1./P1)-(1./P2));
dZdX{iT}=dZdA{iT}.*(bbTemp+(B1(:,1)/2))+(1./P2);
PolyZ=(B1(:,1).^-1).*(P1-P2);    
PolyZC{iT}=PolyZ;
%%% XP=starting X at each time step. XP{1}=X0;
XP{iT}=X0;
X0=dtW0*PolyZ+(1-dtDecay).*X0;
if mod(iT,1/dt)==0
EC{iT*dt}=XC{iT*dt+1}-X0;
end
end

%%%%%%%%%%%%%%
%% GRADIENTS%%
%%%%%%%%%%%%%%

E=EC{end};
W1E=W0'*E;


DpTemp=W1E.*(dZdA{end});
gradW1=(E*PolyZ')/BatchSz;
gradW1Mat(:,:,1)=gradW1;
gradAMat(:,1)=(B1(:,1).^-1).*mean(DpTemp,2);
gradDMat(:,1)=-2*dt*D.*mean(E.*XP{end},2);

%% REMARK: PROBLEM: JACOBIANS dXt/dXt-1 are sample specific!
%% dEdXt is actually product of prod[dXt/dXt-1]'*E
dEdXt=E;%eye(nX);

for i=2:nT

tpIdx=1+nT-i;
    
    dEdXt=(1-dtDecay).*dEdXt+dZdX{tpIdx}.*(dtW0'*dEdXt);
if mod(tpIdx,1/dt)==0
    dEdXt=dEdXt+EC{tpIdx*dt};
end
gradAMat(:,i)=-mean((dtW0'*dEdXt).*dZdA{tpIdx},2)/2;
gradW1Mat(:,:,i)=(dEdXt*PolyZC{tpIdx}')*(dt/BatchSz);
%gradW1Mat(:,:,i)=gradW1;
gradDMat(:,i)=-(2*dt*D).*mean(dEdXt.*XP{tpIdx},2);
%% Remark: Can save gradWk1Mat, gradWk2mat for end
%gradWk1Mat(:,:,i)=gradW1*Wk2';
%gradWk2Mat(:,:,i)=Wk1'*gradW1;
end

gradD=mean(gradDMat,2);
gradW1=dt*squeeze(mean(gradW1Mat,3));
%% Saved multiplying until end
gradA=(A.^1.5).*mean(gradAMat,2);

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
gradW1=gradW1-wCost;


mW1=mu*mW1+(1-mu)*gradW1;
mD=mu*mD+(1-mu)*gradD;
mA=mu*mA+(1-mu)*gradA;
mWk1=mu*mWk1+(1-mu)*gradWk1;
mWk2=mu*mWk2+(1-mu)*gradWk2;

nW1=v*nW1+(1-v)*gradW1.^2;
nD=v*nD+(1-v)*gradD.^2;
nA=v*nA+(1-v)*gradA.^2;
nWk1=v*nWk1+(1-v)*gradWk1.^2;
nWk2=v*nWk2+(1-v)*gradWk2.^2;


gHat=(1-mu)/(1-mu^iBatch);
mHat=mu/(1-mu^(iBatch+1));
v0=1/(1-v^iBatch);


Rtv0=sqrt(v0);

Wk1=Wk1+(((wk1Rate*mHat/Rtv0)*mWk1)+(wk1Rate*gHat/Rtv0)*gradWk1)./(sqrt(nWk1)+(Reg/Rtv0));
Wk2=Wk2+(((wk2Rate*mHat/Rtv0)*mWk2)+(wk2Rate*gHat/Rtv0)*gradWk2)./(sqrt(nWk2)+(Reg/Rtv0));
W1=W1+(((W1Rate*mHat/Rtv0)*mW1)+(W1Rate*gHat/Rtv0)*gradW1)./(sqrt(nW1)+(Reg/Rtv0));
A=((A.^-.5)+(((ARate*mHat/Rtv0)*mA)+(ARate*gHat/Rtv0)*gradA)./(sqrt(nA)+(Reg/Rtv0))).^-2;


if DRegScale>0
D=abs(D+DRegScale*DRate*((mHat*mD)+gHat*gradD)./(sqrt(v0*nD)+DRegScale));
else    
D=abs(D+DRate*((mHat*mD)+gHat*gradD));
end
W1(~Mask)=0;

A=min(max(A,.25*(B1(:,1).^2)+AdiffMin),.25*(B1(:,1).^2)+AdiffMax);
if isfield(ParStr,'Dmax')
D=min(D,sqrt(Dmax)-Dmin);
end

end

%X=sum(HRFcomb.*xHRFpoly,3);
%dX=sum(HRFcomb.*DxHRFpoly,3);

%D=D.*std(dX,[],2);
%W1=W1.*std(dX,[],2);
%Wk1=Wk1.*std(dX,[],2);

X=X0(:,1:end-1);
dX=convn(X0,[1 -1],'valid');
Out.AutoCorr=DiagCorr(X',dX');


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



