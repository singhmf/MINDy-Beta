function [PolyZ,dZdA,dZdX] = MINDy_Tran_Deriv(X,A,B1)
%% Calculate the MINDy_Tran and derivates wrt. A,X
bbTemp=X+B1(:,2);
P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));
dZdA=(B1(:,1).^-1).*((1./P1)-(1./P2));
dZdX=dZdA.*(bbTemp+(B1(:,1)/2))+(1./P2);
PolyZ=(B1(:,1).^-1).*(P1-P2);
end

