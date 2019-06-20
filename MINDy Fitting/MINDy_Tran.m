function [PolyZ] = MINDy_Tran(X,A,B1)
%% Calculate the MINDy_Tran alone (not derivatives)
bbTemp=X+B1(:,2);
P1=sqrt(A+(bbTemp.*(bbTemp+B1(:,1))));
P2=sqrt(A+(bbTemp.*(bbTemp-B1(:,1))));
PolyZ=(B1(:,1).^-1).*(P1-P2);
end

