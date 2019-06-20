function[Out]=MINDyInt(ooP,Start,DnSamp,dt,dW,tEnd,doPlot)
%% Simulate a fitted MINDy model

%% Assumes C=0
%% Always uses W2

%% HARD-CODED for speed

if ~isempty(ooP.Param{5})
    W=ooP.Param{5};
else
    W=ooP.Param{1};
end
A=ooP.Param{2};
D=ooP.Param{6};
B=ooP.Param{3};
A2=A.^2;

A2=A2./(B(:,1).^2);
B5P=(B(:,2)+.5)./B(:,1);
B5N=(B(:,2)-.5)./B(:,1);
%% Factor out the dt terms and add xt
DDT=1-(dt*D);
WDT=dt*W.*(B(:,1)');
if ~isempty(ooP.Param{4})
Cdt=ooP.Param{4}*dt;
else
    Cdt=0;
end
nX=size(Start,1);
nDat=size(Start,2);

tVec=0:(dt*DnSamp):tEnd;
nT=numel(tVec);


if nDat==1
    Out=nan(nX,nT);
    Out(:,1)=Start;
else
    Out=nan(nX,nT,nDat);
    Out(:,1,:)=Start;
end
if dW~=0
    if nDat==1
    Noise=dW*randn(nX,tEnd/dt);
    else
    Noise=dW*randn(nX,nDat,tEnd/dt);
    end
    
    
if nDat==1
    for i=1:(tEnd/dt)
        Start=WDT*(sqrt(A2+(Start+B5P).^2)-sqrt(A2+(Start+B5N).^2))+DDT.*Start+Noise(:,i);
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
        Start=WDT*(sqrt(A2+(Start+B5P).^2)-sqrt(A2+(Start+B5N).^2))+DDT.*Start+Noise(:,:,i)+Cdt;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
    end
end

else

if nDat==1
    for i=1:(tEnd/dt)
        Start=WDT*(sqrt(A2+(Start+B5P).^2)-sqrt(A2+(Start+B5N).^2))+DDT.*Start+Cdt;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp))=Start;
        end
    end
else
    for i=1:(tEnd/dt)
        Start=WDT*(sqrt(A2+(Start+B5P).^2)-sqrt(A2+(Start+B5N).^2))+DDT.*Start;
        if mod(i,DnSamp)==0
            Out(:,1+(i/DnSamp),:)=Start;
        end
    end
end
end

if ~isempty(doPlot)&&strcmpi(doPlot(1),'y')
    figure
    if nDat==1
        plot(tVec,Out);
    else
        plot(tVec,squeeze(Out(:,:,1)));
    end
end
end
