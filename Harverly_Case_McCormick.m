%% Harverly_Case_McCormick
A=sdpvar(1,1,'full');
B=sdpvar(1,1,'full');
C=sdpvar(1,1,'full');
P=sdpvar(1,1,'full');
X=sdpvar(1,1,'full');
Y=sdpvar(1,1,'full');
PX=sdpvar(1,1,'full');
PY=sdpvar(1,1,'full');
CX=sdpvar(1,1,'full');
CY=sdpvar(1,1,'full');
Rprop=sdpvar(1,1,'full');
XP=sdpvar(1,1,'full');
prop=sdpvar(1,1,'full');
profit=9*X+15*Y-6*A-16*B-10*C;

Const=[];
%Á÷Á¿Ô¼Êø
Const=Const+([A,B,C,P,X,Y,PX,PY,CX,CY]>=0);
Const=Const+(X<=100);
Const=Const+(Y<=200);
%Mass Balance
Const=Const+(A+B==P);
Const=Const+(P==PX+PY);
Const=Const+(C==CX+CY);
Const=Const+(X==PX+CX);
Const=Const+(Y==PY+CY);
Const=Const+(XP==3*A+B);
propU=10;
propL=0;
flowU=1000;
flowL=0;
F=[];

%XP=P*prop
F=F+(XP>=flowL*prop+P*propL-flowL*propL);
F=F+(XP>=flowU*prop+P*propU-flowU*propU);
F=F+(XP<=flowU*prop+P*propL-flowU*propL);
F=F+(XP<=flowL*prop+P*propU-flowL*propU);

fy=sdpvar(1,1,'full');
fx=sdpvar(1,1,'full');
XP1=sdpvar(1,1,'full');
XP2=sdpvar(1,1,'full');
XPU=flowU*propU;
XPL=0;
fL=0;
fU=1;
% XP1=fy*XP
F=F+(XP1>=fL*XP+fy*XPL-fL*XPL);
F=F+(XP1>=fU*XP+fy*XPU-fU*XPU);
F=F+(XP1<=fU*XP+fy*XPL-fU*XPL);
F=F+(XP1<=fL*XP+fy*XPU-fL*XPU);
%XP2=fx*XP
F=F+(XP2>=fL*XP+fx*XPL-fL*XPL);
F=F+(XP2>=fU*XP+fx*XPU-fU*XPU);
F=F+(XP2<=fU*XP+fx*XPL-fU*XPL);
F=F+(XP2<=fL*XP+fx*XPU-fL*XPU);
F=F+(fy+fx==1);
%PY=fy*P
F=F+(PY>=fL*P+fy*flowL-fL*flowL);
F=F+(PY>=fU*P+fy*flowU-fU*flowU);
F=F+(PY<=fU*P+fy*flowL-fU*flowL);
F=F+(PY<=fL*P+fy*flowU-fL*flowU);
%PX=fx*P
F=F+(PX>=fL*P+fx*flowL-fL*flowL);
F=F+(PX>=fU*P+fx*flowU-fU*flowU);
F=F+(PX<=fU*P+fx*flowL-fU*flowL);
F=F+(PX<=fL*P+fx*flowU-fL*flowU);
%XP1=PY*prop
F=F+(XP1>=flowL*prop+PY*propL-flowL*propL);
F=F+(XP1>=flowU*prop+PY*propU-flowU*propU);
F=F+(XP1<=flowU*prop+PY*propL-flowU*propL);
F=F+(XP1<=flowL*prop+PY*propU-flowL*propU);
%XP2=PX*prop
F=F+(XP2>=flowL*prop+PX*propL-flowL*propL);
F=F+(XP2>=flowU*prop+PX*propU-flowU*propU);
F=F+(XP2<=flowU*prop+PX*propL-flowU*propL);
F=F+(XP2<=flowL*prop+PX*propU-flowL*propU);
F=F+(XP1+XP2==XP);

Const=Const+(XP1+2*CY-1.5*Y<=0);
Const=Const+(XP2+2*CX-2.5*X<=0);

ops=sdpsettings('verbose',2,'savesolveroutput',1,'savesolverinput',1,'saveyalmipmodel',1,'solver','cplex');
sol=optimize(Const+F,-profit,ops)
iter=1;
valueR=100;
vfy=value(fy);
vfx=value(fx);
vprop=value(prop);
Result=table(A,B,C,P,PX,PY,CX,CY,X,Y,profit,XP,XP1,XP2,vprop,vfy,vfx,iter,valueR);
for i=1:size(Result,2)
   V=table2cell(Result(1,i));
   Result(1,i)=num2cell(value(V{1,1}));
end
Result


while iter<=1||valueR>=0.001
    F_iter=[];
    F_iter=F_iter+(XP==value(XP./P)*P+Rprop);
    F_iter=F_iter+(XP1==value(XP./P)*PY+value(PY)./value(P)*Rprop);
    F_iter=F_iter+(XP2==value(XP./P)*PX+value(PX)./value(P)*Rprop);
    sol=optimize(Const+F_iter,-profit,ops)
    valueR=abs(value(Rprop));
    vprop=XP./P;
    vfy=PY./P;
    vfx=PX./P;
    iter=iter+1;
    Result1=table(A,B,C,P,PX,PY,CX,CY,X,Y,profit,XP,XP1,XP2,vprop,vfy,vfx,iter,valueR);
    for i=1:size(Result1,2)
        V=table2cell(Result1(1,i));
        Result1(1,i)=num2cell(value(V{1,1}));
    end
    Result=[Result;Result1]
    
end
Reuslt_McCormick=Result