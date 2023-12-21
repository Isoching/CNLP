%% test Harverly Case1
% 本测试算例采用非线性求解 求解器使用IPOPT
clear all;clc;

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
p=sdpvar(1,1,'full');
xp=sdpvar(1,1,'full');
yp=sdpvar(1,1,'full');
profit=9*X+15*Y-6*A-16*B-10*C;
Const=[];
%流量约束
Const=Const+([A,B,C,P,X,Y,PX,PY,CX,CY]>=0);
Const=Const+(X<=100);
Const=Const+(Y<=200);
%性质取值scale
Const=Const+(1<=p<=3);
%Mass Balance
Const=Const+(A+B==P);
Const=Const+(P==PX+PY);
Const=Const+(C==CX+CY);
Const=Const+(X==PX+CX);
Const=Const+(Y==PY+CY);
%性质平衡
Const=Const+(3*A+1*B==p*P);
Const=Const+(xp*X==p*PX+2*CX);
Const=Const+(yp*Y==p*PY+2*CY);
%性质约束
Const=Const+(xp<=2.5);
Const=Const+(yp<=1.5);

ops=sdpsettings('verbose',2,'savesolveroutput',1,'savesolverinput',1,'saveyalmipmodel',1,'solver','IPOPT');
sol=optimize(Const,-profit,ops)
Result=table(A,B,C,P,PX,PY,CX,CY,X,Y,p,xp,yp,profit);
for i=1:size(Result,2)
   V=table2cell(Result(1,i));
   Result(1,i)=num2cell(value(V{1,1}));
end
Result