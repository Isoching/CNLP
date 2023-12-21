%% test Harverly Case1

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
over=sdpvar(1,1,'full');
under=sdpvar(1,1,'full');
XP=sdpvar(1,1,'full');
% p=sdpvar(1,1,'full');
% xp=sdpvar(1,1,'full');
% yp=sdpvar(1,1,'full');
profit=9*X+15*Y-6*A-16*B-10*C;
Const=[];
Const0=[];
%流量约束
Const=Const+([A,B,C,P,X,Y,PX,PY,CX,CY,over,under]>=0);
Const=Const+(X<=100);
Const=Const+(Y<=200);
%性质取值scale
% Const=Const+(1<=p<=3);
%Mass Balance
Const=Const+(A+B==P);
Const=Const+(P==PX+PY);
Const=Const+(C==CX+CY);
Const=Const+(X==PX+CX);
Const=Const+(Y==PY+CY);
%性质平衡
p=0.9;
f=0.1;
%Const=Const+(3*A+1*B==p*P+Rprop);
Const0=Const0+(p*PX+p*PY-3*A-B+over-under==0);
Const0=Const0+(p*PX+2*CX-2.5*X+f*(over-under)<=0);
Const0=Const0+(p*PY+2*CY-1.5*Y+(1-f)*(over-under)<=0);


ops=sdpsettings('verbose',2,'savesolveroutput',1,'savesolverinput',1,'saveyalmipmodel',1,'solver','cplex');
sol=optimize(Const+Const0,-profit,ops)
iter=1;
delf=nan;
Result=table(A,B,C,P,PX,PY,CX,CY,X,Y,p,profit,over,under,f,delf,iter);
for i=1:size(Result,2)
   V=table2cell(Result(1,i));
   Result(1,i)=num2cell(value(V{1,1}));
end
Result
delf=1;
while abs(value(over))+abs(value(under))>=10^-5||abs(delf)>=10^-4   %当上下界的绝对值之和大于允许误差
    f0=f;
    p0=p;
    if value(P)~=0
        f1=value(PX./P);    %迭代更新猜值
        f2=1-f1;
        p=p0+value(over./P)-value(under./P);
    else
        f1=0;
        f2=0;
        p=0;
    end
    Const1=[];
    Const1=Const1+(p*PX+p*PY-3*A-B+over-under==0);
    Const1=Const1+(p*PX+2*CX-2.5*X+f1*(over-under)<=0);
    Const1=Const1+(p*PY+2*CY-1.5*Y+f2*(over-under)<=0);
    iter=iter+1;
    sol=optimize(Const+Const1,-profit,ops)
    f=f1;
    delf=f-f0;    %和初始猜值的误差
    Result1=table(A,B,C,P,PX,PY,CX,CY,X,Y,p,profit,over,under,f,delf,iter);
    for i=1:size(Result1,2)
        V=table2cell(Result1(1,i));
        Result1(1,i)=num2cell(value(V{1,1}));
    end
    Result=[Result;Result1]
end