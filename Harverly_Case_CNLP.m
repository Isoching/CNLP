 %% Harverly_Case_CNLP
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
XP1=sdpvar(1,1,'full');
XP2=sdpvar(1,1,'full');
prop=sdpvar(1,1,'full');
profit=9*X+15*Y-6*A-16*B-10*C;
Const=[];
%流量约束
%Const=[Const;([A,B,C,P,X,Y,PX,PY,CX,CY]>=0)];
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
Const=Const+(XP==XP1+XP2);        %XP是 AB混合中 混流的 S含量

Const=Const+(XP1+2*CY-1.5*Y<=0);   %自动得出一个X1 和X2 的初值 得到XP的分配情况
Const=Const+(XP2+2*CX-2.5*X<=0);

ops=sdpsettings('verbose',2,'savesolveroutput',1,'savesolverinput',1,'saveyalmipmodel',1,'solver','cplex');
sol=optimize(Const,-profit,ops)

iter=0;
valueR=100;
prop=nan;
Result=table(A,B,C,P,PX,PY,CX,CY,X,Y,profit,XP,XP1,XP2,prop,iter,valueR);
for i=1:size(Result,2)
   V=table2cell(Result(1,i));
   Result(1,i)=num2cell(value(V{1,1}));
end
Result

while valueR>=10^-3
    F=[];
    F=F+(XP==value(XP./P)*P+Rprop);     % value(XP./P)为上一次计算值，带入更新XP的值
    F=F+(XP1==value(XP./P)*PY+value(PY./P)*Rprop);
    F=F+(XP2==value(XP./P)*PX+value(PX./P)*Rprop);
    iter=iter+1;
    prop=value(XP./P);
    sol=optimize(Const+F,-profit,ops)
    valueR=value(Rprop);
    Result1=table(A,B,C,P,PX,PY,CX,CY,X,Y,profit,XP,XP1,XP2,prop,iter,valueR);
    for i=1:size(Result1,2)
        V=table2cell(Result1(1,i));
        Result1(1,i)=num2cell(value(V{1,1}));
    end
    Result=[Result;Result1]
end
ResultCNLP=Result