%% test McCormick
x1=sdpvar(1,1,'full');
x2=sdpvar(1,1,'full');
x3=sdpvar(1,1,'full');
x4=sdpvar(1,1,'full');
x5=sdpvar(1,1,'full');
x6=sdpvar(1,1,'full');
x7=sdpvar(1,1,'full');
P1=sdpvar(1,1,'full');
P2=sdpvar(1,1,'full');
x11=sdpvar(1,1,'full');
x12=sdpvar(1,1,'full');
x41=sdpvar(1,1,'full');
x42=sdpvar(1,1,'full');
x51=sdpvar(1,1,'full');
x52=sdpvar(1,1,'full');
x71=sdpvar(1,1,'full');
x72=sdpvar(1,1,'full');
XR5=sdpvar(1,1,'full');
XR7=sdpvar(1,1,'full');
ratio5=sdpvar(1,2);
ratio7=sdpvar(1,2);

F0=[];
F1=[];
F0=F0+(0<=x1<=10);
F0=F0+(0<=x2<=30);
F0=F0+(0<=x3<=70);
F0=F0+(0<=x4<=20);
F0=F0+(0<=x6<=100);

F0=F0+(x11>=1);
F0=F0+(x12>=1);
F0=F0+(x41>=1);
F0=F0+(x42>=1);
F0=F0+(x51>=1);
F0=F0+(x52>=1);
F0=F0+(x71>=1);
F0=F0+(x72>=1);

F0=F0+(x5==x2+x3);
F0=F0+(x7==x1+x6);
F0=F0+(P1==x11+x41+x51+x71);
F0=F0+(P2==x12+x42+x52+x72);
F0=F0+(x1==x11+x12);
F0=F0+(x4==x41+x42);
F0=F0+(x5==x51+x52);
F0=F0+(x7==x71+x72);

F1=F1+(XR5==56*x2+78*x3);
F1=F1+(XR7==45*x2+50*x3);
F1=F1+(45*x11+ratio5(1,1)*XR5+98*x41+ratio7(1,1)*XR7>=70*P1);
F1=F1+(45*x12+ratio5(1,2)*XR5+98*x42+ratio7(1,2)*XR7>=85*P2);

OBJ=-1500*x1-1600*x2-2000*x3-2300*x4-1400*x6+2000*P1+2200*P2;
ops=sdpsettings('verbose',2,'savesolveroutput',1,'savesolverinput',1,'saveyalmipmodel',1,'solver','cplex');
sol=optimize(F0,-OBJ,ops)


Rtio50=nan(1,2);
Rtio50(1,1)=value(x51./x5);
Rtio50(1,2)=value(x52./x5);
Rtio70=nan(1,2);
Rtio70(1,1)=value(x71./x7);
Rtio70(1,2)=value(x72./x7);

delratio=1;
delOBJ=100;

record=[1,Rtio50,Rtio70,delratio,value(OBJ),delOBJ];
while delratio>=10^-5||abs(delOBJ)>=10^-2
    OBJ0=value(OBJ);
    Rtio5=nan(1,2);
    Rtio5(1,1)=value(x51./x5);
    Rtio5(1,2)=value(x52./x5);
    Rtio7=nan(1,2);
    Rtio7(1,1)=value(x71./x7);
    Rtio7(1,2)=value(x72./x7);
    
    F2=F1;
    F2=replace(F2,ratio5,Rtio5);
    F2=replace(F2,ratio7,Rtio7);
    
    sol=optimize(F0+F2,-OBJ,ops);
    
    Rtio50=Rtio5;
    Rtio70=Rtio7;
    Rtio51=nan(1,2);
    Rtio51(1,1)=value(x51./x5);
    Rtio51(1,2)=value(x52./x5);
    Rtio71=nan(1,2);
    Rtio71(1,1)=value(x71./x7);
    Rtio71(1,2)=value(x72./x7);
    if all(isnan(Rtio51))
        Rtio51=ones(1,2).*0.5;
    end
    if all(isnan(Rtio71))
        Rtio71=ones(1,2).*0.5;
    end
    
    delratio=sum(abs(Rtio51-Rtio50))+sum(abs(Rtio71-Rtio70))
    OBJ1=value(OBJ)
    
    delOBJ=OBJ1-OBJ0;
    record=[record;[record(end,1)+1,Rtio51,Rtio71,delratio,value(OBJ),OBJ1-OBJ0]];
end
