x=sdpvar(2,1,'full');
f = [-45, -83];
A = [8, 7;7, 20];
b = [56; 70];
Aeq = [];
beq = [];
lbnd = [0; 0];
ubnd = [7; 2];
c=7.1;
F=[];
F=F+(x>=lbnd );
F=F+(x<=ubnd );
F=F+(A*x<=b);
F=F+(x(1)*x(2)==c);
obj=f*x;
ops=sdpsettings('solver','ipopt');
sol=optimize(F,-obj,ops)