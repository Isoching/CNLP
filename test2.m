clear;


%
f = [-45, -83];
A = [8, 7;
    7, 20;];
b = [56; 70];
Aeq = [];
beq = [];
lbnd = [0; 0];
ubnd = [7; 2];

c=7.1;
% x1*x2=c,
[x, fval,exitflag] = birelaxlinprog(f,A,b,Aeq,beq,lbnd, ubnd,c);

%
upperBound=inf;
lowerBound=-inf;
openlist={[x,lbnd,ubnd]};
sorthelplist={[fval,fval]};%[min,this]
count=1;
rst={};
drt=1e-4;
camgap=1e-4;%%用于比较
tic;
while exitflag&&(not(isempty(openlist)))
    
    t=sorthelplist{1};
    fv=t(2);
    m0=inf;
    if(count~=1)
        ca=sorthelplist{2};
        m0=ca(1);
    end
    %if not((lowerBound <= fv) && (fv < upperBound))
    if not(inrange(fv,lowerBound,upperBound,camgap))
        %         disp(1);
        openlist=openlist(2:end);
        sorthelplist=sorthelplist(2:end);
        lowerBound=max(lowerBound,m0);
        count=count-1;
        continue;
    end
    node0=openlist{1};%%取头部元素
    x0=node0(:,1);
    %     fprintf("%f,%f,%f\n", lowerBound,x0);
    lbnd0=node0(:,2);
    ubnd0=node0(:,3);
    %% 选择
    [r,k]=max(ubnd0-lbnd0);
    %     k=1;
    % disp(r);
    r=min(congap(x0,c),r);
    if(r<drt)%% 达到了r-optimal
        openlist=openlist(2:end);
        sorthelplist=sorthelplist(2:end);
        %         lowerBound=max(lowerBound,m0);
        count=count-1;
        upperBound=min(upperBound,fv);
        rst=[x0,rst(:)'];
        continue;
    end
    sp=0.5*(ubnd0(k)+lbnd0(k));
    lbnd1=lbnd0;
    ubnd1=ubnd0;ubnd1(k)=sp;
    lbnd2=lbnd0;lbnd2(k)=sp;
    ubnd2=ubnd0;
    %% 计算
    [x1, fval1,exitflag1] = birelaxlinprog(f,A,b,Aeq,beq,lbnd1, ubnd1,c);
    [x2, fval2,exitflag2] = birelaxlinprog(f,A,b,Aeq,beq,lbnd2, ubnd2,c);
    
    if(exitflag1&&exitflag2)
        node1=[x1,lbnd1, ubnd1];
        node2=[x2,lbnd2, ubnd2];
        if(fval1<=fval2)
            openlist={node1,node2,openlist{2:end}};
            m1=min(fval2,m0);
            m2=min(fval1,m1);
            sorthelplist={[m2,fval1],[m1,fval2],sorthelplist{2:end}};
        else
            openlist={node2,node1,openlist{2:end}};
            m1=min(fval1,m0);
            m2=min(fval2,m1);
            sorthelplist={[m2,fval2],[m1,fval1],sorthelplist{2:end}};
        end
        lowerBound=max(lowerBound,m2);
        %         disp(m2-lowerBound)
        count=count+1;
    elseif(exitflag1)
        node1=[x1,lbnd1, ubnd1];
        m1=min(fval1,m0);
        openlist={node1,openlist{2:end}};
        sorthelplist={[m1,fval1],sorthelplist{2:end}};
        lowerBound=max(lowerBound,m1);
        %         disp(m1-lowerBound)
        %         count=count+1;
    elseif(exitflag2)
        node2=[x2,lbnd2, ubnd2];
        m1=min(fval2,m0);
        openlist={node2,openlist{2:end}};
        sorthelplist={[m1,fval2],sorthelplist{2:end}};
        lowerBound=max(lowerBound,m1);
        %         count=count+1;
    else
        %%无解
        openlist=openlist(2:end);
        sorthelplist=sorthelplist(2:end);
        lowerBound=max(lowerBound,m0);
        count=count-1;
    end
end
toc;
if(not(isempty(rst)))
    xbest=rst{1};
    fprintf("x=[%1.16f,%1.16f],gap=%1.16f,f=%f\n",xbest,abs(xbest(1)*xbest(2)-c),dot(f,xbest));
else
    disp('无解');
end
% inline function
function [x, fval,exitflag] = birelaxlinprog(f,A,b,Aeq,beq,lbnd, ubnd,c)
A2=[lbnd(2),lbnd(1);
    ubnd(2),ubnd(1);
    -lbnd(2),-ubnd(1);
    -ubnd(2),-lbnd(1);];
b2=[lbnd(1)*lbnd(2)+c;ubnd(1)*ubnd(2)+c;-ubnd(1)*lbnd(2)-c;-ubnd(2)*lbnd(1)-c];
%%init;
[x, fval,exitflag0] = linprog(f, [A;A2], [b;b2], Aeq, beq, lbnd, ubnd);
exitflag=(exitflag0==1);
end

function flag=inrange(x,a,b,ep)
flag=(x+ep>=a)&&(x-ep<=b);
end

function gap=congap(x,c)
gap=abs(x(1)*x(2)-c);
end