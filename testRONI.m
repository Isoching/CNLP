% %% testRONI
% for i=1:size(DATA,2)
%     DATA(i).comp_vol=DATA(i).comp_weight./DATA(i).comp_SPG.*1000;
%     DATA(i).comp_vol(isnan(DATA(i).comp_vol))=0;
%     DATA(i).comp_RONI=RON2RONI(DATA(i).comp_RON./100);
%     DATA(i).predict_RON=sum(DATA(i).comp_vol.*DATA(i).comp_RONI)./sum(DATA(i).comp_vol);
%     DATA(i).error=DATA(i).predict_RON-DATA(i).real_RON;
% end
% 
% 
% 
% %% editRONI
% data=DATA;
% abserr=0;
% data(1).edit_RONI=[]
% F=[];
% n=4;
% a=sdpvar(1,n,'full');
% for i=1:size(data,2)
%     data(i).edit_RONI=0;
%     for ii=1:n
%         data(i).edit_RONI=data(i).edit_RONI+a(ii).*data(i).comp_RON.^(ii-1);
%     end
%     data(i).predict_RON=sum(data(i).comp_vol.*data(i).edit_RONI)./sum(data(i).comp_vol);
%     data(i).error=data(i).predict_RON-data(i).real_RON;
%     abserr=abserr+abs(data(i).error);
% end
% ops=sdpsettings('solver','cplex');
% optimize(F,abserr,ops)

%% newtest
F=[];
VOL=W./density.*1000;
VOL(isnan(VOL))=0;
N=5;
a=sdpvar(1,N,'full');
RONI=sdpvar(size(RON,1),size(RON,2),'full');
% F=F+(a(4)==651);
% F=F+(a(3)==-1552.9);
% F=F+(a(2)==1272);
% F=F+(a(1)==-299.5);
for row=1:size(RON,1)
    for col=1:size(RON,2)
        temp=0;
        for n=1:N
            temp=temp+a(n).*RON(row,col).^(n-1);
        end
        F=F+(RONI(row,col)==temp);
    end
end
clear row col n temp
predict=sum(RONI.*VOL,2)./sum(VOL,2);
abserr=sum(abs(predict-result));
ops=sdpsettings('solver','cplex');
optimize(F,abserr,ops)
[value(abserr)]
       
    