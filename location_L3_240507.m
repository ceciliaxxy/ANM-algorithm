clear all;close
L = 3; 
deltaf=1e6;
Dp=4;
c=3e8;
Num_subcarr=64;
Num_pil=Num_subcarr/Dp;
index_Loc=1;
Num_loc=1;
Num_ante=64;
while index_Loc<=Num_loc
    index_Loc
    x(1,index_Loc)=3; 
    y(1,index_Loc)=4;
    x(2,index_Loc)=10; 
    y(2,index_Loc)=11;
    x(3,index_Loc)=x(2,index_Loc)+1/Num_pil/Dp/deltaf*c/8; 
    y(3,index_Loc)=y(2,index_Loc)+1/Num_pil/Dp/deltaf*c/8;
    Dis_TX2BD(:,index_Loc)=sqrt(x(:,index_Loc).^2+y(:,index_Loc).^2);
    Dis_BD2RX(:,index_Loc)=sqrt((x(1,index_Loc)-20).^2+y(:,index_Loc).^2);
    Dis_TX2BD2RX(:,index_Loc)=Dis_TX2BD(:,index_Loc)+Dis_BD2RX(:,index_Loc);
  
    sin_theta(:,index_Loc)=abs(y(:,index_Loc))./sqrt((x(:,index_Loc)-20).^2+y(:,index_Loc).^2);
    index_Loc=index_Loc+1;

end

Dis_TX2RX=sqrt(20^2);

save location_L3_240507.mat