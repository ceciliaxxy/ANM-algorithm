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
    x_test(1)=3; 
    y_test(1)=4;
    x_test(2)=10; 
    y_test(2)=11;
    x_test(3)=x_test(2)+1/Num_pil/Dp/deltaf*c/8; 
    y_test(3)=y_test(2)+1/Num_pil/Dp/deltaf*c/8;
    Dis_TX2BD_test=sqrt(x_test.^2+y_test.^2);
    Dis_BD2RX_test=sqrt((x_test-20).^2+y_test.^2);
    Dis_TX2BD2RX_test=Dis_TX2BD_test+Dis_BD2RX_test;
  
    sin_theta_test=abs(y_test)./sqrt((x_test-20).^2+y_test.^2);

    x(:,index_Loc)=x_test;
    y(:,index_Loc)=y_test;
    Dis_TX2BD(:,index_Loc)=Dis_TX2BD_test;
    Dis_BD2RX(:,index_Loc)=Dis_BD2RX_test;
    sin_theta(:,index_Loc)=sin_theta_test;
    index_Loc=index_Loc+1;

end

Dis_TX2RX=sqrt(20^2);

save location_L3_240507.mat