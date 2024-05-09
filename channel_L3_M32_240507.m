clear all
%clc
load location_L3_240507.mat

Num_ante=32; %%receive antenna at the RX
Num_Loc=1;
N_symbo=64;
Num_symb=64;
N_fft=64;           
Space_pilot=4;       
N_pilot=N_fft/Space_pilot;    
N_sc=N_fft-N_pilot; 
M=4; 
 
 for index_Loc=1:Num_Loc

    PL_TX2BD((1),index_Loc)=1e-3*(Dis_TX2BD((1),index_Loc)).^(-2);
    PL_BD2RX((1),index_Loc)=1e-3*(Dis_BD2RX((1),index_Loc)).^(-2);
    PL_TX2BD((2),index_Loc)=1e-3*(Dis_TX2BD((2),index_Loc)).^(-2);
    PL_BD2RX((2),index_Loc)=1e-3*(Dis_BD2RX((2),index_Loc)).^(-2);
    PL_TX2BD((3),index_Loc)=1e-3*(Dis_TX2BD((3),index_Loc)).^(-2);
    PL_BD2RX((3),index_Loc)=1e-3*(Dis_BD2RX((3),index_Loc)).^(-2);
    
    for index_BD=1:L
        ag(:,index_BD,index_Loc)=exp(1i*2*pi*1/2*sin_theta(index_BD,index_Loc)*[0:1:Num_ante-1].'); 
        g(:,index_BD,index_Loc)=sqrt(PL_BD2RX(index_BD,index_Loc))*ag(:,index_BD,index_Loc); %%M*L
        h(index_BD,index_Loc)=sqrt(PL_TX2BD(index_BD,index_Loc)); %%L*1
    end
    channel_gain(index_Loc)=sum((sqrt(PL_BD2RX(:,index_Loc)).*sqrt(PL_TX2BD(:,index_Loc))).^2); 
    noise(:,:,:,index_Loc)=randn(N_symbo,Num_symb,Num_ante)+1i*randn(N_symbo,Num_symb,Num_ante);
 end

        sin_theta_TX2RX=0;
        PL_TX2RX=1e-4*(Dis_TX2RX).^(-3);
        af=exp(1i*2*pi*1/2*sin_theta_TX2RX*[0:1:Num_ante-1].');
        f=sqrt(PL_TX2RX)*af;
      
 c_data=randi([0 1],Num_symb-1,L);
 data=randi([0 1],1,log2(M)*N_sc*Num_symb); 
 save channel_L3_M32_240507.mat