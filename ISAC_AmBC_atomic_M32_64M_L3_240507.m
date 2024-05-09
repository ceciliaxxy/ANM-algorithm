%clc

%%%�����˵�����ֹ�������޸Ĳ���rho,tau��ο����� superresolution һ��
clear all;
load channel_L3_M32_240507.mat
tic
N_fft=64;           % FFT ����
Space_pilot=4;       %��Ƶ���
N_pilot=N_fft/Space_pilot;    %��Ƶ��Ŀ 32
N_sc=N_fft-N_pilot;            %���ز��� 96
N_cp=N_fft/4;             % ѭ��ǰ׺���ȡ�Cyclic prefix 32
N_symbo=N_fft+N_cp;        % 1������OFDM���ų��� 160
deltaf=1e6;       %���ز����  500K,64M
B=deltaf*N_fft;         %�ܴ���  
ts=1/B;                %OFDM���Ų������
M=4;              %4PSK����

SNR_bB=[20];        %��������� dB
SNR_ratio=10.^(SNR_bB/10);  %��������� 
power_noi_dB=-114;
power_noi=10.^(power_noi_dB/10); %%��������,ÿ�������ϵ��������ʣ�power_sigΪÿ�����������ϵĹ��� 
L=3;
Num_symb=64;              % ÿ֡������OFDM������
c=3*1e8;
Num_BD=L;
factor_oversamp=2^5;
Num_channel=1;
%% �����źŲ���
%�����������ݲ���
data=randi([0 1],1,log2(M)*N_sc*Num_symb);  %%���ɴ��� 96*6*2=768 bit
%%qpsk����
data_resh= reshape(data,log2(M),[])';             %��ÿ��2���ؽ��з��飬M=4, 384*2
data_deci= bi2de(data_resh);                             %������ת��Ϊʮ����, 384*1
data_modu=pskmod(data_deci,M,pi/M);               %4PSK����
%scatterplot(modu_data),grid;                  %����ͼ(Ҳ����ȡʵ����plot����)
%pilot���ݲ���
pilot=ones(1,N_pilot*Num_symb);
data_addpilot=zeros(N_fft*Num_symb,1);
%%���뵼Ƶ
pilot_loc=[]; %ָʾ��Ƶ�����������е�λ��
ip=0;   % ipָʾ�˵�ǰOFDM�������Ѿ�����ĵ�Ƶ������
      for k=1:N_fft*Num_symb % ��Ƶ����ض�λ�ü��뵼Ƶ������
         if mod(k,Space_pilot)==1
           data_addpilot(k) = pilot(floor(k/Space_pilot)+1); pilot_loc = [pilot_loc k]; ip = ip+1;
         else
           data_addpilot(k) = data_modu(k-ip);
         end
      end

%%����ת��
data_seq=reshape(data_addpilot,N_fft,Num_symb); %128*6

%%IFFT�任
%ifft_data=ifft(data_seq,N_fft); %128*6
tau_back=(Dis_TX2BD+Dis_BD2RX)/c; %L*1  


%%BD��������%%��ֱ���
c_data=randi([0 1],Num_symb-1,Num_BD);
c_data_diff_addrefe=ones(Num_symb,Num_BD);%%(Num_symb)*Num_BD
for index_BD=1:Num_BD
    for index_symb=1:Num_symb-1
        c_data_diff_addrefe(index_symb+1,index_BD)=xor(c_data_diff_addrefe(index_symb,index_BD),c_data(index_symb,index_BD));  %%Num_symb*NUm_BD      
    end
end
c_data_minus=(-2)*c_data_diff_addrefe+1;

H_symb=zeros(N_pilot,Num_ante,Num_symb);
Num_disablechannel=0;
Num_usechannel=Num_channel;
ber_s=zeros(length(SNR_ratio),Num_channel);
Ber_c=zeros(length(SNR_ratio),Num_channel);
error_theta=zeros(Num_BD,length(SNR_ratio),Num_channel);
error_theta_sumBD=zeros(length(SNR_ratio),Num_channel);
sintheta_esti=zeros(Num_BD,length(SNR_ratio),Num_channel);
error_R=zeros(Num_BD,length(SNR_ratio),Num_channel);
error_R_sumBD=zeros(length(SNR_ratio),Num_channel);
R_esti=zeros(Num_BD,length(SNR_ratio),Num_channel);
theta_esti=zeros(Num_BD,length(SNR_ratio),Num_channel);
alpha_c_esti=zeros(Num_BD,Num_symb,length(SNR_ratio),Num_channel);
Tx_back_delay=zeros(Num_BD,Num_symb*N_symbo);
row=zeros(Num_BD,length(SNR_ratio),Num_channel);
col=zeros(Num_BD,length(SNR_ratio),Num_channel);
chan_back=zeros(Num_ante,Num_symb*N_symbo,Num_BD);
Rx_back=zeros(Num_ante,Num_symb*N_symbo,Num_BD);
Rx_para=zeros(N_symbo,Num_symb,Num_ante);
Rx_CPremo=zeros(N_fft,Num_symb,Num_ante);
fft_data=zeros(N_fft,Num_symb,Num_ante);
c_esti=zeros(Num_BD,Num_symb-1,length(SNR_ratio),Num_channel);
s_esti=zeros(N_sc,Num_symb);
H_esti=zeros(Num_ante,N_fft,Num_symb);
Loc_channel=[];
Num_Ind_SBL=500;

for index_channel=1:Num_channel 

    index_channel
    
    for index_BD=1:Num_BD
        delay_vec(:,index_BD)=exp(-1i*2*pi*deltaf*tau_back(index_BD,index_channel)*[0:N_fft-1]);  %%128*3
        delay_mat(:,:,index_BD)=diag(delay_vec(:,index_BD));
        data_delay(:,:,index_BD)=delay_mat(:,:,index_BD)*data_seq; %%128*6*3
        ifft_data(:,:,index_BD)=ifft(data_delay(:,:,index_BD),N_fft); %%128*6*3
        %%���뱣�������ѭ��ǰ׺
        data_addCP(:,:,index_BD)=[ifft_data(N_fft-N_cp+1:end,:,index_BD);ifft_data(:,:,index_BD)];%% ���뱣�������ѭ��ǰ׺  %%160*6

        %%����ת��
        Tx_data(index_BD,:)=reshape(data_addCP(:,:,index_BD),1,[]);%���ڴ�����Ҫ 3*960
    
    end
          %%ֱ����·�����ź�
        tau_di=Dis_TX2RX/c;
        delay_vec_di=exp(-1i*2*pi*deltaf*tau_di*[0:N_fft-1]);  %%128*3
        delay_mat_di=diag(delay_vec_di);
        data_delay_di=delay_mat_di*data_seq; %%128*6*3
        ifft_data_di=ifft(data_delay_di,N_fft); %%128*6*3
         %%���뱣�������ѭ��ǰ׺
        data_addCP_di(:,:)=[ifft_data_di(N_fft-N_cp+1:end,:);ifft_data_di(:,:)];%% ���뱣�������ѭ��ǰ׺  %%160*6

        %%����ת��
        Tx_data_di=reshape(data_addCP_di(:,:),1,[]);%���ڴ�����Ҫ 3*960  

    for index_SNR= 1:length(SNR_ratio)
%% �����ŵ�
     
        changainBD=kron(c_data_minus,ones(N_fft+N_cp,1));%%960*3
        for index_BD=1:Num_BD
            chan_back(:,:,index_BD)=h(index_BD,index_channel)*g(:,index_BD,index_channel)*(changainBD(:,index_BD).'); %%M*960*L,���OFDM symbol���ŵ���ͬ�ģ����ź�ʱ�ӣ������ŵ�ϵ����
            Rx_back(:,:,index_BD)=chan_back(:,:,index_BD)*diag(Tx_data(index_BD,:)); %%M*960*L   
        end
        chan_di=f*ones(1,(N_fft+N_cp)*Num_symb);
        Rx_di=chan_di*diag(Tx_data_di);        
        Rx_back_sum=sum(Rx_back(:,:,:),3)+Rx_di; %%����ֱ����·�źſ�����ȥ��
   
%% ���ջ�����
    for index_num_ante=1:Num_ante
        %%����ת��
        Rx_para(:,:,index_num_ante)=reshape(Rx_back_sum(index_num_ante,:),N_fft+N_cp,[]); %%160*6*Num_ante
        
        %%ȥ�����������ѭ��ǰ׺
        Rx_CPremo(:,:,index_num_ante)=Rx_para(N_cp+1:end,:,index_num_ante); %%128*6*Num_ante
        
        %%FFT
        fft_data_noifree(:,:,index_num_ante)=fft(Rx_CPremo(:,:,index_num_ante)); %128*6*Num_ante              
    end
    
    power_noi(index_SNR,index_channel)=channel_gain(index_channel)/SNR_ratio(index_SNR);  % ���ݷ�����·SNR������������
    fft_data=fft_data_noifree+sqrt(power_noi(index_SNR,index_channel))/sqrt(2)*noise(:,:,:,index_channel);
    
    delay_vec_di_knw=exp(-1i*2*pi*Dp*deltaf*tau_di*[0:N_pilot-1]);
    Rx_di_knw=f*delay_vec_di_knw;
    %%λ�� �ŵ�����
    for index_sym=1:Num_symb
        for index_Num_ante=1:Num_ante        
            H_symb_backdi(:,index_Num_ante,index_sym)=fft_data(pilot_loc(1:N_pilot),index_sym,index_Num_ante);  %%64*6*Num_anteֻ���˵�һ��symbol����λ���ŵ�����            
        end
        H_symb_back(:,:,index_sym)=H_symb_backdi(:,:,index_sym)-Rx_di_knw.';
    end   
 
    H_symb=H_symb_back(:,:,1);
    
   %% atomic alogrithm
   Num_loop_atom=100;
   H_symb_atom=reshape(H_symb,[],1); 
   rho_atom=0.05;
   tau_atom=sqrt(power_noi(index_SNR,index_channel))*sqrt(Num_ante*N_pilot*log(Num_ante*N_pilot));   %%atomic norm weight 
   %%Initialize
   xi_atom=zeros(Num_ante*N_pilot,1);
   U_atom=zeros(2*N_pilot-1,2*Num_ante-1);
   Toep_U=zeros(Num_ante*N_pilot,Num_ante*N_pilot);
   t_atom=0;
   Upsilon_atom=zeros(Num_ante*N_pilot+1,Num_ante*N_pilot+1);
   Upsilon_atom_0=Upsilon_atom(1:Num_ante*N_pilot,1:Num_ante*N_pilot);
   Upsilon_atom_1=Upsilon_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1);
   Upsilon_atom_bar=Upsilon_atom(Num_ante*N_pilot+1,Num_ante*N_pilot+1);
   
   Lambda_atom=zeros(Num_ante*N_pilot+1,Num_ante*N_pilot+1);
   Lambda_atom_0=Lambda_atom(1:Num_ante*N_pilot,1:Num_ante*N_pilot);
   Lambda_atom_1=Lambda_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1);
   Lambda_atom_bar=Lambda_atom(Num_ante*N_pilot+1,Num_ante*N_pilot+1);
     
   for index_atom=1:Num_loop_atom   
       
       %%xi_atom update
       xi_atom_old=xi_atom;
       xi_atom=1/(1+2*rho_atom)*(2*rho_atom*(Upsilon_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1)+Upsilon_atom(Num_ante*N_pilot+1,1:Num_ante*N_pilot)')/2 ...
                      +2*(Lambda_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1)+Lambda_atom(Num_ante*N_pilot+1,1:Num_ante*N_pilot)')/2+H_symb_atom);
       %Lagaugument=1/2*norm(H_symb_atom-xi_atom)^2+tau_atom/2*U_atom(N_pilot,Num_ante)+tau_atom/2*t_atom+real(trace((Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom])*Lambda_atom'))...
          % + rho_atom/2*norm(Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom],'fro')^2 
       % norm(xi_atom-H_symb_atom)
       
       if norm(xi_atom_old-xi_atom)/norm(xi_atom)<1e-2 || index_atom>=Num_loop_atom
           break;
       end
              
       %%t_atom update
       t_atom=1/rho_atom*(Lambda_atom_bar-tau_atom)+Upsilon_atom_bar;
       %Lagaugument=1/2*norm(H_symb_atom-xi_atom)^2+tau_atom/2*U_atom(N_pilot,Num_ante)+tau_atom/2*t_atom+real(trace((Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom])*Lambda_atom'))...
           %+ rho_atom/2*norm(Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom],'fro')^2 
       
       %%U_atom update
       for k_atom=-N_pilot+1:N_pilot-1
           for m_atom=-Num_ante+1:Num_ante-1                             
               if m_atom<0
                   for l_atom=0:(Num_ante-abs(m_atom)-1)  %%l_atom������ U_atom(k_atom+N_pilot,m_atom+Num_ante)�ĵ�l_atom���Ӿ���|l_atom|�ܹ��� Num_ante-abs(m_atom)
                       trace_Upsilon(l_atom+1)=sum(diag(Upsilon_atom_0(l_atom*N_pilot+1:(l_atom+1)*N_pilot,(l_atom+abs(m_atom))*N_pilot+1:(l_atom+abs(m_atom)+1)*N_pilot),-k_atom));
                       trace_Lambda(l_atom+1)=sum(diag(Lambda_atom_0(l_atom*N_pilot+1:(l_atom+1)*N_pilot,(l_atom+abs(m_atom))*N_pilot+1:(l_atom+abs(m_atom)+1)*N_pilot),-k_atom));
                   end
               end
               
               if m_atom>=0
                   for l_atom=0:(Num_ante-abs(m_atom)-1)
                       trace_Upsilon(l_atom+1)=sum(diag(Upsilon_atom_0((l_atom+m_atom)*N_pilot+1:(l_atom+m_atom+1)*N_pilot,l_atom*N_pilot+1:(l_atom+1)*N_pilot),-k_atom));
                       trace_Lambda(l_atom+1)=sum(diag(Lambda_atom_0((l_atom+m_atom)*N_pilot+1:(l_atom+m_atom+1)*N_pilot,l_atom*N_pilot+1:(l_atom+1)*N_pilot),-k_atom));
                   end                 
               end
               trace_Upsilon_sum=sum(trace_Upsilon(1:Num_ante-abs(m_atom)));
               trace_Lambda_sum=sum(trace_Lambda(1:Num_ante-abs(m_atom)));
               U_atom(k_atom+N_pilot,m_atom+Num_ante)=1/(rho_atom*(Num_ante-abs(m_atom))*(N_pilot-abs(k_atom)))*(rho_atom*trace_Upsilon_sum+trace_Lambda_sum);               
           end
       end
       U_atom(N_pilot,Num_ante)= U_atom(N_pilot,Num_ante)-1/(rho_atom*Num_ante*N_pilot)*tau_atom;
       for l_atom_x=1:Num_ante
           for l_atom_y=1:Num_ante
                Toep_U((l_atom_x-1)*N_pilot+1:l_atom_x*N_pilot,(l_atom_y-1)*N_pilot+1:l_atom_y*N_pilot)=toeplitz(U_atom(N_pilot:2*N_pilot-1,l_atom_x-l_atom_y+Num_ante),U_atom([N_pilot:-1:1],l_atom_x-l_atom_y+Num_ante));
           end
       end 
      % Lagaugument=1/2*norm(H_symb_atom-xi_atom)^2+tau_atom/2*U_atom(N_pilot,Num_ante)+tau_atom/2*t_atom+real(trace((Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom])*Lambda_atom'))...
           %+ rho_atom/2*norm(Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom],'fro')^2 
       
      
       %%Upsilon_atom update
%        for l_atom_x=1:Num_ante
%            for l_atom_y=1:Num_ante
%                 Toep_U((l_atom_x-1)*N_pilot+1:l_atom_x*N_pilot,(l_atom_y-1)*N_pilot+1:l_atom_y*N_pilot)=toeplitz(U_atom(N_pilot:2*N_pilot-1,l_atom_x-l_atom_y+Num_ante),U_atom([N_pilot:-1:1],l_atom_x-l_atom_y+Num_ante));
%            end
%        end      
        Q_atom=[Toep_U xi_atom;xi_atom' t_atom]-1/rho_atom*Lambda_atom;
        [vec_Upsilon_atom_inte,eig_Upsilon_atom_inte]=eig(Q_atom/2+Q_atom'/2);
        eig_Upsilon_atom_inte(eig_Upsilon_atom_inte<0)=0;
        Upsilon_atom=vec_Upsilon_atom_inte*eig_Upsilon_atom_inte*vec_Upsilon_atom_inte';
        Upsilon_atom=Upsilon_atom/2+Upsilon_atom'/2;
        Upsilon_atom_0=Upsilon_atom(1:Num_ante*N_pilot,1:Num_ante*N_pilot);
        Upsilon_atom_1=Upsilon_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1);
        Upsilon_atom_bar=Upsilon_atom(Num_ante*N_pilot+1,Num_ante*N_pilot+1);
        
       %Lagaugument=1/2*norm(H_symb_atom-xi_atom)^2+tau_atom/2*U_atom(N_pilot,Num_ante)+tau_atom/2*t_atom+real(trace((Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom])*Lambda_atom'))...
           %+ rho_atom/2*norm(Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom],'fro')^2 
      
        %%Lambda_atom update
        Lambda_atom=Lambda_atom+rho_atom*(Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom]);
        Lambda_atom_0=Lambda_atom(1:Num_ante*N_pilot,1:Num_ante*N_pilot);
        Lambda_atom_1=Lambda_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1);
        Lambda_atom_bar=Lambda_atom(Num_ante*N_pilot+1,Num_ante*N_pilot+1);
        
       %Lagaugument=1/2*norm(H_symb_atom-xi_atom)^2+tau_atom/2*U_atom(N_pilot,Num_ante)+tau_atom/2*t_atom+real(trace((Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom])*Lambda_atom'))...
           %+ rho_atom/2*norm(Upsilon_atom-[Toep_U xi_atom;xi_atom' t_atom],'fro')^2  
       
   end
     dualsolu=-1/4*(Lambda_atom(1:Num_ante*N_pilot,Num_ante*N_pilot+1)+Lambda_atom(Num_ante*N_pilot+1,1:Num_ante*N_pilot)');   
    
    Num_ante_over=factor_oversamp*Num_ante;
    i_ante=[0:Num_ante_over-1]';
    j_ante=[0:Num_ante-1]';
    ij_ante=j_ante*i_ante'; %128*512
    WN_AoA=exp(1i*2*pi/Num_ante_over);
    DFTmatrix_AoA=WN_AoA.^ij_ante;  %128*512
    
    Num_pilot_over=factor_oversamp*N_pilot;
    i_pilot=[0:Num_pilot_over-1]';
    j_pilot=[0:N_pilot-1]';
    ij_pilot=j_pilot*i_pilot'; %32*128
    WN_delay=exp(-1i*2*pi/Num_pilot_over);
    IDFTmatrix_delay=WN_delay.^ij_pilot;  %32*128
    
    search_atom=kron(DFTmatrix_AoA,IDFTmatrix_delay); %%������Num_ante��N_pilot
    for index_search=1:Num_ante_over*Num_pilot_over
        modu_search(index_search)=abs(search_atom(:,index_search)'*dualsolu);
    end
    search_twodomain=reshape(modu_search,Num_pilot_over,Num_ante_over);   
   % modu_search_tau =1./abs(modu_search-tau_atom);
                        
    %% 3Dͼ    
    R_esti_vetex=c/Space_pilot/deltaf/Num_pilot_over*(0:N_pilot*factor_oversamp-1);  %128    
    sintheta_esti_vetex=2/Num_ante_over*(0:Num_ante_over-1);   % 512/2+1=(Num_ante*factor_oversamp/2+1)
    theta_esti_vetex=asin(sintheta_esti_vetex)*180/pi;
    theta=asin(sin_theta)*180/pi;
    amp = 1/tau_atom;
    search_twodomain_norm=amp*search_twodomain;
    figure(1)
    surf(real(theta_esti_vetex),real(R_esti_vetex),abs(search_twodomain_norm));
    colorbar
    ylabel('range (m)');
    xlabel('AoA (��)')
    zlabel('Intensity')
        hold on
    stem3(rad2deg(asin(sin_theta(:,1))),Dis_TX2BD(:,1)+Dis_BD2RX(:,1),max(max(abs(search_twodomain_norm)))*ones(Num_BD,1),'r*')
     
     %% search BD point   
       position= FastPeakFind(abs(search_twodomain_norm));
   if length(position)<6 && length(position)>=4
      position([5 6])=[Num_ante_over/2, Num_pilot_over/2];
   elseif length(position)<4 && length(position)>=2
       position([3 4 5 6])=[Num_ante_over/2, Num_pilot_over/2,Num_ante_over/2, Num_pilot_over/2];
   elseif length(position)<2
       position([1 2 3 4 5 6])=[Num_ante_over/2, Num_pilot_over/2,Num_ante_over/2, Num_pilot_over/2,Num_ante_over/2, Num_pilot_over/2];
   elseif length(position)>6      %%����ֵ����Ŀ����6����ѡ������������ֵ��
       position_resh=reshape(position,2,[]).';
       position_resh_value=search_twodomain_norm(position_resh(:,2),position_resh(:,1));
       [max_value,maxposition]=sort(diag(position_resh_value));
       iddidv=position_resh(maxposition(end-2:end),:).';
       position=reshape(iddidv,[],1);
   end
      idv=position([1 3 5]);
      idd=position([2 4 6]);
      
     sintheta_half_esti(:,index_SNR,index_channel)=((idv-1).'/Num_ante_over);  
     sintheta_esti(:,index_SNR,index_channel)=2*sintheta_half_esti(:,index_SNR,index_channel);
     tau_esti(:,index_SNR,index_channel)=1/deltaf/Space_pilot*((idd-1).'/Num_pilot_over);
     R_esti(:,index_SNR,index_channel)=tau_esti(:,index_SNR,index_channel)*c;
     
    end
end 

toc









