%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mohamed Shajith K %
% This is the proposed code for the project based assessment for Computer Networks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%% Initial parameters %%%%
% Dimensions of the network range (in metres) %
x_max=100
y_max=100
x=0
y=0

% x coordinates and y coordinates of the base station %
base_x=0.5*x_max
base_y=0.5*y_max

% Total number of nodes in the field %
n=20; % can also be n=100
% Count of Dead Nodes in the initial state %
nodes_dead=0;

%%% Values for the energy parameters %%%
% Initial Energy for a node (unit in Joules)% 
Eo=0.5;
% Energy required for communication ( units in Joules/bit)%
Eelec=50*10^(-9); 
ETX=50*10^(-9);
ERX=50*10^(-9);

% Energy for the transmitting amplifier(units in Joules/bit/m^2) %
Efs=10*10^(-12)
Emp=0.0013*10^(-12); %  (Energy spent by the amplifier to transmit the bits) %
% Data Aggregation Energy (units in Joules/bit) %
EDA=5*10^(-9); 
% Size of data package (units in bits)%
k=4000;
% Probability for the cluster head approval %
p=0.1; % i.e. 10 % of the total nodes are assumed to give good results
% Number of Clusters %
Nc=p*n; 
% Present count of operating Nodes %
nodes_being_operated=n;
temp_val=0;
yd1=33
%Percentage of the advanced nodes in the network
m=0.2;
%a is alpha (times the energy of advanced nodes are greater than that of normal nodes)
a=1;
% Maximum number of rounds %
rmax=9000
%%%% End of Parameters %%%%
%Computation of do
do=sqrt(Efs/Emp);
Et=0;


            %%% Simulation of a WSN %%%

% WSN plot %
for i=1:1:n
    
    Simu(i).id=i;	% sensor ID number
    Simu(i).xd=rand(1,1)*x_max;	% X-axis coordinates of sensor node
    Simu(i).yd=rand(1,1)*y_max;	% Y-axis coordinates of sensor node
    Simu(i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
    Simu(i).rop=0;	% number of rounds node was operational
    Simu(i).rleft=0;  % rounds left for node to become available for Cluster Head election
    Simu(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    Simu(i).cluster=0;	% the cluster which a node belongs to
    Simu(i).tel=0;	% states how many times the node was elected as a Cluster Head
    Simu(i).rn=0;     % round node got elected as cluster head
    Simu(i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to
    Simu(i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
    Simu(i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
    Simu(i).dts=0;    % nodes distance from the base station
    

    
    hold on;
    figure(1)
    plot(x,y,x_max,y_max,Simu(i).xd,Simu(i).yd,'ob',base_x,base_y,'*r');
    title 'The simulated Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    
end
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LEACH Implementation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:n
    LEACH(i).xd=rand(1,1)*x_max;
    XR(i)=LEACH(i).xd;
    LEACH(i).yd=rand(1,1)*y_max;
    YR(i)=LEACH(i).yd;
    LEACH(i).G=0;
    %initially there are no cluster heads only nodes
    LEACH(i).type='N';
    
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1)
        LEACH(i).E=Eo;
        LEACH(i).ENERGY=0;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)
        LEACH(i).E=Eo*(1+a)
        LEACH(i).ENERGY=1;
    end
end

LEACH(n+1).xd=base_x;
LEACH(n+1).yd=base_y;

%counter for CHs
CH_count_LEACH=0;
%counter for CHs per round
CH_count_perround_LEACH=0;
cluster3=1;

CH_count_LEACH;
CH_count_perround_LEACH=CH_count_perround_LEACH+CH_count_LEACH;
flag_first_dead_LEACH=0;
alive_LEACH=n;
packets_sent_TO_BS_LEACH=0;
packets_TO_CH_LEACH=0;
for r=0:1:rmax    
    %Operation for epoch
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            LEACH(i).G=0;
            LEACH(i).cl=0;
        end
    end        
    
    %Number of dead nodes
    dead_LEACH=0;
    %Number of dead Advanced Nodes
    dead_adv_LEACH=0;
    %Number of dead Normal Nodes
    dead_nrml_LEACH=0;
    tic
    for i=1:1:n
        %checking if there is a dead node
        if (LEACH(i).E<=0)
            dead_LEACH=dead_LEACH+1;
            if(LEACH(i).ENERGY==1)
                dead_adv_LEACH=dead_adv_LEACH+1;
            end
            if(LEACH(i).ENERGY==0)
                dead_nrml_LEACH=dead_nrml_LEACH+1;
            end
        end
        if LEACH(i).E>0
            LEACH(i).type='N';
            if (LEACH(i).ENERGY==0)
            end
            if (LEACH(i).ENERGY==1)
            end
        end
        
        
        STATS.DEAD_LEACH(r+1)=dead_LEACH;
        STATS.ALIVE_LEACH(r+1)=alive_LEACH-dead_LEACH;
    end            
    %When the first node dies
    if (dead_LEACH==1)
        if(flag_first_dead_LEACH==0)
            first_dead_LEACH=r
            flag_first_dead_LEACH=1;
        end
    end    
    CH_count_LEACH=0;
    cluster3=1;
    for i=1:1:n
        if(LEACH(i).E>0)
            temp_rand=rand;
            if ( (LEACH(i).G)<=0)
                
                %Election of cluster3 Heads
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    CH_count_LEACH=CH_count_LEACH+1;
                    packets_sent_TO_BS_LEACH=packets_sent_TO_BS_LEACH+2;
                    % packets_sent_TO_BS_LEACH(r+1)=packets_sent_TO_BS_LEACH;
                    
                    LEACH(i).type='C';
                    LEACH(i).G=round(1/p)-1;
                    C(cluster3).xd=LEACH(i).xd;
                    C(cluster3).yd=LEACH(i).yd;                    
                    distance=sqrt( (LEACH(i).xd-(LEACH(n+1).xd) )^2 + (LEACH(i).yd-(LEACH(n+1).yd) )^2 );
                    C(cluster3).distance=distance;
                    C(cluster3).id=i;
                    X(cluster3)=LEACH(i).xd;
                    Y(cluster3)=LEACH(i).yd;
                    cluster3=cluster3+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        LEACH(i).E=LEACH(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        LEACH(i).E=LEACH(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                    end
                end
                
            end
        end
    end
    
    STATS.CLUSTERHEADLEACH(r+1)=cluster3-1;
    CLUSTERHS(r+1)=cluster3-1;
    
    %Election of Associated Cluster Head for Normal Nodes
    for i=1:1:n
        if ( LEACH(i).type=='N' && LEACH(i).E>0 )
            if(cluster3-1>=1)
                min_dis=sqrt( (LEACH(i).xd-LEACH(n+1).xd)^2 + (LEACH(i).yd-LEACH(n+1).yd)^2 );
                min_dis_cluster=1;
                for c=1:1:cluster3-1
                    temp=min(min_dis,sqrt( (LEACH(i).xd-C(c).xd)^2 + (LEACH(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end
                
                %Energy dissipated by associated Cluster Head
                min_dis;
                if (min_dis>do)
                    LEACH(i).E=LEACH(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    LEACH(i).E=LEACH(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                %Energy dissipated
                if(min_dis>0)
                    LEACH(C(min_dis_cluster).id).E = LEACH(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                    packets_TO_CH_LEACH(r+1)=n-dead_LEACH-cluster3+1;
                end
                
                LEACH(i).min_dis=min_dis;
                LEACH(i).min_dis_cluster=min_dis_cluster;
                
            end
        end
    end   
    
    packets_sent_TO_BS_LEACH=packets_sent_TO_BS_LEACH+1;
    STATS.packets_sent_TO_BS_LEACH(r+1)=packets_sent_TO_BS_LEACH;
    
    CH_count_LEACH;
    CH_count_perround_LEACH=CH_count_perround_LEACH+CH_count_LEACH;
    
end
toc
EndtoEndDelay_LEACH=toc/STATS.packets_sent_TO_BS_LEACH(r+1)
Throughput_LEACH=STATS.packets_sent_TO_BS_LEACH(r+1)/EndtoEndDelay_LEACH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Z - SEP Implementation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1);
for i=1:1:n
    if (i<=5)
        
        Z_SEP(i).xd=rand(1,1)*x_max;
        Z_SEP(i).yd=rand(1,1)*yd1;
        Z_SEP(i).E=Eo*(1+a);
        Z_SEP(i).ENERGY=1;
        Z_SEP(i).type='A';
    end
    if (i>10 &&  i<=100)
        Z_SEP(i).xd=rand(1,1)*x_max;
        Z_SEP(i).yd=33+rand(1,1)*yd1;
        Z_SEP(i).E=Eo;
        Z_SEP(i).ENERGY=0;
        Z_SEP(i).type='N';
    end
    if(i<=10 && i>5)
        Z_SEP(i).xd=rand(1,1)*x_max;
        Z_SEP(i).yd=66+rand(1,1)*yd1;
        Z_SEP(i).E=Eo*(1+a);
        Z_SEP(i).ENERGY=1;
        Z_SEP(i).type='A';
    end
end
S2(n+1).xd=base_x;
S2(n+1).yd=base_y;
flag_first_dead_Zsep=0;
alive_Zsep=n;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS12=0;
packets_TO_BS22=0;
packets_TO_BS_ZSEP=0;
packets_TO_CH2=0;

for r=0:1:rmax
    %Election Probability for Advanced Nodes
    padv= ( p*(1+a)/(1+(a*m)) );
    
    %Operations for sub-epochs
    if(mod(r, round(1/padv) )==0)
        for i=1:1:n
            if(Z_SEP(i).ENERGY==1)
                Z_SEP(i).G=0;
            end            
        end
    end    
    dead_Zsep=0;
    %Number of dead Advanced Nodes
    dead_adv_Zsep=0;
    %Number of dead Normal Nodes
    dead_nrml_Zsep=0;
    tic
    for i=1:1:n
        %checking if there is a dead node
        if (Z_SEP(i).E<=0)            
            dead_Zsep=dead_Zsep+1;
            
            if(Z_SEP(i).ENERGY==1)
                dead_adv_Zsep=dead_adv_Zsep+1;
            end
            if(Z_SEP(i).ENERGY==0)
                dead_nrml_Zsep=dead_nrml_Zsep+1;
            end
        end
        if (Z_SEP(i).E>0)
            Z_SEP(i).type='N';
            if (Z_SEP(i).ENERGY==0)
            end
            if (Z_SEP(i).ENERGY==1)
            end
        end

        STATS.DEAD_ZSEP(r+1)=dead_Zsep;
        STATS.ALIVE_ZSEP(r+1)=alive_Zsep-dead_Zsep;
        
    end
    
    %When the first node dies
    if (dead_Zsep==1)
        if(flag_first_dead_Zsep==0)
            first_dead_Zsep=r
            flag_first_dead_Zsep=1;
        end
    end
    for(i=1:1:n)
        if(Z_SEP(i).E>=0)
            if(Z_SEP(i).type=='N')
                distance=sqrt( (Z_SEP(i).xd-(S2(n+1).xd) )^2 + (Z_SEP(i).yd-(S2(n+1).yd) )^2 );
                if (distance>do)
                    Z_SEP(i).E=Z_SEP(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                end
                if (distance<=do)
                    Z_SEP(i).E=Z_SEP(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                end
                packets_TO_BS12=packets_TO_BS12+1;
            end
        end
        
    end
    
    
    
    countCHs2=0;
    cluster2=1;
    
    for i=1:1:n
        if(Z_SEP(i).E>=0)
            if(Z_SEP(i).G<=0)
                if(Z_SEP(i).type=='A')
                    temp_rand=rand;
                    if( ( Z_SEP(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
                        
                        countCHs2=countCHs2+1;
                        packets_TO_BS22=packets_TO_BS22+1;                                                
                        Z_SEP(i).type='C';
                        Z_SEP(i).G=(1/padv)-1;
                        C(cluster2).xd=Z_SEP(i).xd;
                        C(cluster2).yd=Z_SEP(i).yd;                        
                        distance=sqrt( (Z_SEP(i).xd-(S2(n+1).xd) )^2 + (Z_SEP(i).yd-(S2(n+1).yd) )^2 );
                        C(cluster2).distance=distance;
                        C(cluster2).id=i;
                        X(cluster2)=Z_SEP(i).xd;
                        Y(cluster2)=Z_SEP(i).yd;
                        cluster2=cluster2+1;
                        
                        %          Calculation of Energy dissipated
                        distance;
                        if (distance>do)
                            Z_SEP(i).E=Z_SEP(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                        end
                        if (distance<=do)
                            Z_SEP(i).E=Z_SEP(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                        end
                        
                    end
                end
                
            end
            for i=1:1:n
                if ( Z_SEP(i).type=='A' && Z_SEP(i).E>0 )
                    if(cluster2-1>=1)
                        min_dis=inf;
                        min_dis_clustera=1;
                        
                        for c=1:1:cluster2-1
                            temp=min(min_dis,sqrt( (Z_SEP(i).xd-C(c).xd)^2 + (Z_SEP(i).yd-C(c).yd)^2 ) );
                            
                            if ( temp<min_dis )
                                min_dis=temp;
                                min_dis_cluster2=c;
                            end
                            
                        end
                        
                        %Energy dissipated by  Cluster menmber for transmission of packet
                        if (min_dis>do)
                            Z_SEP(i).E=Z_SEP(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis));
                        end
                        if (min_dis<=do)
                            Z_SEP(i).E=Z_SEP(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                        end
                        %Energy dissipated by clustre head in receving
                        if(min_dis>0)
                            S2(C(min_dis_cluster2).id).E = S2(C(min_dis_cluster2).id).E- ( (ERX + EDA)*4000 );
                            PACKETS_TO_CH(r+1)=n-dead-cluster2+1;
                        end
                        
                        Z_SEP(i).min_dis=min_dis;
                        Z_SEP(i).min_dis_cluster2=min_dis_cluster2;
                        
                    end
                end
            end
        end
        
        CLUSTERHS(r+1)=cluster2+1;
        
        STATS.CLUSTERHEADS2(r+1)=cluster2+1;
        
        
    end
    
    
    packets_TO_BS_ZSEP=packets_TO_BS12+packets_TO_BS22;
    STATS.PACKETS_TO_BS_ZSEP(r+1)=packets_TO_BS_ZSEP;
    
end
toc
EndtoEndDelay_ZSEP=toc/STATS.PACKETS_TO_BS_ZSEP(r+1)
Throughput_ZSEP=STATS.PACKETS_TO_BS_ZSEP(r+1)/EndtoEndDelay_ZSEP
%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%
figure(2);
r=0:rmax;
plot(r,STATS.ALIVE_LEACH(r+1),'-r',r,STATS.ALIVE_ZSEP(r+1),'-g')
legend('LEACH','ZSEP');
xlabel('Number of Rounds');
ylabel('Count of Alive Nodes');

figure(3);
plot(r,STATS.DEAD_LEACH(r+1),'-r',r,STATS.DEAD_ZSEP(r+1),'-g')
legend('LEACH','ZSEP');
xlabel('Number of Rounds');
ylabel('Count of Dead Nodes');

figure(4);
plot(r,STATS.packets_sent_TO_BS_LEACH(r+1),'-r',r,STATS.PACKETS_TO_BS_ZSEP(r+1),'-g')
legend('LEACH','ZSEP');
xlabel('Number of Rounds');
ylabel('Number of Packets that were sent to BS');  

figure(5);
X=categorical({'LEACH','ZSEP'});
X=reordercats(X,{'LEACH','ZSEP'});
Y=[Throughput_LEACH,Throughput_ZSEP]
b=bar(X,Y)
b.FaceColor='flat';
b.CData(2,:)=[.5 0 .5]
legend('LEACH');
legend('ZSEP');
xlabel('Number of Rounds');
ylabel('Throughput');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
