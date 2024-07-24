%% LR-FHSS Satellite-IoT Uplink Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Enrico Testi, Enrico Paolini
% Updated: 19/07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please reference:
% [E. Testi and E. Paolini, “Packet Collision Probability of
% Direct-to-Satellite IoT Systems,” IEEE Internet of Things Journal,
% under review.]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc


%% Main parameters
MC = 2000; %Monte Carlo iterations
Nstart = 0; % starting total number of interferers
Nstep = 500; % simulation step
Nend = 4000; % ending total number of nodes interferers
min_el = 55; % minimum elevation angle [deg]
d = 600e3; % sat orbit [m]
v = 7.5e3; % sat speed [m/s]
PL = 100;   % Data Payload in bytes

%% Reference node is in position (a,0), with a<L
n = Nstart:Nstep:Nend;
L = d*cot(min_el*pi/180); % max spot radius
a=0; % origin of the xOy plane
% Alternatives are:
% a=L/4;
% a=L/2;
% a=sqrt(L^2-(v^2*Time_on_Air^2)/4)-1; % near the boundary of the spot (minimum contact time)
target_theta = acos(a/L); % angle between target node and the radius parallel to x axis
Ac = (4*sin(target_theta)+pi)*L^2; % Surface covered by the circular spot during the pass
lambda = n./Ac; % Densities of Poisson process

%% Modulations Settings
N_replica = 2; % N° of headers replicas
T_H = 233.472e-3; % TimeOnAir of 1 header in seconds 
T_F = 102.4e-3; % TimeOnAir of 1 Payload fragment in seconds 
Tlast = 102.4e-3; % TimeOnAir of the last Payload fragment in seconds [not used for now]
CR = 2; % Code Rate (1 for code rate 1/3, 2 for code rate 2/3)
Ch = 35; % N° of available channels    
[Time_on_Air,Nfrag,Nhops,gamma] = lrfhssToA(N_replica,T_H,T_F,PL,Tlast,CR);

%% Output
outputFolder = 'resultsLRFHSS/';

%%% Configure output
outputFileName = [outputFolder 'a_' num2str(round(a)) '_Nh_' num2str(N_replica) '_CR_' num2str(CR) '_Ch_' num2str(Ch) '.txt'];
fileID = fopen(outputFileName, 'w');

%% SIMULATION
% Loop over the densities
for i=1:length(n)
    
    %% Initialization
    collision = [];
    collision_H = [];
    collision_F = [];


    %% Monte Carlo loop
    for mc=1:MC
        % Poisson inteferers generation
        N = poissrnd(n(i)); % Draw from Poisson distribution

        %%% Initialize nodes
        pos = [];
        WToA = [];
        txs_time = [];
        colliders = [];
        txs_time_t_frag = zeros(1,Nhops+1);
        txs_time_frag = zeros(N,Nhops+1);
        sel_txs_time = [];
        sel_txs_time_shift = [];
        head_coll = zeros(1,N_replica);
        frag_coll = zeros(1,Nfrag);

        %Target node
        target_init = [L*cos(target_theta) 0]; % starting position of the target 
        l_t = 2*L*sin(target_theta); % length of the segment associated to the target
        Tl_t = l_t/v; % visibility time of the target (time necessary to travel along the whole segment)
        
        if Tl_t-Time_on_Air>=0 % Check if the target has enough time to transmit its packet
            txs_time_t = ((2*L*sin(target_theta))/v-Time_on_Air).*rand()-(L*sin(target_theta))/v; % tx start time for target node
            %Interferers
            pos = throw_nodes_circle_spot(N,L,target_theta,0); % starting position of the interferers   
            % calculate the length of the segment associated to each interferer
            WToA = [(pos(:,2)-sqrt(L^2-pos(:,1).^2))/v (pos(:,2)+sqrt(L^2-pos(:,1).^2))/v-Time_on_Air];
            txs_time = (WToA(:,2)-WToA(:,1)).*rand(N,1)+WToA(:,1); % tx start time for each interferer

            % Detection of a collision for FHSS
            % assign a channel to each header replica and fragment of the reference user
            ref_channels = randi(Ch,[1 Nhops]);

            % vector containing the tx starting times for each
            % fragment of reference user
            ref_idxs = 1:1:length(ref_channels);
            txs_time_t_frag(1:N_replica) = txs_time_t+(ref_idxs(1:N_replica)-1)*T_H;
            txs_time_t_frag(N_replica+1:Nhops) = txs_time_t_frag(N_replica)+T_H+(ref_idxs(N_replica+1:end)-N_replica)*T_F;
            txs_time_t_frag(Nhops+1) = txs_time_t_frag(Nhops)+T_F;

            % assign a channel to each header replica and
            % fragment of the interferers
            int_channels = randi(Ch,[N Nhops]);

            % vector containing the tx starting times for each
            % fragment of interferers
            int_idxs = 1:1:size(int_channels,2);
            txs_time_frag(:,1:N_replica) = txs_time+(int_idxs(1:N_replica)-1)*T_H;
            txs_time_frag(:,N_replica+1:Nhops) = txs_time_frag(:,N_replica)+T_H+(int_idxs(N_replica+1:end)-N_replica)*T_F;
            txs_time_frag(:,Nhops+1) = txs_time_frag(:,Nhops)+T_F;

            % check for each header replica and fragment of reference node if
            % there is a collision
            for r=1:Nhops
               sel_txs_time_t = txs_time_t_frag(r);
               sel_end_time_t = txs_time_t_frag(r+1);
               sel_txs_time = txs_time_frag(int_channels==ref_channels(r));
               sel_txs_time_shift = txs_time_frag(logical(circshift([int_channels==ref_channels(r) zeros(N,1)],[0 1])));
               if r<=N_replica
                   h_rep_coll = nnz(and((sel_txs_time_shift>=sel_txs_time_t),(sel_txs_time<=(sel_end_time_t)))); 
                   if h_rep_coll>0
                       head_coll(r)=1;  
                   else
                       head_coll(r)=0;
                   end
                   if sum(head_coll)==N_replica
                       break;
                   end
               else
                   f_coll = nnz(and((sel_txs_time_shift>=sel_txs_time_t),(sel_txs_time<=(sel_end_time_t)))); 
                   if f_coll>0
                       frag_coll(r-N_replica)=1;
                   else
                       frag_coll(r-N_replica)=0;
                   end
                   if sum(frag_coll>=(Nfrag-gamma))
                       break;
                   end
               end
            end
  
            % decide if packet has been received correctly (at least 1 of the headers + gamma fragments)
            if sum(head_coll)==N_replica
                collision_H(mc)=true;
            else
                collision_H(mc)=false;
            end
            if sum(frag_coll>=(Nfrag-gamma))
                collision_F(mc)=true;
            else
                collision_F(mc)=false;
            end
            if collision_H(mc)+collision_F(mc)>0
                collision(mc)=true;
            else
                collision(mc)=false;
            end              
        else
            collision(mc) = true;
        end
    end

    %% Simulation Results
    % Header no collision prob.
    PH(i) = 1-nnz(collision_H)/length(collision_H);
    % Fragment no collision prob.
    PF(i) = 1-nnz(collision_F)/length(collision_F);
    % Collision prob.
    PC(i) = nnz(collision)/length(collision);
    % No collision prob.
    PnoC(i) = 1-PC(i);
    fprintf(fileID, 'Nodes: %i ; PnoC: %f\n', n(i),PnoC(i));
    fprintf('Nodes: %i ; PnoC: %f\n', n(i),PnoC(i));
end

%% Analytical model
S1 = 2*v*L*(T_H*(2*N_replica+Nfrag)+T_F*Nfrag)/(Ac*Ch);
S2 = 2*v*L*(T_H*(N_replica+2*Nfrag-3)+T_F*(5-3*Nfrag))/(Ac*Ch^2);
Theta = 2*S2/S1-floor(2*S2/S1);
alpha = 1-(Theta*S1^2)/((2-Theta)*S1+2*S2)-((1-Theta)*S1^2)/((1-Theta)*S1+2*S2);
if N_replica==2
    PSanalytical = -exp(n*(alpha^2-1))+2*exp(n*(alpha-1));
elseif N_replica==3
    PSanalytical = -exp(n*(alpha^3-1))-3*exp(n*(alpha^2-1))+3*exp(n*(alpha-1));
end

%% Store simulation results
save([outputFolder '_a_' num2str(round(a)) '_PL_' num2str(PL) '_Nh_' num2str(N_replica) '_CR_' num2str(CR) '_Ch_' num2str(Ch) '.mat'],...
        'PnoC','PSanalytical','lambda','Nstart','Nend','Nstep','n','a','Ac','L',...
        'gamma','Nfrag','N_replica','T_H','T_F','Ch','Time_on_Air','v','PH','PF')
fclose(fileID);

%% Plot results
plot(n,PSanalytical,'b')
hold on
plot(n,PnoC,'*r')
legend('Analytical','Monte Carlo')
xlabel('Average number of interferers')
ylabel('P(S)')

%% Functions
function [positions] = throw_nodes_circle_spot(N,L,theta,distr)
        positions = zeros(N,2);
        if distr==0
            for n=1:N
                positions(n,:) = [2*L*rand()-L 2*L*(1+sin(theta))*rand()-L*(1+sin(theta))];
                while (positions(n,2)>L*sin(theta) && sqrt((positions(n,2)-L*sin(theta))^2+positions(n,1)^2)>L)||(positions(n,2)<-L*sin(theta) && sqrt((positions(n,2)+L*sin(theta))^2+positions(n,1)^2)>L)
                    positions(n,:) = [2*L*rand()-L 2*L*(1+sin(theta))*rand()-L*(1+sin(theta))];
                end
            end  
        else
            positions(:,1)=2*L*rand([N,1])-L;
            positions(:,2)=((2*L*sin(theta)+2*sqrt(L^2-positions(:,1).^2)).*rand([N,1])-L*sin(theta)-sqrt(L^2-positions(:,1).^2));
        end
end
