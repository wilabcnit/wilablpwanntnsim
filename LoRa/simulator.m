%% LoRa Satellite-IoT Uplink Simulator
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
MC = 1e4; %Monte Carlo iterations
Nstart = 100; % starting total number of interferers
Nstep = 100; % simulation step
Nend = 4000; % ending total number of nodes interferers
min_el = 55; % minimum elevation angle [deg]
d = 600e3; % sat orbit [m]
v = 7.5e3; % sat speed [m/s]
PL = 56;   % Data Payload in bytes

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
symPre=8; % number of preamble symbols (8 in LoRaWAN)
symHead=20; % number of header symbols (20 if header is variable, 0 if its fixed)
bitCRC=16; % number of bit of CRC (16 bit if CRC on, 0 otherwise)
SF=7; % spreading factor (5 to 12)
CR=1; % coding rate (CR=1,2,3,4 for rates 4/5, 4/6, 4/7 and 4/8 respectively)
B=125; % bandwidth (kHz)
Time_on_Air = LoraToA(symPre,symHead,PL,bitCRC,SF,CR,B); % % Time on Air for target node (depends on BW and SF) (s)

%% Output
outputFolder = 'resultsLoRa/'; 

%%% Configure output
outputFileName = [outputFolder 'SF_' num2str(SF) '_CR_' num2str(CR) '_B_' num2str(B) '_a_' num2str(round(a)) '.txt'];
fileID = fopen(outputFileName, 'w');


%% SIMULATION
% Loop over the densities
for i=1:length(n)
    collision = zeros(1,MC); % Initialize collision counter 
    %% Monte Carlo loop
    for mc=1:MC

        % Poisson inteferers generation
        N = poissrnd(n(i)); % Draw from Poisson distribution

        %%% Initialize nodes
        pos = [];
        WToA = [];
        txs_time = [];
        colliders = [];

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
            colliders = txs_time(and((txs_time>=txs_time_t-Time_on_Air),(txs_time<=(txs_time_t+Time_on_Air))));
            collision(mc) = length(colliders)>0;          
        else
            collision(mc) = true;
        end

        % Store partial results
        PERpart = nnz(collision(1:mc))/mc; % collision probability in current MC iteration
        PSpart = 1-PERpart;
    end

    %% Simulation Results
    % Collision prob.
    PC(i) = nnz(collision)/length(collision);
    % No collision prob.
    PnoC(i) = 1-PC(i);
    fprintf(fileID, 'Nodes: %i ; PnoC: %f\n', n(i),PnoC(i));
    fprintf('Nodes: %i ; PnoC: %f\n', n(i),PnoC(i));
end

%% Analytical model
PSanalytical = exp(-4*L*Time_on_Air*v*lambda);

%% Store simulation results
save([outputFolder 'SF_' num2str(SF) '_CR_' num2str(CR) '_B_' num2str(B) '_a_' num2str(round(a)) '.mat'],'PnoC','PSanalytical','lambda','Nstart','Nend','Nstep','n','a','Ac','L','Time_on_Air','v');
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