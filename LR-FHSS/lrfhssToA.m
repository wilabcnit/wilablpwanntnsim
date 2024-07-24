function [ToA,Nfrag,Nhops,gamma] = lrfhssToA(N_H,T_H,T_F,PL,Tlast,CR)
% LoraToA --> Calculate the time-on-air of LoRa packet

% List of inputs:
%
% N_H : N° of headers replicas
% T_H : time-on-air of 1 header (s) 
% T_F : time-on-air of 1 Payload fragment in seconds 
% bytePL : number of byte of payload
% Tlast : time-on-air of the last Payload fragment in seconds [not used for now]
% CR : code rate (1 for code rate 1/3, 2 for code rate 2/3)


% Outputs
%
% ToA : time-on-air (s)
% Nfrag : number of payload fragments
% Nhops : number of frequency hops
% gamma : number of fragments that have to be received without collisions


% Calculate the number of fragments
switch CR
    case 1
        M=2;
        Nfrag=ceil((PL+3)/M);
        gamma = ceil(1/3*Nfrag); 
    case 2
        M=4;
        Nfrag=ceil((PL+3)/M);
        gamma = ceil(2/3*Nfrag);
    otherwise
        error('Code rate not supported.')
end


% Calculate the time-on-air
ToA = (N_H*T_H+Nfrag*T_F);

% Calculate the number of hops
Nhops = N_H+Nfrag;

end

