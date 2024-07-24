function [ToA] = LoraToA(symPre,symHead,bytePL,bitCRC,SF,CR,B)
% LoraToA --> Calculate the time-on-air of LoRa packet

% List of inputs:
%
% symPre : number of preamble symbols (8 in LoRaWAN)
% symHead : number of header symbols (20 if header is variable, 0 if its fixed)
% bytePL : number of byte of payload
% bitCRC : number of bit of CRC (16 bit if CRC on, 0 otherwise)
% SF : spreading factor (5 to 12)
% CR : coding rate (CR=1,2,3,4 for rates 4/5, 4/6, 4/7 and 4/8 respectively)
% B : bandwidth (kHz)

% Output
%
% ToA : time-on-air (s)


% Calculate the number of symbols according to SX1280 datasheet
if SF<7
    Nsym = symPre+6.25+8+ceil(max(8*bytePL+bitCRC-4*SF+symHead,0)/(4*SF)*(CR+4));
else 
    if SF<=10
        Nsym = symPre+4.25+8+ceil(max(8*bytePL+bitCRC-4*SF+symHead+8,0)/(4*SF)*(CR+4));
    else
        % with low data rate optimization 
        Nsym = symPre+4.25+8+ceil(max(8*bytePL+bitCRC-4*SF+symHead,0)/(4*(SF-2))*(CR+4));
    end
end

% Calculate the time-on-air
ToA = (2^SF)/B*Nsym*1e-3;

end

