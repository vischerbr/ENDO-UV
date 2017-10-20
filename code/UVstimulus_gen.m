% This function generates the square wave intensity profile that our
% adaption model uses as a stimulus.  
function [UVstimulus] = UVstimulus_gen(duty, period, reps) % unitless inputs

UVstimulus = zeros(1, round(period)); % initialize to off
UVstimulus(1,1:round(period*duty))= 1; % set when intensity is on

UVstimulus = repmat(UVstimulus, 1, reps); % repeat for total length required
end