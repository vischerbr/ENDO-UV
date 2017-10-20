% This function creates the force due to UV excitation. We use a model that
% includes adaption based on observed cell response. These equations were
% adapted from [citation]. We posit that the force the particle in the
% double well feels from the UV light will be proportional to the rate of
% population change in a two-level system, where we think of the levels as
% "unbleached" and "bleached". Solutions to this system have the correct
% profile; with a duty cycle of unity the force decays to zero [etc].

function [forceprofile] = forceprofile_gen(pars)

% Initialize paramaters to default in case they are missing. Greek letters
% are constants to be fitted to data to get correct force profile. 

alpha = 1;
beta=1;
gamma=.005;
delta = .05;
inhib_threshhold=1.0;
iters = 10000;
dt = .01;
duty = 1;
reps = 4;
period = iters/reps;

% Read in parameters

if(isfield(pars,'alpha')) 
    alpha = pars.alpha;
end

if(isfield(pars,'beta')) 
   beta = pars.beta;
end

if(isfield(pars,'gamma')) 
    gamma = pars.gamma;
end

if(isfield(pars,'delta')) 
    delta = pars.delta;
end

if(isfield(pars,'inhib_threshhold')) 
    inhib_threshhold= pars.inhib_threshhold;
end

if(isfield(pars,'iters')) 
    iters= pars.iters;
end

if(isfield(pars, 'dt'))
    dt = pars.dt;
end

if(isfield(pars, 'duty'))
    duty = pars.duty;
end

if(isfield(pars, 'reps'))
    reps = pars.reps;
end

if(isfield(pars, 'period'))
    period = pars.period;
end

% Initialize all arrays we iterate over

forceprofile = zeros(1,iters);
M = zeros(1,iters); % number of bleached molecules
Mdots = zeros(1,iters); % rate of bleaching, force is proportional to Mdot
I = zeros(1,iters); % resistance of molecules to UV excitation
Idots = zeros(1,iters); % rate of resistance change

% Generate the intensity profile. Dimensionless square wave UV, ranges from 0 to 1 in strength.
UVstimulus = UVstimulus_gen(duty, period, reps);

j=2; 
while j<= iters 
    % iterate from t \in [0, dt*iters]
    dI = (alpha*UVstimulus(j-1) - beta*I(j-1))*dt;
    dM = (gamma*UVstimulus(j-1)*heaviside(inhib_threshhold - I(j-1)) - delta*M(j-1))*dt;
    
    Idots(j) = dI/dt;
    I(j) = I(j-1)+dI;
    
    Mdots(j) = dM/dt;
    M(j) = M(j-1)+dM;
    
    j=j+1;
end

% The force is proportional to \dot M, but always positive. 

forceprofile = Mdots;
forceprofile(forceprofile<0) = 0;
end