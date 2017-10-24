% This is the main simulation. Here, we setup a dynamical simulation in the
% overdamped regime so that \dot x \propto F. We simulate for a fixed
% number of iterations, repeating this simulation for each sampling. Our
% desired calculations are the particle trajectories under these forces and
% the times at which each particle first passes the barrier.
function [vcaltotal,firstpasstimes,potentialspline] = model_virtualdatagen_doublewell(UVprofile,pars)

samples=2;
dt = 0.01;
framerate = 100;
kbt = 1.5;

% Drag coeffecient. Appears to have little impact on 
% simulation, kept for dimensional reasoning.
damping = 1.0; 
iters=100000;

% Read in paramters.

% Thermal energy.
if(isfield(pars,'kbt'))
    kbt = pars.kbt;
end

% Position sampling frequency, unitless. Should be integer multiple of
% iterations.
if(isfield(pars,'framerate')) 
    framerate = pars.framerate;
end

% Time step, seconds.
if(isfield(pars,'dt'))
    dt = pars.dt;
end

% Number of iterations, unitless.
if(isfield(pars,'iters'))
    iters = pars.iters;
end

% Number of total simulations, unitless.
if(isfield(pars,'samples'))
    samples = pars.samples;
end


% These are the best fitting parameters so far. May require further
% optimizing.

peakheight = 3.0;
secondwellbottom = 2.5;
k1 = .03;
k2 = .02;
k3 = .0025;
xmin1 = 0;
xmax = 70;
xmin2 = 80;

% Output arrays
vcaltotal = zeros(samples, iters/framerate); % spatial sampling should be position vs. time
firstpasstimes = zeros(1, samples);

% Construct the potential profile and the force due to the double well.
potentialspline = doublewellpotential_gen(peakheight,secondwellbottom,k1,k2,k3,xmin1,xmin2,xmax);
cellforcespline = fnder(potentialspline);

samplesize = 1;

% Do multiple trials for each iteration, average
while samplesize <= samples

    % Initialization, 
    tnow = 0;
    %initial position centered at minimum of left well
    xnow = -xmin1;

    % track framerate
    tempcount = 0;

  
    % boolean first pass vcaltotacheck 
    firstpass = 0;

    % initialize trajectory array
    vcal = zeros(1,iters/framerate);

    i=1;
    while i<=iters

        % Three forces act on the particle: the force from the incident UV
        % light; the force from the doublewell potential; and thermal
        % noise.
        UVforce = UVprofile(i);
        cellforce = fnval(cellforcespline,xnow);

        %something like Ornstein-Uhlenbeck process
        dx = (UVforce-cellforce)*(dt/damping)+sqrt(2*kbt)*randn(1,1)*sqrt(dt/damping);

        % Check to see if the particle has passed the barrier yet

        if xnow> xmax && firstpass == 0
            firstpass = 1;
            firstpasstimes(samplesize) =  i*dt;
        end

        xnow = xnow+dx;
        i = i+1;
        tempcount = tempcount+1;

        % Store trajectories every framerate of steps
        if(tempcount == framerate)
            vcal(floor(i/framerate)) = xnow;
            tempcount = 0;
        end

    end
    
    % If the particle failed to get over the barrier during the simulation,
    % store something noticeable.
    if firstpass == 0
        firstpasstimes(samplesize) = -1;
    end
    vcaltotal(samplesize,:) = vcal;
    samplesize = samplesize+1
end
end