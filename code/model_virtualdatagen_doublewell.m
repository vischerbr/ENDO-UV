% This is the main simulation. Here, we setup a dynamical simulation in the
% overdamped regime so that \dot x \propto F. We simulate for a fixed
% number of iterations, repeating this simulation for each sampling. Our
% desired calculations are the particle trajectories under these forces and
% the times at which each particle first passes the barrier.
function [vcaltotal,firstpasstimes,potentialspline] = model_virtualdatagen_doublewell(UVprofile,pars)
%UVprofile in seconds



maxruns = 1;
samplegoal =10;

totalfirstpass = [];

dt = 0.01;
framerate = 100;
kbt = 1.5;

% Drag coeffecient. Appears to have little impact on 
% simulation, kept for dimensional reasoning.
damping = 1; 
iters=100000;

% Read in paramters.

% Thermal energy.
if(isfield(pars,'kbt'))
    kbt = pars.kbt;
end

% Position sampling frequency, unitless.
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

totalsimlen = iters;%unitless

% Output arrays
vcaltotal = [];
firstpasstimes = [];

% Construct the potential profile and the force due to the double well.
potentialspline = doublewellpotential_gen(peakheight,secondwellbottom,k1,k2,k3,xmin1,xmin2,xmax);
cellforcespline = fnder(potentialspline);

samplesize = 0;

% Do multiple trials for each iteration, average
while samplesize < samplegoal

    % Initialization, 
    tnow = 0;
    %initial position centered at minimum of left well
    xnow = -xmin1;

    % track framerate
    tempcount = 0;

    %

    %the effective time of the force
    %tforce = 0;

    % boolean first pass vcaltotacheck 
    firstpass = 0;

    % initialize trajectory array
    vcal = [];
    % do the simulation
    i=1;
    while i<=iters

        UVforce = UVprofile(i);
%             if UVforce>1e-2
%                 tforce = tforce+UVforce*tforcegrowrate*dt;
%             else
%                 tforce = tforce-tforcedecayrate*dt;
%                 tforce = max([tforce,0]);
%             end
% 
% 
%             if tforce > tdecay
%                 decayfac = exp((tforce-tdecay)/tdecaytimescale);
%                 UVforce = UVforce/decayfac;
% 
%             end

        %approximately 2*ksecondwell(xnow-1) within second well
        cellforce = fnval(cellforcespline,xnow);

        %something like Ornstein-Uhlenbeck process
        dx = (UVforce-cellforce)*(dt/damping)+sqrt(2*kbt)*randn(1,1)*sqrt(dt/damping);

         %activate force decay%
    %    if((xnow<=Xb) && ((xnow+dx)>Xb))    
    %        decayflag = 1; 
    %    end

         %reset force decay
    %    if((xnow>=Xa) && ((xnow+dx)<Xa))
    %        decayflag = 0;
    %        tdecay = 0;
    %    end 

        % check for first pass

        if xnow> xmax && firstpass == 0
            firstpass = 1;
            firstpasstimes = [firstpasstimes, i*dt];
        end

        xnow = xnow+dx;
        i = i+1;
        tempcount = tempcount+1;

        if(tempcount == framerate)
            vcal = [vcal,xnow];
            tempcount = 0;

        end

    end
    if firstpass == 0
        firstpasstimes = [firstpasstimes, 0]
    end
    samplesize = samplesize+1
    vcaltotal = [vcaltotal; vcal];
end
totalfirstpass = [totalfirstpass;firstpasstimes]
totalruns = totalruns+1;

end