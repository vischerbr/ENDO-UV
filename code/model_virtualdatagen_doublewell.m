function [vcaltotal,firstpasstimes,potentialspline] = model_virtualdatagen_doublewell(UVprofile,pars)
%UVprofile in seconds



% setup iteration over vals
totalruns = 1;
maxruns = 1;
samplegoal =10;

totalfirstpass = [];

dt = 0.01;
framerate = 100;
kbt = 1.5;
damping = .01; %drag coeff
omega = 0.5;
iters=100000;


if(isfield(pars,'kbt'))
    kbt = pars.kbt;
end

%simulated trajectory will only be reported every framerate-th of time steps

if(isfield(pars,'framerate')) 
    framerate = pars.framerate;
end

if(isfield(pars,'dt'))
    dt = pars.dt;
end

if(isfield(pars,'iters'))
    iters = pars.iters;
end

totalsimlen = iters;%unitless


peakheight = 3.0;
secondwellbottom = 2.5;
k1 = .03;
k2 = .02;
k3 = .0025;
xmin1 = 0;
xmax = 70;
xmin2 = 80;
% kfirstwell = 20; 
% ksecondwell = 20;

if(isfield(pars,'peakheight'))%peakheight of potential
    peakheight = pars.peakheight;
end


if(isfield(pars,'curvature'))%curviness of double well
    curvature = pars.curvature;
end

% if(isfield(pars,'secondwellbottom'))%height of second well bottom
%     secondwellbottom = pars.secondwellbottom;
% end
% 

% if(isfield(pars,'kfirstwell'))%slope of first well 
%     kfirstwell = pars.kfirstwell;
% end
% dx = (UVforce-cellforce)*(dt/damping)+sqrt(2*kbt)*randn(1,1)*sqrt(dt/damping);
% if(isfield(pars,'ksecondwell'))%slope of second well 
%     ksecondwell = pars.ksecondwenowll;
%end


% tdecay = 600; % after tforce>tdecay, force will exponentially decay
% tdecaytimescale = 100; %in seconds
% tforcegrowrate = 1; %the growth rate of tforce
% tforcedecayrate = 0.5; % the decrease rate of tforce
% 
% 
% if(isfield(pars,'tdecay'))
%     tdecay = pars.tdecay;
% end
% 
% if(isfield(pars,'tdecaytimescale'))%decay rate of the force
%     tdecaytimescale = pars.tdecaytimescale;
% end
% 
% if(isfield(pars,'tforcegrowrate'))
%     tforcegrowrate = pars.tforcegrowrate;
% end
% 
% if(isfield(pars,'tforcedecayrate'))
%     tforcedecayrate = pars.tforcedecayrate;
% end



vcaltotal = [];


%iterate over kvals, dampingvals, or omegavals
while totalruns<maxruns+1
    
    %keep track of the run
    %totalruns
    
    samplesize = 0;
    firstpasstimes = [];
    
    % change the parameter value based on iteration number
    %dampingnow = damping*totalruns;
    %curvnow = curvature*totalruns;
    %omeganow = omega*totalruns;
    %forcegrownow = tforcegrowrate*totalruns
    
    % position of the well bottom
    %xmin = 2*sqrt(peakheight)/(curvnow);
    
    %potentialspline = doublewellpotential_gen(peakheight,curvnow);
    potentialspline = doublewellpotential_gen(peakheight,secondwellbottom,k1,k2,k3,xmin1,xmin2,xmax);
    cellforcespline = fnder(potentialspline);
    
    % do multiple trials for each iteration, average
    while samplesize < samplegoal
        
        %initialization, equipartition theorem
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
end