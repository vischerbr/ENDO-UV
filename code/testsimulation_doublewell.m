% Main script file. This script generates a UV intensity spectrum to be
% passed to the main simulation, which dynamically simulates particles in a double well
% potential acting under thermal effects and a UV driving force.

pars.nothing = 1;
pars.kbt = .8;
pars.iters = 100000;
pars.samples = 2;
pars.framerate = 100;
pars.dt = .01;
pars.duty = .95;
pars.reps = 3;
pars.period = pars.iters/pars.reps;
pars.alpha = 1.0;
pars.beta = 1.0;
pars.gamma = 0.1;
pars.delta = .01;
pars.inhib_threshhold = 1.0;

% Generate the UV force profile

UVprofile = forceprofile_gen(pars); 

% Simulate particles in the double well under the generated force  

[vcaltotal,firstpasstimes,potential] = model_virtualdatagen_doublewell(UVprofile,pars);
%meanvcal = mean(vcaltotal);
%plot(meanvcal);
%plot(firstpasstimes);
%figure,fnplt(potential);
    
%hold on;
%plot(vcal,0.2*(1:length(vcal)));
%xlim([0,max(vcal)]);