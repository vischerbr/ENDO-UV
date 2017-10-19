% #NOTES: currently changing wrong k. Need to alter k of barrier, not k of
% #second well. Moreover, first pass time should only depend on k1k2
% #product, so i can probably just iterate over one and fix the other. check
% #to make sure that this works by fixing the first and iterating the other
% #direction

% #new note: spline is almost certainly incorrect. arbitrary values are
% #causing problems. probably a good idea to run this stuff by bo, maybe he
% #can help fix the spline
%
% use a phi - 4 potential to model double well - ax^2  + bx^4, average
% first pass time between multiple runs. look up stochastic resonance -
% driving in a double well potential causes odd features in the first
% passage times! Try to simulate a driven stocastic system, check how the
% driving frequency changes the mean passage time, or check for periodic
% motion of the particle. 

% note - for double well, both curvature constants are controlled by the
% same paramater, the scaling on the the quadratic term. implementation NOT
% complete in loops, beware

pars.nothing = 1;
pars.kbt = .8;
pars.gamma = 1;
pars.iters = 100000;
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


UVprofile = forceprofile_gen(pars);

[vcaltotal,firstpasstimes,potential] = model_virtualdatagen_doublewell(UVprofile,pars);
%meanvcal = mean(vcaltotal);
%plot(meanvcal);
%plot(firstpasstimes);
%figure,fnplt(potential);
    
%hold on;
%plot(vcal,0.2*(1:length(vcal)));
%xlim([0,max(vcal)]);