% This function takes in parameters that define the double well. We use a
% cubic spline interpolation to both have freedom in shaping the well 
% (as opposed to, say, a \Phi-4 potential) and to make taking derivatives
% later simpler. Coefficients are not robust, a later implementation may
% make this easier to use. 
function potentialspline = doublewellpotential_gen(barrierheight,secondwellbottom,k1,k2,k3,xmin1,xmin2,xmax)

% Shaping the first well
firstwellx = -0.1*(xmax-xmin1):0.1:(xmin1 + 0.01*(xmax-xmin1));
firstwelly = k1*(firstwellx-xmin1).^2;

% Creating the barrier (or peak)
barrierx = (xmin1 +(xmax-xmin1)*0.995):0.1:(xmax+ (xmin2-xmax)*0.4);
barriery = -k2*(barrierx-xmax).^2+barrierheight;

% Shaping the second well.
secondwellx = (xmax+(xmin2-xmax)*0.6):0.1:1.4*xmin2;
secondwelly = k3*(secondwellx-xmin2).^2+secondwellbottom;

% Use cubic spline to meld each piece together.
allx = [firstwellx,barrierx,secondwellx];
ally = [firstwelly,barriery,secondwelly];
potentialspline = csapi(allx,ally);
end