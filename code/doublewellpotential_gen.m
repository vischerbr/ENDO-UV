%function potentialspline = doublewellpotential_gen(peakheight,curvnow)
%function potentialspline = doublewellpotential_gen(peakheight,secondwellbottom,kfirstwell,ksecondwell)
function potentialspline = doublewellpotential_gen(peakheight,secondwellbottom,k1,k2,k3,xmin1,xmin2,xmax)

% phi-4 potential
% xmin = 2*sqrt(peakheight)/(curvnow);
% xs = -1.5*xmin:0.1/curvnow:1.5*xmin;
% 
% potentialspline = csapi(xs,-0.5*curvnow^2*xs.^2 + 0.25*(curvnow^4/(4*peakheight))*xs.^4 + peakheight);


% full control of potential %
% first well
firstwellx = -0.1*(xmax-xmin1):0.1:(xmin1 + 0.01*(xmax-xmin1));
firstwelly = k1*(firstwellx-xmin1).^2;

%peak
peakx = (xmin1 +(xmax-xmin1)*0.995):0.1:(xmax+ (xmin2-xmax)*0.4);
peaky = -k2*(peakx-xmax).^2+peakheight;

%second well
secondwellx = (xmax+(xmin2-xmax)*0.6):0.1:1.4*xmin2;
secondwelly = k3*(secondwellx-xmin2).^2+secondwellbottom;

allx = [firstwellx,peakx,secondwellx];
ally = [firstwelly,peaky,secondwelly];
potentialspline = csapi(allx,ally);

% original implementation
% firstwellx = -1:0.1:0.1;
% firstwelly = kfirstwell*firstwellx.^2;
% %peak
% 
% peakx = 0.3:0.02:0.7;%0.3:0.02:0.6;
% peaky = -ksecondwell*(peakx-0.5).^2+peakheight;%-40*(peakx-0.5).^2+peakheight;
% %second well
% 
% secondwellx = 1.0:0.05:10;%0.9:0.05:10;
% secondwelly = 5*(secondwellx-1).^2+secondwellbottom;
% allx = [firstwellx,peakx,secondwellx];
% ally = [firstwelly,peaky,secondwelly];
% potentialspline = csapi(allx,ally);


end