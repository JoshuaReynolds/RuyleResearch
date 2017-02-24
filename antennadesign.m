% Alexander Moreno
% 3-18-2015
%
% Computes dimenesions for plannar antenna design for given user inputs
%
% Input parameters
% er: dieletric constant
% f: frequency
% Z0: impedance
% maxL: maximum length of the antenna
% maxW: maximum width of the antenna
% h: height of the substrate
% user: array of values
%   [restr, wslb, seplb, wsub, sepub]
%   restr: 1 or 0[boolean]. 
%   1 indicates user will input values:
%   wslb, seplb, wsub, sepub and 0 indicates user will use default values
%   0 indicates user will not input values
%
%
% Ouput parameters
%
% error: array of values, intial value is 0
%   [0, 1, 2]
%   indciates the error that has occur within the program
%    0: no error(s)
%    1: spiral dimensions exceed user's max dimensions
%    2: VSWR of serrations exceeds recommend value, change max dimensions
%
% spiralDim: outputs dimensions require for spiral
%   array of values
%   [ws1,sep1,N1,fval1,thetastopcurve1,spiralwidth1,lengthlo1,ws2,sep2,N2,fval2,thetastopcurve2,spiralwidth2, lengthlo2]
%
%   ws1: width of spiral for first spiral
%   sep1: space between loops within first spiral
%   N1: number of turns for first spiral
%   fval1:
%   thetastopcurve1: angle of spiral [Degrees] for first spiral
%   spiralwidth1: total width of first spiral
%   lengthlo1: length from slot to spiral width becomes constant for first
%   spiral
%   ws2: width of spiral for second spiral
%   sep2: space between loops within second spiral
%   N2: number of turns for second spiral
%   fval2:
%   thetastopcurve2:  angle of spiral [Degrees] for second spiral
%   spiralwidth2: total width of second spiral
%   lengthlo2:length from slot to spiral width becomes constant for second
%   spiral
%
% dstandoff: array of values
%   [dstandoffl,dstandoffw]    
%   dstandoffl: scaled spacing between spiral and serrations  
%   dstandoffw: scaled spacing between spiral and antenna edge
%
% serrDim: output dimensions require for serrations
%   array of values
%   [VSWR,Wserr,Dserr]
%   [VSWR, Wserr, Dserr, Zc, Zin]
%   VSWR: VSWR between input impedance of serration and ideal
%   Wserr: width of middle serration(s)
%   Dserr: depth of middle serration(s)
%
% antennaDim: 
%
% La: length of aperture(slot)
%
% wa: width of aperture(slot)
%
% VSWR: VSWR of ideal input impedance and calucated input impedance of the
% parallel plate and serrations

function [error, spiralDim, dstandoff,  Wserr, Dserr, antennaWidth, antennaLength, wa, La, VSWR] = antennadesign(er,f,Z0,maxL,maxW,h,user)
% intial output values 
% so Matlab does not crash when leaving out of function earlt
error     = 0;
spiralDim = 0;
dstandoff = 0;
Wserr = 0;
Dserr = 0;
antennaWidth=0; 
antennaLength=0; 
wa = 0;
La = 0;
VSWR = -123;

% user input values or calucated from user's input values
c = 3e8;          %[m/s] speed of light
lambda0 = c/f;
mu0= 1.256637e-6; % mu0 
e0 = 8.8542e-12;  % es0

% Following parameters @ f=433[MHz] are inital values
fi = 433e6;         %[Hz]
lambda0i = c/fi;
La_int   = 15.28e-3;%[m] length of slot
wa_intial= 4e-3;    %[m] width of slot
wslbi  = 0.3e-3;    %[m] lower bound width of spiral
seplbi = 0.1e-3;    %[m] lower bound spacing of loop(s) in spiral
wsubi  = 2e-3;      %[m] upper bound width of spiral
sepubi = 1.5e-3;    %[m] upper bound spacing of loop(s) in spiral
ws0i   = wsubi/2;   %[m] intial width of spiral
sep0i  = sepubi/2;  %[m] inital spacing of loop(s) in spiral
dstandoffl = 7.13e-3; %[m] space from spiral to edge of pp
dstandoffw = 1.33e-3; %[m] space from spiral to serrations

% Scaling intial values 
wa = (wa_intial*lambda0)/lambda0i; % [m] slot width scaled by user's input
dstandoffl = (dstandoffl*lambda0)/lambda0i;
dstandoffw = (dstandoffw*lambda0)/lambda0i;
dstandoff = [dstandoffl, dstandoffw];

% Compute effective lambda using Cohn's
[~,lambda_eff_intial,~,~,~]=slotcalcscohn(f,er,wa_intial,h); % intial lambda
[~,lambda_eff_user,~,~,~]=slotcalcscohn(f,er,wa,h);          % scaled lambda
% Scale La by user's input
La = (La_int*lambda_eff_user)/lambda_eff_intial;

% restrictions values
% lower & upper bound of N number of loops within spiral(s)
N0 = 1;
Nlb = 1;
Nub = 8;

if(user(1)==0)
    % user will be using our restrictions
    % scale to lambda0
    wslb  = (wslbi*lambda0)/lambda0i;  %[m]
    seplb = (seplbi*lambda0)/lambda0i; %[m]
    wsub  = (wsubi*lambda0)/lambda0i;  %[m]
    sepub = (sepubi*lambda0)/lambda0i; %[m]
    ws0   = wsub/2;                    %[m]
    sep0  = sepub/2;                   %[m]
else
    % user inputs restrictions
    wslb   = user(2); %[m]
    seplb  = user(3); %[m]
    wsub   = user(4); %[m]
    sepub  = user(5); %[m]
    ws0    = wsub/2;  %[m]
    sep0   = sepub/2; %[m]
end

% functions
%% Hank's code [works~]
% computes the N-turns (for ideal transfomer) between slot and ground plane
Nturns = transformer(f,er,wa,h);
%a='hank'
%% Ruyle's code [works~]
% Compute impedance(s) (left & right) of the slot
[ZL1, ZL2] = ZLRimpedanceslot(La,f,Z0,wa,er,h,maxW,Nturns);
%a='slot'
%% Ruyle's code [works~]
% Compute the dimensions for spiral
[ws1,sep1,N1,fval1,thetastopcurve1,spiralwidth1, lengthlo1]=findoptimumspiraldimstomatchimp(ws0,sep0,N0,wslb,seplb,Nlb,wsub,sepub,Nub,er,h,f,ZL1,ZL2,Z0,wa);
[ws2,sep2,N2,fval2,thetastopcurve2,spiralwidth2, lengthlo2]=findoptimumspiraldimstomatchimp(ws0,sep0,N0,wslb,seplb,Nlb,wsub,sepub,Nub,er,h,f,ZL2,ZL1,Z0,wa);
spiralDim = [ws1,sep1,N1,fval1,thetastopcurve1,spiralwidth1,lengthlo1,ws2,sep2,N2,fval2,thetastopcurve2,spiralwidth2, lengthlo2];
%a='spiral'
%% Check if dimenesions are below user's maximum width and length [works~]
lpi = 2*dstandoffl + spiralwidth1 + spiralwidth2 - wa;
%wp  = 2*dstandoffw + spiralwidth1 + spiralwidth2 + lengthlo1(numel(lengthlo1)) + lengthlo2(numel(lengthlo2)) + La;
wp  = 2*dstandoffw + spiralwidth1 + spiralwidth2 + lengthlo1 + lengthlo2 + La;

if((wp<=maxW) && (lpi<maxL))   
    % dimensions are within given max dimensions
    % find Lpp and Dmax
    % do return function
    Dmax = (maxL-lpi)/2;    
    Lpp = lpi/2
    error = 0;
else
    % dimensions are larger than given max dimensions
    error = 1;
    return;
end

%% Alex's code
% compute serration's dimensions (Dserr & Wserr)
%{
lb  = [0.1*Dmax, 0.1*maxW];    % lower bounds
ub  = [Dmax, maxW];            % upper bounds
%x0   = [0.1*Dmax, 0.1*maxW];    % stating point
Din = abs(ub(1)-lb(1))/2; 
Win = abs(ub(2)-lb(2))/2; 
x0  =  [Din, Win];
%}

% lb: lower bounds
% ub: upper boounds
%   [Serration Depth, Serration Width]
lb  = [0.1*Dmax, 0.1*wp]; % lower bounds
ub  = [Dmax, wp];         % upper bounds
Din = abs(ub(1)-lb(1))/2;   % starting point, Dserr
Win = abs(ub(2)-lb(2))/2;   % starting point, Wserr          
x0  = [Din, Win];          % starting points
antennaDim = [lpi, wp];
a='serr'
% older ver
%[VSWR, Wserr, Dserr] = findoptimalserration(f,h,er,Lpp,maxL,maxW,x0,lb,ub);
%[VSWR, Wserr, Dserr, Zc, Zin] = findoptimalserration(f,h,er,Lpp,maxL,maxW,x0,lb,ub);
[VSWR, Wserr, Dserr, Zc, Zin] = findoptimalserration(f,h,er,Lpp,maxL,wp,x0,lb,ub);
antennaLength = lpi+2*Dserr;
antennaWidth  = wp;
% check if VSWR is within bounds
if (VSWR <= 7)
    % VSWR is within bound
    error = 0;
else
    % invalid VSWR
    error = 2;
    return;
end

end