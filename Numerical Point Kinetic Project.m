% Charles Goodman
% 20 November 2018
% NE 301
% Project 2: Numerical Point Kinetics
clear; close; clc;

%% Obtain input values
maxErr = .0001;
finalval = -1;
fname = input('File name of test case?: ','s');
data = dlmread([fname,'.csv']);
num = data(1,1);
li = log(2)./data(2:num+1,1);
Bi = data(2:num+1,2);
L = data(num+2,1);
rho = data(num+3,1);
tf = data(num+4,1);
steps = data(num+5,1);

%% Setup varibles
n = zeros(1,steps+1); % setup matrix for n over time
    
B = sum(Bi); % calculate beta
rho = rho * B; % convert from dollars to unitless


while abs(n(end)-finalval) > maxErr
    finalval = n(end);
    
    %% Numerical Solve
    
    n = zeros(1,steps+1); % setup matrix for n over time
    Ci = zeros(num,steps+1); % setup matrix for Ci over time
    n(1) = 1; % set initial value
    Ci(:,1) = Bi./li./L.*n(1);

    delt = tf/steps; % time step
    
    for s = 1:steps
        n(s+1) = (1+(rho-B)/L*delt).*n(s) + sum(li.*Ci(:,s))*delt;
        Ci(:,s+1) = Ci(:,s) + Bi.*n(s)./L.*delt - li.*Ci(:,s).*delt;
    end

    %% Output Data file
    out = [0:delt:tf;n;Ci];
    csvwrite(['soln_',fname,'_',num2str(steps),'.csv'],out);
    
    steps = steps*10;

end

