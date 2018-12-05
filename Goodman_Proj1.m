% Charles Goodman
% 21 September 2018
% NE 301
% Project 1: Numerical Decay Chain
clear; close; clc;

%% Obtain input values
num = input('Number of Nuclides? ');
T = zeros(num,1);
N = zeros(num,1);
for n = 1:num
    T(n) = input(['Half-life of Nuclide ', num2str(n), '? ']);
end
for n = 1:num
    N(n) = input(['Initial value of Nuclide ', num2str(n), '? (atoms) ']);
end
t0 = input('Starting time? ');
tf = input('Ending time? ');
steps = input('Number of time steps? ');
A0 = input('Source strength of the first isotope? ');
A0t = input('Duration of source? ');

%% Setup varibles
l = log(2) ./ T; % decay constants
delt = (tf-t0)/steps; % time step
fprintf("max lambda*deltaT = %f",max(delt*T));

%setup source vector
%syms t
%sources = symfun([A0*heaviside(A0t-t);zeros(num-1,1)],t);


%% Numerical Solve
for s = 1:steps
    t = t0 + delt*s;
    N(:,end+1) = N(:,end) - l.*N(:,end)*delt + [0; l(1:end-1).*N(1:end-1,end)*delt];
    if t <= A0t
        N(1,end) = N(1,end) + A0*delt;
    end
end

%% Output Data file
name = ['dat',num2str(num),'_',num2str(steps),'.csv'];
out = [t0:delt:tf;N];
csvwrite(name,out);

%two nuclide check error
%A5 = 1.39761082;
%B5 = (2^-2.5-2^-5)/log(2)+2/log(2)*(1-2^-2.5);
%err = N(:,1+ceil(5/delt))-[A5;B5];
%max(abs(err))

%format long;
%disp(N(:,1+ceil(30/delt)));

