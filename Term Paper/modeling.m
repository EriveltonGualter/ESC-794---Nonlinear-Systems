% Erivelton Gualter dos Santos
%
% Ri, S., Huang, J., Wang, Y., Kim, M., & An, S. (2014). 
%    Terminal sliding mode control of mobile wheeled inverted pendulum 
%    system with nonlinear disturbance observer. Mathematical Problems 
%    in Engineering, 2014.

clear all
close all
clc

syms r l THb THw alpha

xb = l*sin(THb) + r*THw*cos(alpha);
yb = l*cos(THb) + r*THw*sin(alpha);

xw = r*THw*cos(alpha);
yw = r*THw*sin(alpha);




