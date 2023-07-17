function f = pstar(P)
%
%           f = pstar(P)
%
% It gives back a first guess of the function f(p) which root gives the
% pressure in the W* region.
% 
% Please note that this P* is mainly used to determine what is happening on
% the left and on the right of the discontinuity (i.e. shock or rarefation)
%

global gamma u

% Constants
A_l = 2/((gamma+1)*u(1)); % A_l = 2/((gamma+1)*rho_l)
A_r = 2/((gamma+1)*u(4)); % As above
B_l = ((gamma-1)/(gamma+1))*u(2); % B_l = mu * P_l
B_r = ((gamma-1)/(gamma+1))*u(5); % as above

% % Sound speed
% a_l = (gamma*u(2)/u(1))^0.5; % Sound speed right
% a_r = (gamma*u(5)/u(4))^0.5; % Sound speed right

% Functions
f_l = (P-u(2))*(sqrt((A_l)/(P+B_l)));
%f_l(2) = (2*a_l)/(gamma-1)*((P/u(2))^((gamma-1)/(2*gamma))-1);

f_r = (P-u(5))*(sqrt((A_r)/(P+B_r)));
%f_r(2) = (2*a_r)/(gamma-1)*((P/u(5))^((gamma-1)/2*gamma)-1);

du = u(6)-u(3);

f = f_l+f_r+du;