function f = StatusCheck(P)
% 
%       f = StatusCheck(P)
% 
% Read which status we are in (i.e. 'shock', 'rarefaction') and calculates
% accordingly the correct function for P* and u*.
%


global gamma N u LeftStatus RightStatus counter

% Constants
A_l = 2/((gamma+1)*u(1)); % A_l = 2/((gamma+1)*rho_l)
A_r = 2/((gamma+1)*u(4)); % As above
B_l = ((gamma-1)/(gamma+1))*u(2); % B_l = mu * P_l
B_r = ((gamma-1)/(gamma+1))*u(5); % as above

% Speed terms
a_l = (gamma*u(2)/u(1))^0.5; % Sound speed right
a_r = (gamma*u(5)/u(4))^0.5; % Sound speed right

du = u(6)-u(3);

% Functions
f_l =  zeros(N, 1);
f_r =  zeros(N, 1);


% Checking...
switch LeftStatus
    case 'shock'
        f_l = (P-u(2))*(sqrt((A_l)/(P+B_l)));
    case 'rarefaction'
        f_l = (2*a_l)/(gamma-1)*((P/u(2))^((gamma-1)/(2*gamma))-1);
end

switch RightStatus
    case 'shock'
        f_r = (P-u(5))*(sqrt((A_r)/(P+B_r)));
    case 'rarefaction'
        f_r = (2*a_r)/(gamma-1)*((P/u(5))^((gamma-1)/2*gamma)-1);
end

if(counter)
    f = 0.5*(u(3)+u(6))+0.5*(f_r - f_l); % In this case we want to calculate u_star
else
    f = f_l + f_r + du; % If we want to find P_star
end