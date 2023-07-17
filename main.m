% Riemann's solver.
% It solves the main cases of the shock tube problem.
% 
% Bibliography: Toro, Eleuterio F., and Riemann Solvers. "Numerical Methods for Fluid Dynamics: A Practical Introduction." (1999).
% 
%   ________Left______x0______Right_____________
% Cases: 
%        'case_1' : Left Rarefaction - Right Shock - Sod Conditions
%        'case_2' : Left Rarefaction - Right Rarefaction - 123 Problem ***
%        'case_3' : Left Rarefaction - Right Shock - Woodward Colella blast
%        'case_4' : Left Shock - Right Rarefaction - Woodward Colella blast*
%        'case_5' : Left shock - Right Shock - Mix case_4 & case_3
%
%  *** It doesn't currently work with fzero. p* = 0.00189, u* = 0.00000
%
%  * fzero fails to find the correct p*. p* = 46.0950

global gamma N u LeftStatus RightStatus counter

optini = 'case_1'; % Cases: 'See Help.
N = 400; % Set by the user

% Physical conditions

gamma = 1.4;

initial_condition = zeros(7,1);
initial_condition = initial(initial_condition, optini); % It automatically get the correct initial conditions
u = initial(initial_condition, optini);

rho_l = initial_condition(1);
P_l = initial_condition(2);
u_l = initial_condition(3);
a_l = (gamma*P_l/rho_l)^0.5; % Sound speed left

rho_r = initial_condition(4);
P_r = initial_condition(5);
u_r = initial_condition(6);
a_r = (gamma*P_r/rho_r)^0.5; % Sound speed right

% Grid
x_min = 0; x_max = 1;
x0 = (x_max-x_min)/2;
T_fin = initial_condition(7);

dx = (x_max-x_min)/N;
dt = T_fin/N;
CFL =dt/dx;

if (CFL > 1) % Check if the CFL condition is satisfied
    fprintf('Error, CFL is higher than 1.\n')
    return;
end

%% Main body of the code:
% First guess at P: f(P) = 0
P_star = fzero('pstar', [0, 100000]);

% Determine which case are we in:
if (P_star > P_l)
    LeftStatus = 'shock';
elseif (P_star < P_l)
    LeftStatus = 'rarefaction';
end

if (P_star > P_r)
    RightStatus = 'shock';
elseif (P_star < P_r)
    RightStatus = 'rarefaction';
end

% Finally determining P* and u*
counter = 0;
P_star = fzero('StatusCheck', P_star);

counter = 1;
u_star = StatusCheck(P_star);

% Calculating rho* for the all different cases
rho_lstar = rho_l*(((P_star/P_l)+(gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*(P_star/P_l)+1));
rho_lfan = rho_l*(P_star/P_l)^(1/gamma); % For the rarefacted left region

rho_rstar = rho_r*(((P_star/P_r)+(gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*(P_star/P_r)+1));
rho_rfan = rho_r*(P_star/P_r)^(1/gamma); % For the rarefacted left region

% Vectors
x = x_min:dx:x_max;

data.x = x;
data.rho = zeros(N+1, 1);
data.P = zeros(N+1, 1);
data.u = zeros(N+1, 1);
data.e = zeros(N+1, 1);


% Cycle
again = 1;
t = 0;
while (t<T_fin+dt)
    if t+dt > T_fin && again
        t = T_fin - dt;
        again = 0;
    end
    
    x_middle = x0 + u_star*t;
    for i=1:N+1
        switch LeftStatus
            case 'shock'
                % Calculating main constants
                S_L = u_l - a_l*((gamma+1)/(2*gamma)*P_star/P_l + (gamma-1)/(2*gamma))^0.5; % Left shock speed
                x_ls = x0 + S_L*t; % position of Left Shock
    
                % Updating the grid
                if data.x(i) <= x_ls % Unshocked region
                    data.rho(i) = rho_l;
                    data.P(i) = P_l;
                    data.u(i) = u_l;        
                elseif (x_ls <= data.x(i) && data.x(i) <= x_middle) % Shocked region
                    data.rho(i) = rho_lstar;
                    data.P(i) = P_star;
                    data.u(i) = u_star;
                end            
            case 'rarefaction'
                % Calculating main constants
                a_lstar = a_l*(P_star/P_l)^((gamma-1)/(2*gamma)); % Sound speed behind the rarefaction 
                S_HL = u_l - a_l; % Head of rarefaction wave
                S_TL = u_star - a_lstar; % Tail of rarefaction wave
    
                x_hl = x0 + S_HL*t; % Head Position
                x_tl = x0 + S_TL*t; % Tail Position
    
                % Updating the grid
                if (data.x(i) <= x_hl) % Unperturbed Region
                    data.rho(i) = rho_l;
                    data.P(i) = P_l;
                    data.u(i) = u_l;  
                elseif (x_hl <= data.x(i) && data.x(i) <= x_tl) % Rarefacted region
                    data.rho(i) = rho_l*((2/(gamma+1)+(gamma-1)/((gamma+1)*a_l)*(u_l-(data.x(i)-x0)/t)))^(2/(gamma-1));
                    data.P(i) = P_l*((2/(gamma+1)+(gamma-1)/((gamma+1)*a_l)*(u_l-(data.x(i)-x0)/t)))^(2*gamma/(gamma-1));
                    data.u(i) = (2/(gamma+1))*(+a_l+((gamma-1)/2)*u_l+(data.x(i)-x0)/t);
                elseif (x_tl <= data.x(i) && data.x(i) <= x_middle) % Shocked region
                    data.rho(i) = rho_lfan;
                    data.P(i) = P_star;
                    data.u(i) = u_star;
                end
        end
    
        switch RightStatus
            case 'shock'
                % Calculating main constants
                S_R = u_r + a_r*((gamma+1)/(2*gamma)*P_star/P_r + (gamma-1)/(2*gamma))^0.5; % Right shock speed
                x_rs = x0 + S_R*t; % position of right shock
    
                % Updating the grid
                if (data.x(i) > x_rs) % Unperturbed Region
                    data.rho(i) = rho_r;
                    data.P(i) = P_r;
                    data.u(i) = u_r;
                elseif (x_middle <= data.x(i) && data.x(i) <= x_rs)
                    data.rho(i) = rho_rstar;
                    data.P(i) = P_star;
                    data.u(i) = u_star;
                end
            case 'rarefaction'
                % Calculating main constants
                a_rstar = a_r*(P_star/P_r)^((gamma-1)/(2*gamma)); % Sound speed behind the rarefaction 
                S_HR = u_r - a_r; % Head of rarefaction wave
                S_TR = u_star - a_rstar; % Tail of rarefaction wave
    
                x_hr = x0 + S_HR*t; % Head Position
                x_tr = x0 + S_TR*t; % Tail Position
                
                % Updating the grid
                if (x_hr >= data.x(i)) % Unperturbed Region
                    data.rho(i) = rho_r;
                    data.P(i) = P_r;
                    data.u(i) = u_r;
                elseif (x_tr <= data.x(i) && data.x(i) <= x_hr) % Fan region
                    data.rho(i) = rho_r*((2/(gamma+1)-(gamma-1)/((gamma+1)*a_r)*(u_r-(data.x(i)-x0)/t)))^(2/(gamma-1));
                    data.P(i) = P_r*((2/(gamma+1)-(gamma-1)/((gamma+1)*a_r)*(u_r-(data.x(i)-x0)/t)))^(2*gamma/(gamma-1));
                    data.u(i) = (2/(gamma+1))*(-a_r+((gamma-1)/2)*u_r+(data.x(i)-x0)/t);
                elseif (x_middle <= data.x(i) && data.x(i) <= x_tr)
                    data.rho(i) = rho_rfan;
                    data.P(i) = P_star;
                    data.u(i) = u_star;
                end
    
        end 
        % Internal Energy
        data.e(i) = data.P(i)/((gamma - 1)*data.rho(i));
    end

    % Plotting
    figure(1)

    subplot(4,1,1)
    plot(data.x,data.rho)
    title('density at time=', t)
    xlabel('x')
    ylabel('rho')

    subplot(4,1,2)
    plot(data.x, data.P)
    title('pressure')
    xlabel('x')
    ylabel('p')

    subplot(4,1,3)
    plot(data.x, data.u)
    title('velocity')    
    xlabel('x')
    ylabel('u')

    subplot(4,1,4)
    plot(data.x, data.e)
    title('internal energy')
    xlabel('x')
    ylabel('e')
    drawnow

    % In order to cut down wait time...
    fprintf('We are at %.2f %%\n', 100*t/T_fin)

    % Increment time
    t = t+dt;
    

end


