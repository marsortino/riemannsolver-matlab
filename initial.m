function u = initial(u, optini)
%
%           u = initial(u, optini)
%
% Gives the different initial condition for the Riemann's problem.
% There are 5 different options:
% Cases: 
%        'case_1' : Left Rarefaction - Right Shock - Sod Conditions
%        'case_2' : Left Rarefaction - Right Rarefaction - 123 Problem ***
%        'case_3' : Left Rarefaction - Right Shock - Woodward Colella blast
%        'case_4' : Left Shock - Right Rarefaction - Woodward Colella blast*
%        'case_5' : Left shock - Right Shock - Mix case_4 & case_3
%
% *** It doesn't currently work with fzero.
% * It doesn't currently work.

switch optini
    case 'case_1'  % L rw - R sh
        u(1) = 1.0; % rho_l   / left density 
        u(2) = 1.0; % P_l     / left pressure
        u(3) = 0.0; % u_l     / left velocity

        u(4) = 0.125; % rho_r / right density
        u(5) = 0.1; % P_r     / right pressure
        u(6) = 0.0; % u_r     / right velocity

        u(7) = 0.25; % T_fin

    case 'case_2'  % L rw - R rw
        u(1) = 1.0; % rho_l   / left density 
        u(2) = 0.4; % P_l     / left pressure
        u(3) = -2.0; % u_l     / left velocity

        u(4) = 1.0; % rho_r / right density
        u(5) = 0.4; % P_r     / right pressure
        u(6) = 2.0; % u_r     / right velocity

        u(7) = 0.15; % T_fin
        
    case 'case_3'  % L rw - R sh
        u(1) = 1.0; % rho_l   / left density
        u(2) = 1000.0; % P_l  / left pressure
        u(3) = 0.0; % u_l     / left velocity

        u(4) = 1.0; % rho_r   / right density
        u(5) = 0.01; % P_r    / right pressure
        u(6) = 0.0; % u_r     / right velocity

        u(7) = 0.012; % T_fin
        
    case 'case_4'
        u(1) = 1.0; % rho_l   / left density 
        u(2) = 0.01; % P_l    / left pressure
        u(3) = 0.0; % u_l     / left velocity

        u(4) = 1.0; % rho_r   / right density
        u(5) = 100.0; % P_r   / right pressure
        u(6) = 0.0; % u_r     / right velocity

        u(7) = 0.0035; % T_fin
        
    case 'case_5'
        u(1) = 5.99924; % rho_l   / left density 
        u(2) = 460.894; % P_l     / left pressure
        u(3) = 19.5975; % u_l     / left velocity

        u(4) = 5.99242; % rho_r   / right density
        u(5) = 46.0950; % P_r     / right pressure
        u(6) = -6.19633; % u_r    / right velocity

        u(7) = 0.035; % T_fin
        
end