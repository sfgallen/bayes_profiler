function [Z0] = river_incision_forward_model(S,S_DA,S_U,Ui,K,m,n,run_time,dt)

% This function used the implicit finite difference scheme to solve the one
% dimensional stream power model (E = K*A^m*S^n) on a TopoToolbox STREAMobj
% using the  scheme of Braun and Willet (2013). 
% Some of the codes are modified and simplified versions of the solvers 
% available with the TTLEM (Campforts et al., 2017, E-Surf). 
% 
% Inputs:
% (1) S - Matlab structure formated the same as a STREAMobj from TopoToolbox
% (2) S_DA - Drainage area vector mapped to STREAMobj
% (3) S_U - the final uplift m/yr [Can be a scalar or vector]
% (4) Ui - the intial uplift rate in m/yr [Can be scalar or vector]
% (5) K - the erodibility constant [Can be scalar or vector]
% (6) m - drainage area exponent (scalar)
% (7) n - slope exponent (scalar)
% (8) run_time - length of model run time in years (scalar)
% (9) dt - time step used in the model in years (scalar)
%
% Output:
% (1) Z0 - topologically sorted river network elevations at the end of the
% model run
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 09/28/2020

% Run the fastscape data preperation function
[d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, K, m);

% Calculate Chi for the river network
mn = m/n;

% Make sure uplift doesn't go negative (note, you can comment this out if
% you don't care about this condition
S_U(S_U <= 0) = 1e-7;

% calculate initial steady-state river elevations
S_Zi = calculate_z(S,S_DA,Ui,K,mn,n);

% get time loop prepped
tsteps = round(run_time/dt);

% allocate memory before running the model and calc. initial erosion rate
% and sediment flux
Z0 = S_Zi;

% run the forward model
for t = 1:tsteps
    
    % update uplift field as needed
    Z1 = fastscape_eroder_outlets(Z0, n, dt, A, d, r, dx, S_U,outlet_nodes);
        
    % make sure river doesn't flow backwards [note you can comment this out
    % if you don't care about this condition]
    Z1 = check_z(S,Z1);

    % update Z0 data with new elevations
    Z0 = Z1;
    
end
end

function z = calculate_z(S,S_DA,S_U,K,mn,n)
% function calculate the steady state elevation of the river network
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 09/28/2020

    z = zeros(size(S.distance));
    Six = S.ix;                         % donors
    Sixc = S.ixc;                       % recievers
    Sd = S.distance;                    % distance from mouth
    Sa = (S_U./K).^(1/n).*(1./(S_DA)).^mn;       % chi transformation variable
    
    % calculating selevation for the entire river network
    for lp = numel(Six):-1:1
        z(Six(lp)) = z(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sd(Sixc(lp))-Sd(Six(lp))));
    end
end

function z = check_z(S,z)
% this function makes sure the river doesn't flow backwards
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 09/28/2020
    Six = S.ix;                         % donors
    Sixc = S.ixc;                       % recievers
    
    % ma
    for lp = numel(Six):-1:1
        z_pre = z(Sixc(lp));
        z_cur = z(Six(lp));
        if z_cur <= z_pre
            z(Six(lp)) = z_pre+0.1;
        end
    end
end

function [d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, S_K, m)

% simplified from Campforts and Schwanghart TTLEM updateDrainDir.m function
% to handel incision along STREAMobj only
%
% Inputs:
%   S               STREAMobj
%   S_DA            drainage area in meters^2 along the stream network
%   S_K             A scalar or vector of the erodibility coefficent
%   m               The drainage area exponent
%
% Outputs:
%   d               donors
%   r               recievers
%   A               the "velocity field" of the steam power incision model
%   dx              distance between nodes along network
%   outlet_nodes	locations of outlets
%
% Modified by: Sean F. Gallen

% get donor and recievers from STREAMobj
d = S.ix;
r = S.ixc;

% get velocity field
A = S_K.*S_DA.^m;

% Find outlet nodes in STREAMobj
outlet_nodes = S.outlets;
for i = 1:length(outlet_nodes)
    outlet_nodes(i) = find(S.IXgrid == outlet_nodes(i));
end

% get distance between nodes along profile
Sd = S.distance;
dx = abs(Sd(d) - Sd(r));

end

function S_Z = fastscape_eroder_outlets(S_Z, n, dt, A, d, r, dx, U, outlets)

% simplified from Campforts and Schwanghart TTLEM funerosion_implin.m and
% funerosion_impnlin.m to model incision along STREAMobj only using the
% detachment limited stream power incision model. There is a data
% preperation frunction, fastscape_eroder_data_prep.m that should be run
% before one goes into the time evolution forloop where this function is
% applied. Modification also fixes all outlet elevations
%
% Inputs:
%   S_Z     River network elevations in meters sorted as STREAMobj
%   n       Slope exponent in stream power incision model
%   dt      time step in years
%   A       the "velocity field" of the steam power incision model
%   d       STEAMobj donors
%   r       STREAMobj recievers
%   dx      distance between nodes in STREAMobj
%   U       scalar or vector of uplift rate in meters per year
%   outlets outlet positions along the stream network
%
% Outputs
%   S_Z     Updated river network elevations in meters sorted as STREAMobj
%
% Modified by: Sean F. Gallen

time=dt;
dte = dt;

while time>0
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    
    S_Z=S_Z+dte.*U; % add uplift to elevation field
    if isscalar(U)
        S_Z(outlets) = S_Z(outlets)-dte.*U; %except for the outlets
    else
        S_Z(outlets) = S_Z(outlets)-dte.*U(outlets);
    end
    
    if n == 1
        for j = numel(d):-1:1
            tt      = A(d(j))*dte/(dx(j));
            S_Z(d(j)) = (S_Z(d(j)) + S_Z(r(j))*tt)./(1+tt);
        end
    else
        for j = numel(d):-1:1
            
            tt      = A(d(j))*dte/(dx(j));
            % z_t
            zt      = S_Z(d(j));
            % z_(t+dt) of downstream neighbor
            ztp1d   = S_Z(r(j));
            % dx
            dx_n      = dx(j);
            
            % initial value for finding root
            if ztp1d < zt
                ztp1    = newtonraphson(zt,ztp1d,dx_n,tt,n);
            else
                ztp1    = zt;
            end
            
            if ~isreal(ztp1) || isnan(ztp1)
               % disp('Non real solutions converted to real')
                ztp1=real(ztp1);
            end
            S_Z(d(j))=ztp1;
        end
    end
end
    function ztp1 = newtonraphson(zt,ztp1d,dx,tt,n)
        
        tempz   = zt;
        tol = inf;
        
        while tol > 1e-3
            % iteratively approximated value
            ztp1  =  tempz - (tempz-zt + ...
                (tt*dx) * ...
                ((tempz-ztp1d)./dx)^n) / ...
                (1+n*tt*((tempz-ztp1d)./dx)^(n-1));
            tol   = abs(ztp1-tempz);
            tempz = ztp1;
        end
    end
end