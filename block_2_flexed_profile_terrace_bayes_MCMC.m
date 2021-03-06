% This is the master script used to run the Bayesian Markov chain Monte
% Carlo (MCMC) inversion of river profiles and marine terrace-derived rock
% uplift rates assuming a two-stage block-to-flexural (1D) uplift history
% using the stream power incision model as presented in Gallen and
% Fernández-Blanco, in review, "A New Data-driven Bayesian Inversion of Fluvial 
% Topography Clarifies the Tectonic History of the Corinth Rift and Reveals 
% a Channel Steepness Threshold", JGR-Earth Surface
%
% Author: Sean F. Gallen
% contact: sean.gallen[at]colostate.edu
% date modified: 09/28/2020

% clear workspace and command window
clear
clc

% point to current directory for saving files
cur_dir = cd;

% pick a file identifier for output file naming (called fileTag here)
fileTag = 'test_model';

% define timestep (in years) of forward model used in Bayesian inversion
dt = 25000;

% define burnin and number of model interations (total interation is the
% sum of both values)
burn_in = 1e4;
n_runs = 1e5;
% Note that in Gallen and Fernández-Blanco, in review, we used a burn in of 3e5
% and 3e6 post burnin interations, which takes ~21 - 28 days to run on a
% normal computer. Above we set the burn in to 1e4 and the post-burn in
% iterations to 1e5, which will take ~16-24 hours to run on a typical
% computer. Because model initiates around the MAP solution, this should
% give the user a good idea of how the MCMC operates

%% load prepped stream data and unpack
% data was prepared using TopoToolbox (TTB). All data is topologically
% sorted
stream_data = load('stream_data_east.mat');
stream_data = stream_data.stream_data;
S = stream_data.S;              % TTB STREAMobj
Sz = stream_data.Sz;            % elevations of river network from DEM in meters
S_eud = stream_data.S_eud;      % stream node euclidian distance from fault in meters
S_DA = stream_data.S_DA ;       % upstream drainage area in m^2
Serr = 5.*ones(size(Sz));       % error on stream network elevations in meters

%% load marine terrace data from de Gelder et al. (2019)
% add in terrace data from Corinth
[dat,txt,raw] = xlsread('corinth_terrace_elevation_age_de_Gelder_2019.xlsx');

mis = dat(:,1);         % assigned Marine Isotope State (MIS)
el = dat(:,2);          % elevation of terrace inner shoreline elevaton (ISE) in meters
esd = dat(:,3);         % standard deviation of terrace ISE
age = dat(:,4);         % assigned age of terrace based on sea level correlations
asd = dat(:,5);         % standard deviation of age estimate
sl = dat (:,6);         % sea level elevation at the time of terrace fromation
ssd = dat(:,7);         % standard deviation of sea level elevation
x = dat(:,8);           % UTM easting (meters) grid zone 34N
y = dat(:,9);           % UTM northing (meters) grid zone 34N
f_dist = dat(:,10);     % Euclidian distance from approx. fault trace (see Figure 2a in G & F-B)

% make basic calculations from the terrace data
t_up = el - sl;         % Total rock uplift in meters
z_err = esd+ssd;        % uncertainty on rock uplift 

terrace_unit_age = unique(age); % ID unique terrace levels based on age

%% setup river parameters
% fold in a weighting factor for the stream data to equally weight the
% stream data and terrace data based on number of points. Note that this
% formulations stems from weighted least squares where sum(w*(x-y)^2) is
% minimized. Assuming we are minimizing sum((x-y/sig)^2), one can properly
% pack the weighting into the error (sig) by multiplying it by the square
% root of the inverse of the weight.
Ws = length(z_err)/length(Serr);
Ws =sqrt(1/Ws);
Serr = Serr.*Ws;

%% place all obervations and errors into vectors
t_obs = [Sz;t_up];          % observations
t_err = [Serr;z_err];       % observation errors

%%
% set flexural parameters and define uplift pattern
rhom_rhoc=3300 - 2700;      % difference in density between mantle and crust (kg/m^3
g = 9.81;                   % acceleration due to gravity (m/s^2)
E = 70e9;                   % Young's Modulus (kg/m/yr^2)
v = 0.25;                   % Poisson's ratio

%% set up priors. 
% Note we assume unifrom, uniformative priors and thus the values represent
% parameter ranges
Ui = [0.1, 0.5].*1e-3;              % initial uplift rate in m/yr
Uf = [0.5, 3].*1e-3;                % final uplift rate in m/yr
Kexp = [-80,-4];                    % log10(K) -- erodibility.
Te = [1000,10000];                  % effective elastic thickness in m
fault_init_time = [5e5, 1.5e6];     % fault initation time in years
n = [0.5, 20];                      % slope exponent on stream power model
m = [0.1, 20];                      % drainage area exponent on stream power model

% place priors in a matrix for to make things easier
prior_bounds = [Ui; Uf; Kexp; Te; fault_init_time; n; m];

%% set up matrices to catch relevant data before entering time MCMC loop
% initialize matrix to catch parameters accepted during the MCMC
params = nan(burn_in+n_runs,7);

% allocate memory to catch likelihoods, acceptance rate 
like_data = nan(burn_in+n_runs,3);

% make sure new matlab session use different random numbers
rng shuffle

%% set up the initial model
% generate initial model parameters [in this case they are picked about the
% approximate MAP solution based on previous model iterations
params(1,1) = 0.12e-3+0.1e-3*randn(1);      % Ui
params(1,2) = 1.65e-3+0.1e-3*randn(1);      % Uf
params(1,3) = -6.9+1e-1*randn(1);           % Kexp
params(1,4) = 2.350e3+1e2*randn(1);         % Te
params(1,5) = 6.00e5+1e4*randn(1);          % fault_init_time
params(1,6) = 7+1e-1*randn(1);              % n
params(1,7) = 2+1e-1*randn(1);              % m

%% run initial model
% cast parameters for model run (_t stands for temporary)

% initial uplift
Ui_t = params(1,1);

% final uplift at fault
Uf_t = params(1,2);

% flexed uplift 
Te_t = params(1,4);
D = E.*Te_t.^3./(12.*(1-v^2));              % flexural rigidity
lambda = (4.*D./((rhom_rhoc)*g)).^(1/4);    % flexural parameter (broken plate)
S_U = Uf_t.*exp(-S_eud./lambda).*cos(S_eud./lambda);  % flexed uplift rate pattern for stream

K_t = (10^(params(1,3))); % steam erodibility

n_t = params(1,6);      % slope exponent
m_t = params(1,7);      % DA exponent
t_int = params(1,5);    % fault initiation time

%% run the forward models
% stream incision model
[Z0] = river_incision_forward_model(S,S_DA,S_U,Ui_t,K_t,m_t,n_t,t_int,dt);
% terrace uplift model
ter_z = terrace_uplift_forward_model(age,f_dist,lambda,t_int,Ui_t,Uf_t);

% store model data in vector to compare to observed data
t_mod = [Z0;ter_z];

%% calculate relevant log-likehoods for the model parameters and results
% log-likelihood of model parameter given the priors
lpcurrent = logprior_uniform(params(1,:),prior_bounds);

% calculate the log-likelihood of the model given the data
llcurrent = mod_loglikelihood(t_obs,t_mod,t_err);

%% cast relevant variables before entering the loop
current = params(1,:);
lMAP = -Inf;
mMAP = current;
nacc = 0;

%% decide on transition step size (this is a tedious process that takes trial and error)
% vector as set up such that [Ui Uf Kexp Te Timing n m];
%p_steps = [0.2e-6, 2e-6, 2.5e-2, 5, 1000, 2.5e-2, 6e-4]; % lower accetance rate that allows one to search parameter space, but not optimally tuned
p_steps = [0.2e-6, 2e-6, 1e-3, 5, 1000, 1e-3, 2.5e-4]; % higher acceptance rate (~20-60%), suggesting tuned properly

%% time to run the MCMC loop
h = waitbar(0,'running MCMC...');
for i = 2:burn_in+n_runs
   
    % pick candidate parameters
    candidate = current+p_steps.*randn(1,7);
   
    % set candidate parameters for model run
    Ui_t = candidate(1);
    Uf_t = candidate(2);
    
    Te_t = candidate(4);
    D = E.*Te_t.^3./(12.*(1-v^2));
    lambda = (4.*D./((rhom_rhoc)*g)).^(1/4);
    
    S_U = Uf_t*exp(-S_eud./lambda).*cos(S_eud./lambda);
    
    K_t = (10^(candidate(3)));
    n_t = candidate(6);
    m_t = candidate(7);
    t_int = candidate(5);
    
    % run forward models and put results ito vector
    [Z0] = river_incision_forward_model(S,S_DA,S_U,Ui_t,K_t,m_t,n_t,t_int,dt);
    ter_z = terrace_uplift_forward_model(age,f_dist,lambda,t_int,Ui_t,Uf_t);
    t_mod = [Z0;ter_z];
    
    % calculate model, priors and step for log of acceptance ratio
    lpcandidate = logprior_uniform(candidate,prior_bounds);
    llcandidate = mod_loglikelihood(t_obs,t_mod,t_err);
    
    % transition probabilities
    lr1 = logproposal(candidate,current,p_steps);
    lr2 = logproposal(current,candidate,p_steps);
    
    % add all the probabilities together for the acceptance ratio
    logalpha = lpcandidate + llcandidate + lr1 - lpcurrent - llcurrent - lr2;
    
    % Take the minimum of the log(alpha) and 0.
    if (logalpha > 0)
        logalpha = 0;
    end
    
    % Generate a U(0,1) random number and take its logarithm.
    logt = log(rand());
    
    % Accept or reject the step.
    if (logt < logalpha)
        
        % Accept the step.
        current = candidate;
        %if i > burn_in
            nacc = nacc + 1;
        %end
        
        % Update the MAP solution if this one is better.
        if ((lpcandidate + llcandidate) > lMAP)
            lMAP = lpcandidate + llcandidate;
            mMAP = candidate;
        end
        
        % update current likelihoods if needed
        lpcurrent = lpcandidate;
        llcurrent = llcandidate;
      
    else
        % reject the step
    end
    
    % store the chain paths and other data
    params(i, :) = current;
    like_data(i,1) = llcurrent;
    like_data(i,2) = lpcurrent; 

    % calculate acceptance ratio if past burnin
    if i > burn_in
        accrate = nacc / (n_runs);
        like_data(i,3) = nacc/(i)*100;
    end

    waitbar(i/(burn_in+n_runs),h)
end
close(h)

% save the data
save([cur_dir,'/params_',fileTag,'.mat'],'params');
save([cur_dir,'/mMAP_',fileTag,'.mat'],'mMAP');
save([cur_dir,'/like_dat_',fileTag,'.mat'],'like_data');