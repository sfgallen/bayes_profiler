% this is just a basic script to plot the chain paths, likelihood path,
% acceptance rate history and best-fit results
% 
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 09/28/2020

% clear workspace and command window
clear
clc

% point to current directory for saving files
cur_dir = cd;

% use the same file identifier used for output file naming in the MCMC
% script
fileTag = 'test_model';

% define the timestep in years
dt = 25000;

%% load the results
% parameter values
params = load(['params_',fileTag,'.mat']);
params = params.params;

% likelihood history
like_data = load(['like_dat_',fileTag,'.mat']);
like_data = like_data.like_data;

% maximum a posteriori (MAP) solution
mMAP = load(['mMAP_',fileTag,'.mat']);
mMAP = mMAP.mMAP;

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

%% run the best-fit forward model
% set flexural parameters and define uplift pattern
rhom_rhoc=3300 - 2700;      % difference in density between mantle and crust (kg/m^3
g = 9.81;                   % acceleration due to gravity (m/s^2)
E = 70e9;                   % Young's Modulus (kg/m/yr^2)
v = 0.25;                   % Poisson's ratio

% initial uplift
Ui_t = mMAP(1);

% final uplift at fault
Uf_t = mMAP(2);

% flexed uplift 
Te_t = mMAP(4);
D = E.*Te_t.^3./(12.*(1-v^2));              % flexural rigidity
lambda = (4.*D./((rhom_rhoc)*g)).^(1/4);    % flexural parameter (broken plate)
S_U = Uf_t.*exp(-S_eud./lambda).*cos(S_eud./lambda);  % flexed uplift rate pattern for stream

K_t = (10^(mMAP(3))); % steam erodibility

n_t = mMAP(6);      % slope exponent
m_t = mMAP(7);      % DA exponent
t_int = mMAP(5);    % fault initiation time

%% run the forward models
% stream incision model
[Z0] = river_incision_forward_model(S,S_DA,S_U,Ui_t,K_t,m_t,n_t,t_int,dt);
% terrace uplift model
ter_z = terrace_uplift_forward_model(age,f_dist,lambda,t_int,Ui_t,Uf_t);

%% plot the results
% best-fit forward model
figure(1)
subplot(1,2,1)
plot(S.distance./1e3,Sz,'ko','markerfacecolor',[0.8 0.8 0.8],'markersize',3); hold on
plot(S.distance./1e3,Z0,'ko','markerfacecolor',[0.1 0.2 0.9],'markersize',3); hold on
xlabel('Distance (km)'); ylabel('Elevation (m)');

subplot(1,2,2)
plot(f_dist./1e3,t_up,'ko','markerfacecolor',[0.8 0.8 0.8],'markersize',6); hold on
plot(f_dist./1e3,ter_z,'ko','markerfacecolor',[0.1 0.2 0.9],'markersize',6); hold on
xlabel('Distance from fault (km)'); ylabel('Total uplift (m)');
legend('observed','modeled')

% chain paths
figure(2)
max_it = length(params(:,1));
y_labs = {'Ui', 'Uf', 'log10(K)','Te','Timing','n','m'};
for jj = 1:7
    subplot(7,1,jj)
    plot(params(:,jj),'k-'); hold on
    ylabel(y_labs{jj})
    xlim([0 max_it])
end
xlabel('Iteration')

% likelihood and acceptance rate history
figure(3)
subplot(2,1,1)
plot(like_data(:,1),'k-');
xlim([0 max_it])
ylabel('log-likelihood');
subplot(2,1,2)
plot(like_data(:,3),'k-');
xlim([0 max_it])
ylabel('acceptance rate (%) after burn-in');
xlabel('Iteration');