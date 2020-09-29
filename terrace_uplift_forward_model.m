function [ter_z] = terrace_uplift_forward_model(age,f_dist,lambda,fault_init_time,Ui,Uf)
% this function calculates the total uplift experenced by a terrace since
% the time of terrace formation based on a simple two stage uplift history.
% The first stage of uplift is assumed to be due to spatially uniform
% uplift. The second stage is modeled a broken plate flexure and uplift is
% parameterized as terrace Euclidean distance from the broken segement
% (e.g. a normal fault).
%
% Inputs:
% age - terrace ages in kyr
% f_dist - Euclidean distance from broke segment
% lambda - flexural parameter
% fault_init_time - fault initiation time in years
% Ui - initial block uplift rate in m/yr
% Uf - uplift rate at broken segment in m/yr
%
% Outputs:
% ter_z - total rock uplift of terrace
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date Modified: 09/28/2020

terrace_unit_age = unique(age);
ter_z = nan(size(f_dist));

for i = 1:length(terrace_unit_age)
    age_inds = find(age == terrace_unit_age(i));
    terrace_age = terrace_unit_age(i).*1e3;
    
    if fault_init_time(1) < terrace_age
        dUpt = terrace_age - fault_init_time;
        tU1 = Ui.*dUpt;
        ter_z(age_inds) = (Uf*exp(-f_dist(age_inds)'./lambda).*cos(f_dist(age_inds)'./lambda))*fault_init_time+tU1;
    else
        ter_z(age_inds) = (Uf*exp(-f_dist(age_inds)'./lambda).*cos(f_dist(age_inds)'./lambda))*terrace_age;
    end
end
end