function lp = logprior_uniform(candidate,prior_bounds)
% log prior assuming uniform priors
log_test = zeros(size(candidate));
for i = 1:length(candidate)
    if candidate(i) >= prior_bounds(i,1) &&  candidate(i) <= prior_bounds(i,2)
        log_test(i) = 1;
    end
end

if sum(log_test) == length(candidate)
    lp=0;
else
    lp = -inf;
end
end