function l = mod_loglikelihood(obs,mod,sig)
% gaussian log-likelihood function
err = obs-mod;
l = (-1/2).*sum((err./sig).^2);
end