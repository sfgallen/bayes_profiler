function [lr] = logproposal(p_1,p_2,step)
% log gaussian transition probability
lr = (-1/2)*sum((p_1-p_2).^2./step.^2);
end