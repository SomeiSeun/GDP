function [v_sqrt] = VisViva(r,a,mu)
% Outputs square root of velocity using vis-viva relation
v_sqrt = sqrt(mu.*((2./r) - (1./a)));

end