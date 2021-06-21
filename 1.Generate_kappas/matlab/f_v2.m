function y= f_v2(z1,u,z0,z,pars,p1,C)
% z1 is the integral up boundary
% u is random number
% z is hidden varible of the super node
% pars=[alpha,beta,c,mu] in the stable distribution.
y=integral(@(x) C*stable_pdfC(x, p1,1).*stable_pdfC(z-x, p1,1)/stable_pdfC(z, pars,1), z0, z1)-u;
    
    