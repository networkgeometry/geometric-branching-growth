function [zuplayer]= cal_zplus_zminus_z_p_split_v1(X,pars,z0,ps)
%split supernode as a probability ps 
% X is the emparical dataset-z; i.e., the hidden varibles of the super nodes
% pars=[alpha,beta,c,mu] in the stable distribution.
tic;%record the start time 
N=length(X);
zuplayer=[];
sum=0;
subnode=0;
sum_C=0;
for i=1:N
%      disp(i);
    z=X(i);    
    if rand<ps        
        u=rand;    
        if z0>=z-z0
            z_subnode1=0.5*z;
            z_subnode2=0.5*z;
            zuplayer=[zuplayer;[subnode i-1 z_subnode1]];
            subnode=subnode+1;
            zuplayer=[zuplayer;[subnode i-1 z_subnode2]];
            subnode=subnode+1;
            sum=sum+1;
        else
            %zplus(i)=solve_bisection(u,z0,z,pars);
            alpha1=pars(1);
            beta1=pars(2);
            c1=pars(3)/((1+ps)^(1/pars(1)));
            mu1=pars(4)/(1+ps);
            p1= [alpha1,beta1,c1,mu1];
            
            if stable_pdfC(z, pars,1)<1e-20 && z<4*abs(pars(4)+1)
                sum_C=sum_C+1;
                z_subnode1=0.5*z;
                z_subnode2=0.5*z;
            else            
                C=1.0/integral(@(x) stable_pdfC(x, p1,1).*stable_pdfC(z-x, p1,1)/stable_pdfC(z, pars,1),z0, z-z0);
                if C>1e15 || C<1e-15
                    fprintf('C=%f, ',C);
                    fprintf('z=%f, ',z);
                    fprintf('z-z0=%f \n',z-z0);
                    sum_C=sum_C+1;
                    z_subnode1=0.5*z;
                    z_subnode2=0.5*z;
                else
                    [z_subnode1, fval, exitflag, output]=fzero(@(z1) f_v2(z1,u,z0,z,pars,p1,C),[z0,z-z0]);%find the root in range [0,z]
                    if exitflag==1
                        z_subnode2=z-z_subnode1; 
                    else
                        disp('can not converge!')
                    end 
                end
            end
            zuplayer=[zuplayer;[subnode i-1 z_subnode1]];
            subnode=subnode+1;
            zuplayer=[zuplayer;[subnode i-1 z_subnode2]];
            subnode=subnode+1;            
        end        
    else
        zuplayer=[zuplayer;[subnode i-1 z]];
        subnode=subnode+1;
    end
end 
disp(sum);
fprintf('Num zero pdf=%f \n',sum_C);
elapsed = toc;
fprintf('Time: %g\n', elapsed); %record the end time 
