function [pars_fit]=Find_kappa(filepath, layers,ps)
%nohup matlab -nodisplay -nosplash -r 'Find_kappa('./data/TW_65',10,0.5)' > outfile1.txt &
    loadlibrary('libstable','stable.h');
    %filepath='./data/TW_65'; 
   
    data=load([filepath '_coordinates.txt']);
    parameters=load([filepath '_parameters.txt']);    
    %parameters=[N, beta, mu]
    beta=parameters(2);    
    kappa=data(:,2);
    
    
    z=kappa.^beta;%kappa has been sorted in data, and so z also sorted  
    
    [pars_fit]=stable_fit_mle2dC(z, [], 1);%alpha stabel fitting parameters
    disp(pars_fit);
    
    
    %output fitting results
    Pc_data=Empirical_ccdf(z);
    ycdf = stable_cdfC(Pc_data(:,1), pars_fit, 1);
    Pc_fit =(1-ycdf(1:end))';
    
    fid10= fopen([filepath '_fitting_performance.txt'], 'w');  
    for j= 1:length(z)      
        fprintf(fid10, '%10f %10f %10f\n', [Pc_data(j,1) Pc_data(j,2) Pc_fit(j)]);
    end
    fclose(fid10);      
    
    fid0= fopen([filepath '_fitting_stable_parameter.txt'], 'w');          
    fprintf(fid0, '%10f %10f %10f %10f\n', [pars_fit(1) pars_fit(2) pars_fit(3) pars_fit(4)]);
    fclose(fid0); 
    
    
    pars=pars_fit;
    %disp(length(z));    
    %output z and kappa distribution in layer 0
    l=0;
    Pc3=Empirical_ccdf(z);
    Pc4=Empirical_ccdf(kappa);
    fid2= fopen([filepath '_z_kappa_l_' num2str(l) '.txt'], 'w');  
    for j= 1:length(Pc3)      
        fprintf(fid2, '%10f %10f %10f %10f\n', [Pc3(j,1) Pc3(j,2) Pc4(j,1) Pc4(j,2)]);
    end
    fclose(fid2);    
    
    %cal z uplayr and the related distribution
    for l=1:layers
        
        %z0=min(z)*2^(-beta/(xgamma-1.0));
        if pars(1)>1.0
            z0=0.4*min(z);
        else
        	z0=min(z)*2^(-1/(pars(1)));
        end
        [zuplayer]= cal_zplus_zminus_z_p_split_v1(z,pars,z0,ps);       
        kappa=zuplayer(:,3).^(1.0/beta);
        disp(l);
        
        %output z and kappa dist in uplayer
        Pc3=Empirical_ccdf(zuplayer(:,3));
        Pc4=Empirical_ccdf(kappa);      
        fid2= fopen([filepath '_z_kappa_l_' num2str(l) '.txt'], 'w');  
        for j= 1:length(Pc3)      
            fprintf(fid2, '%10f %10f %10f %10f\n', [Pc3(j,1) Pc3(j,2) Pc4(j,1) Pc4(j,2)]);
        end
        fclose(fid2); 

        %very important output!! [subnode supernode kappa]
        fid= fopen([filepath '_kappa_l_' num2str(l) '.txt'], 'w');  
        for j= 1:length(kappa)      
            fprintf(fid, '%10d %10d %16f\n', [zuplayer(j,1) zuplayer(j,2) kappa(j)] );
        end
        fclose(fid); 
        
        
        z=zuplayer(:,3);    
        alpha1=pars(1);
        beta1=pars(2);
        c1=pars(3)/((1+ps)^(1/pars(1)));
        mu1=pars(4)/(1+ps);
        pars=[alpha1,beta1,c1,mu1];
    end
quit;

