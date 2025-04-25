function[gamma1,gamma2,rho2,nu,LL_Matern]=EstimateLocalMaternCovariance(M_latentData,V_Lat,V_Lon,V_a0,nb_max_days,nb_max_dryGauges)

%parameter estimation - maximize pairwise log-likelihood
param0=zeros(3,1);
param0(1)=50; %gamma1 (rq: rho1=sqrt(gamma1^2+gamma1^2))
param0(2)=50; %gamma2
param0(3)=80; %rho2
param0(4)=0.5;
rho1_ini=sqrt(param0(1)^2+param0(2)^2);

A=[[0 0 -1 0];[0 0 1 0];[0 0 0 -1];[0 0 0 1]];
b=[-15; 150; -0.3; 2.0];

fun = @(param)(-1*log_likelihood_mattern(M_latentData,V_Lon,V_Lat,V_a0,param(1),param(2),param(3),param(4),nb_max_days,nb_max_dryGauges));
myconfun=@(param)(mycon(param(1),param(2),param(3),15,150,5));
param = fmincon(fun,param0,A,b,[],[],[],[],myconfun);

gamma1=param(1);
gamma2=param(2);
rho2=param(3);
nu=param(4);

[LL_Matern]=log_likelihood_mattern(M_latentData,V_Lon,V_Lat,V_a0,gamma1,gamma2,rho2,nu,nb_max_days,nb_max_dryGauges);

%Pairwise log-likelihood
    function[LL_Matern]=log_likelihood_mattern(M_data,V_Lon,V_Lat,V_a0,my_gamma1,my_gamma2,my_rho2,my_nu,nb_max_days,nb_max_dryGauges)
        LL_Matern=0;

        my_d=sqrt(my_gamma1^2+my_gamma2^2);
        M_Lambda=[[my_d*my_d 0];[0 my_rho2*my_rho2]];
        M_Gamma=[[my_gamma1/my_d -my_gamma2/my_d];[my_gamma2/my_d my_gamma1/my_d]];
        Sig_aniso=M_Gamma*M_Lambda*M_Gamma';
        inv_Sig_aniso=inv(Sig_aniso);

        for my_day=1:min(size(M_data,2),nb_max_days)

            V_obs_latent=M_data(:,my_day);
            
            inds_keep=~isnan(V_obs_latent);
            V_obs_latent=V_obs_latent(inds_keep);
            V_Lon_day=V_Lon(inds_keep);
            V_Lat_day=V_Lat(inds_keep);
            V_a0_day=V_a0(inds_keep);
            
            inds_Rain=find(V_obs_latent>-5);
            inds_NoRain=find(V_obs_latent==-5);
            inds_NoRain=inds_NoRain(1:min(nb_max_dryGauges,length(inds_NoRain)));
            
            %pair with 2 dry
            for i=1:length(inds_NoRain)
                for j=i+1:length(inds_NoRain)
                    V_Coord_diff=[V_Lon_day(inds_NoRain(i))-V_Lon_day(inds_NoRain(j));V_Lat_day(inds_NoRain(i))-V_Lat_day(inds_NoRain(j))]*111;
                    h=sqrt(V_Coord_diff'*inv_Sig_aniso*V_Coord_diff);
                    if h<30
                        rho=1/(gamma(my_nu)*2^(my_nu-1))*((h)^(my_nu))*besselk(my_nu,h);
                        Sigma=[[1 rho];[rho 1]];
                        LL_Matern=LL_Matern+max(-30,log(mvncdf([V_a0(inds_NoRain(i));V_a0(inds_NoRain(j))],[0;0],Sigma)));
                    end
                end
            end
            %pair with 1 dry and 1 wet
            for i=1:length(inds_NoRain)
                for j=1:length(inds_Rain)
                    Zj=V_obs_latent(inds_Rain(j));
                    V_Coord_diff=[V_Lon_day(inds_NoRain(i))-V_Lon_day(inds_Rain(j));V_Lat_day(inds_NoRain(i))-V_Lat_day(inds_Rain(j))]*111;
                    h=sqrt(V_Coord_diff'*inv_Sig_aniso*V_Coord_diff);
                    if Zj<3.5 && h<30
                        rho=1/(gamma(my_nu)*2^(my_nu-1))*((h)^(my_nu))*besselk(my_nu,h);
                        Sigma=[[1 rho];[rho 1]];
                        LL_Matern=LL_Matern+max(-30,log(normcdf(V_a0_day(inds_NoRain(i)))-mvncdf([V_a0_day(inds_NoRain(i));V_a0_day(inds_NoRain(i))],[0;0],Sigma)))+log(normpdf(Zj));
                    end
                end
            end
            %pair with 2 wet
            for i=1:length(inds_Rain)
                for j=i+1:length(inds_Rain)
                    Zi=V_obs_latent(inds_Rain(i));
                    Zj=V_obs_latent(inds_Rain(j));
                    V_Coord_diff=[V_Lon_day(inds_Rain(i))-V_Lon_day(inds_Rain(j));V_Lat_day(inds_Rain(i))-V_Lat_day(inds_Rain(j))]*111;
                    h=sqrt(V_Coord_diff'*inv_Sig_aniso*V_Coord_diff);
                    if Zi<3.5 && Zj<3.5 && h<30
                        rho=1/(gamma(my_nu)*2^(my_nu-1))*((h)^(my_nu))*besselk(my_nu,h);
                        Sigma=[[1 rho];[rho 1]];
                        LL_Matern=LL_Matern+max(-30,log(mvnpdf([Zi;Zj],[0;0],Sigma)));
                    end
                end
            end
        end
        
end %end function LL_Matern


    function [c,ceq] = mycon(gam1,gam2,rho2,rho1_min,rho1_max,rho_ratio)
        c(1)=gam1^2+gam2^2-rho1_max^2;
        c(2)=-gam1^2-gam2^2+rho1_min^2;
        c(3)=sqrt(gam1^2+gam2^2)-rho2; %rho1 <= rho2
        c(4)=rho2/rho_ratio-sqrt(gam1^2+gam2^2);%rho2 <= rho_ratio*rho1
        ceq=[];
    end

end