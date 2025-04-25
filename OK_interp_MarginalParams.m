function[V_targetpoints]=OK_interp_MarginalParams(X_datapoints, Y_datapoints, V_datapoints, X_targetpoints, Y_targetpoints)


%Experimental variogram
pts=[X_datapoints Y_datapoints V_datapoints]; %check the distance unit is km
[meanh,gammah,~]=Variogram_expe_2d(pts,15,40,1);
V_lags_vario=meanh;
V_values_vario=gammah;
sig2_emp=var(V_datapoints);

%Fit variogram (Matern model)
param0=zeros(4,1);
param0(1)=1;%nu
param0(2)=20;%r
%param0(4)=V_values_vario(1)+0.01;%Nugget
param0(4)=sig2_emp*0.05;
param0(3)=sig2_emp-param0(4);%sill

A=[[-1 0 0 0];[0 -1 0 0];[1 0 0 0];[0 1 0 0];[0 0 1 1];[0 0 0 1];[0 0 0 -1]];
%b=[-0.3;-5;2; 150;sig2_emp;sig2_emp/2;-V_values_vario(1)*0.9];
b=[-0.3;-5;2; 150;sig2_emp;sig2_emp*0.1;-sig2_emp*0.01];

fun = @(param)(RMSE_vario_Matern(V_lags_vario,V_values_vario,param(1),param(2),param(3),param(4)));
param = fmincon(fun,param0,A,b);

nu=param(1);
r=param(2);
sig2=param(3);
Nugget=param(4);

%build covariance matrix
MVX1=repmat(X_datapoints,1,length(X_datapoints));
MVX2=repmat(X_datapoints',length(X_datapoints),1);
MVY1=repmat(Y_datapoints,1,length(Y_datapoints));
MVY2=repmat(Y_datapoints',length(Y_datapoints),1);

M_ds=sqrt((MVX2-MVX1).^2+(MVY2-MVY1).^2);

clear MVX1
clear MVX2
clear MVY1
clear MVY2

Sigma=2^(1-nu)/gamma(nu).*((M_ds/r).^(nu)).*besselk(nu,M_ds/r);
Sigma(M_ds==0)=1;
Cab=Nugget*double(M_ds==0)+sig2*Sigma;
Cab=[[Cab ones(size(Cab,1),1)];[ones(1,size(Cab,2)) 0]];

MVX1=repmat(X_datapoints,1,length(X_targetpoints));
MVX2=repmat(X_targetpoints',length(X_datapoints),1);
MVY1=repmat(Y_datapoints,1,length(Y_targetpoints));
MVY2=repmat(Y_targetpoints',length(Y_datapoints),1);

M_ds=sqrt((MVX2-MVX1).^2+(MVY2-MVY1).^2);

Sigma=2^(1-nu)/gamma(nu).*((M_ds/r).^(nu)).*besselk(nu,M_ds/r);
Sigma(M_ds==0)=1;
Cav=Nugget*double(M_ds==0)+sig2*Sigma;
Cav=[Cav;ones(1,size(Cav,2))];

Za=V_datapoints;
Za=[Za;0];

Lambda_a=inv(Cab)*Cav;
Z_est=Za'*Lambda_a;

V_targetpoints=Z_est';

function[RMSE]=RMSE_vario_Matern(V_lags_vario,V_values_vario,my_nu,my_r,sig2,Nugget)
    V_vario_Mod=Nugget+sig2*(1-2^(1-my_nu)/gamma(my_nu).*((V_lags_vario/my_r).^(my_nu)).*besselk(my_nu,V_lags_vario/my_r));
    V_delta=V_values_vario-V_vario_Mod;
    RMSE=sqrt(mean(V_delta.^2));
end

end