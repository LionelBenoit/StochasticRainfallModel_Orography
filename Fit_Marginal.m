function[V_Par]=Fit_Marginal(V_data)

prop_dry=sum(V_data==0)/length(V_data);
V_data_nonzero=V_data(V_data>0);

%initialization
a0=min(max(norminv(prop_dry),-4),4);
paramIntensity0(1)=1;
paramIntensity0(2)=1;

A=[[-1 0];[1 0];[0 -1]];
b=[-0.01;10;-2.01];

%fit
fun = @(paramIntensity)(-1*log_likelihood_marginalModel(V_data_nonzero,paramIntensity));
paramIntensity = fmincon(fun,paramIntensity0,A,b);
V_Par(1)=a0;
V_Par(2)=paramIntensity(1);
V_Par(3)=paramIntensity(2);

    function[LL_marginal]=log_likelihood_marginalModel(V_Robs,V_Par)
        if length(V_Par)==2
            k=V_Par(1);
            teta=V_Par(2);
            LL_marginal=0;
            for i=1:length(V_Robs)
                Ri=V_Robs(i);
                liky=1/(gamma(k)*(teta^k))*Ri^(k-1)*exp(-Ri/teta);
                LL_marginal=LL_marginal+log(max(liky,0.000001));
            end
        else
            disp('wrong number of parameters')
            LL_marginal=NaN;
        end
    end

end %end function