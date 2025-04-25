function[V_Y]=Transform_R_to_Y(V_R,V_a0,V_k,V_teta)
V_Y=NaN(size(V_R));
for my_gauge=1:length(V_R)
    if V_R(my_gauge)==0
        V_Y(my_gauge)=-5.0;
    else

    end

    if V_R(my_gauge)>0
        my_prop_zeros=normcdf(V_a0(my_gauge));
        my_uniform=gamcdf(V_R(my_gauge),V_k(my_gauge),V_teta(my_gauge));
        my_gaussian=norminv(max(0.0001,min(my_prop_zeros+(1-my_prop_zeros)*my_uniform,0.9999)));
        V_Y(my_gauge)=my_gaussian;
    end
end

end %end function