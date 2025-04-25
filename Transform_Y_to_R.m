function[V_R]=Transform_Y_to_R(V_Y,V_a0,V_k,V_teta)

V_R=NaN(size(V_Y));
for my_gauge=1:length(V_Y)
    if V_Y(my_gauge) > V_a0(my_gauge)
        p=(normcdf(V_Y(my_gauge))-normcdf(V_a0(my_gauge)))/(1-normcdf(V_a0(my_gauge)));
        V_R(my_gauge)=gaminv(p,V_k(my_gauge),V_teta(my_gauge));
    else
        V_R(my_gauge)=0;
    end
end

end %end function
