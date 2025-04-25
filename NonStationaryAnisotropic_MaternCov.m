function[M_Cov]=NonStationaryAnisotropic_MaternCov(V_Lon1,V_Lon2,V_Lat1,V_Lat2,V_g1_1,V_g1_2,V_g2_1,V_g2_2,V_rho2_1,V_rho2_2,V_nu_1,V_nu_2)

M_Cov=NaN(length(V_Lon1),length(V_Lon2));
for i=1:size(M_Cov,1)
    for j=1:size(M_Cov,2)

        if sqrt((V_Lon1(i)*111-V_Lon2(j)*111)^2+(V_Lat1(i)*111-V_Lat2(j)*111)^2)<0.01 
            M_Cov(i,j)=1;
        else

            di=sqrt(V_g1_1(i)^2+V_g2_1(i)^2);
    
            M_Lambda_i=[[di*di 0];[0 V_rho2_1(i)*V_rho2_1(i)]];
            M_Gamma_i=[[V_g1_1(i)/di -V_g2_1(i)/di];[V_g2_1(i)/di V_g1_1(i)/di]];
            Sig_i=M_Gamma_i*M_Lambda_i*M_Gamma_i';

            dj=sqrt(V_g1_2(j)^2+V_g2_2(j)^2);
            M_Lambda_j=[[dj*dj 0];[0 V_rho2_2(j)*V_rho2_2(j)]];
            M_Gamma_j=[[V_g1_2(j)/dj -V_g2_2(j)/dj];[V_g2_2(j)/dj V_g1_2(j)/dj]];
            Sig_j=M_Gamma_j*M_Lambda_j*M_Gamma_j';

            Delta_X=(V_Lon1(i)-V_Lon2(j))*111;
            Delta_Y=(V_Lat1(i)-V_Lat2(j))*111;
            Qij=[Delta_X Delta_Y]*inv((Sig_i+Sig_j)/2)*[Delta_X;Delta_Y];

            V_nu_1(i)=0.5;
            V_nu_2(j)=0.5;
            
            nu_ij=(V_nu_1(i)+V_nu_2(j))/2;

            phy_xy=(det(Sig_i)^(1/4))*(det(Sig_j)^(1/4))*(det((Sig_i+Sig_j)/2)^(-1/2));

            M_Cov(i,j)=2^(1-nu_ij) / sqrt(gamma(V_nu_1(i))*gamma(V_nu_2(j))) * phy_xy *(sqrt(Qij)^nu_ij)*besselk(nu_ij,sqrt(Qij));
            
        end
    end
end

end