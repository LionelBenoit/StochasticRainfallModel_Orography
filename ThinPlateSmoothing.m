function[V_Z_targetpoints]=ThinPlateSmoothing(V_X_datapoints,V_Y_datapoints,V_Z_datapoints,V_X_targetpoints,V_Y_targetpoints,alpha_thin_plate)

M1=eye(length(V_X_datapoints))*alpha_thin_plate;

for i=1:length(V_X_datapoints)
    for j=1:length(V_X_datapoints)
        if i~=j
            rij=sqrt((V_X_datapoints(j)-V_X_datapoints(i))^2+(V_Y_datapoints(j)-V_Y_datapoints(i))^2);
            if rij>0.01
                aij=(rij^2)*log(rij);
                M1(i,j)=aij;
                M1(j,i)=aij;
            else
                M1(i,j)=1;
                M1(j,i)=1;
            end
        end
    end
end
M2=[M1 ones(size(V_X_datapoints)) V_X_datapoints V_Y_datapoints];
M3=[[ones(size(V_X_datapoints))'; V_X_datapoints'; V_Y_datapoints'], zeros(3,3)];
M4=[M2;M3];

MZ=[V_Z_datapoints;0;0;0];

M_lambda=inv(M4)*MZ;
%--
V_Z_targetpoints=NaN(size(V_X_targetpoints));
for i=1:length(V_Z_targetpoints)
    tmp=0;
    for j=1:length(V_Z_datapoints)
        rij=sqrt((V_X_datapoints(j)-V_X_targetpoints(i))^2+(V_Y_datapoints(j)-V_Y_targetpoints(i))^2);
        if rij>0.01
            aij=(rij^2)*log(rij);
        else
            aij=1;
        end
        tmp=tmp+M_lambda(j)*aij;
    end
    %V_Z_targetpoints(i)=tmp+sum(M_lambda(end-2:end));
    V_Z_targetpoints(i)=tmp+M_lambda(end-2)+M_lambda(end-1)*V_X_targetpoints(i)+M_lambda(end)*V_Y_targetpoints(i);
end

end