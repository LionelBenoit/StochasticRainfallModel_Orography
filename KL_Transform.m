function[M_EigVect,M_PCs]=KL_Transform(M_Data_latent,M_EigVect)
X=M_Data_latent';
if isempty(M_EigVect)
    [W,~] = eig(X'*X); %eigenvector decomposition of variance-covariance matrix
    M_EigVect=fliplr(W);
end
M_PCs=X*M_EigVect;
M_EigVect=M_EigVect(:,1:3);

end