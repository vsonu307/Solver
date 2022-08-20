function A = Gaussseidel(Nx,Ny,H,psi,L2,w)
psi_n = zeros(Nx+2) ;
while L2 > 10^(-5)
    for i = 2:Nx+1
        for j = 2:Ny+1
            psi_n(i,j) = (0.25)*(w(i,j)/H + (psi_n(i,j-1) + psi_n(i-1,j) + psi(i+1,j) + psi(i,j+1))) ;
        end  
    end
%     psi_n(Nx+2,:) = psi_n(2,:) ;     psi_n(:,Ny+2) = psi_n(:,2) ;
%     psi_n(1,:) = psi_n(Nx+1,:) ;     psi_n(:,1) = psi_n(:,Ny+1) ;
    L2 = sqrt(sum(sum((psi_n(2:Nx+1,2:Ny+1) - psi(2:Nx+1,2:Ny+1)).^2))) ;
    psi = psi_n ;
    A = psi;
end
