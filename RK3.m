function a = RK3(u,v,w,Nx,Ny,Re,H,dt)
K1 = zeros(Nx+2);
for i = 2:Nx+1
    for j = 2:Ny+1
        K1(i,j) = -u(i,j)*(w(i+1,j) - w(i-1,j))*(0.5*Nx) - v(i,j)*(w(i,j+1) - w(i,j-1))*(0.5*Ny) + ...
            (H\Re)*(w(i+1,j)+w(i-1,j) - 4*w(i,j) + w(i,j+1)+w(i,j-1)) ;
    end
end
w = w + (dt/3)*K1 ;
for i = 2:Nx+1
    for j = 2:Ny+1
        K1(i,j) = -u(i,j)*(w(i+1,j) - w(i-1,j))*(0.5*Nx) - v(i,j)*(w(i,j+1) - w(i,j-1))*(0.5*Ny) + ...
            (H\Re)*(w(i+1,j)+w(i-1,j) - 4*w(i,j) + w(i,j+1)+w(i,j-1)) - (5/9)*K1(i,j) ;
    end
end
w = w + (15*dt/16)*K1 ;
for i = 2:Nx+1
    for j = 2:Ny+1
        K1(i,j) = -u(i,j)*(w(i+1,j) - w(i-1,j))*(0.5*Nx) - v(i,j)*(w(i,j+1) - w(i,j-1))*(0.5*Ny) + ...
            (H\Re)*(w(i+1,j)+w(i-1,j) - 4*w(i,j) + w(i,j+1)+w(i,j-1)) - (153/128)*K1(i,j) ;
    end
end
a = w + 8*dt*K1/15 ;