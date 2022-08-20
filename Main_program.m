clear ;
Nx = 64;    dx = 1/Nx ;
Ny = 64;    dy = 1/Ny ;
X = linspace(0,1,Nx+2) ;
Y = linspace(0,1,Nx+2) ;
H = 1/(dx^2);    K = 1/(dy^2);  % H = K
L2 = 1; dt = 0.0001 ;
t = (0:dt:0.01);
Re = 1000 ;
psi = rand(Nx+2);   % Initializing psi
psi_n = zeros(Nx+2);
%% Calculation of u & v
u = zeros(Nx+2) ;
v = zeros(Ny+2) ;
for i = 2:Nx+1
    for j = 2:Ny+1
        u(i,j) = (psi(i,j+1)-psi(i,j-1))/(2*dy) ;
        v(i,j) = (psi(i-1,j)-psi(i+1,j))/(2*dx) ;
    end
end
%% Initial calculation of w
w = zeros(Nx+2) ;
for i = 2:Nx+1
    for j = 2:Ny+1
        w(i,j) = (0.5*Nx)*(v(i+1,j) - v(i-1,j)) - (0.5*Ny)*(u(i,j+1) - u(i,j-1)) ;
    end
end


%% calculations

for k = 1:length(t)
    a = RK3(u,v,w,Nx,Ny,Re,H,dt);
    w = a ;
    A = Gaussseidel(Nx,Ny,H,psi,L2,w);
    psi = A ;
    for i = 2:Nx+1
        for j = 2:Ny+1
            u(i,j) = (0.5*Ny)*(psi(i,j+1)-psi(i,j-1)) ;
            v(i,j) = (0.5*Nx)*(psi(i-1,j)-psi(i+1,j)) ;
        end
    end
    
    
    %% Boundary Condtion Update
    u(Nx+2,:) = u(2,:) ;    u(1,:) = u(Nx+1,:) ;    u(:,Ny+2) = u(:,2) ;    u(:,1) = u(:,Ny+1) ;
    v(Nx+2,:) = v(2,:) ;    v(1,:) = v(Nx+1,:) ;    v(:,Ny+2) = v(:,2) ;    v(:,1) = v(:,Ny+1) ;
    w(Nx+2,:) = w(2,:) ;    w(1,:) = w(Nx+1,:) ;    w(:,Ny+2) = w(:,2) ;    w(:,1) = w(:,Ny+1) ;
    psi(Nx+2,:) = psi(2,:) ;    psi(1,:) = psi(Nx+1,:) ;    psi(:,Ny+2) = psi(:,2) ;    psi(:,1) = psi(:,Ny+1) ;
    
    
    %% plot and Animation
    clf
    hold on
    contourf(X,Y,w) 
    shading interp
    colorbar
    title(' 2D Turbulence ')
    drawnow
    MovieVector(k) = getframe ;
end
mywriter = VideoWriter('Vorticity','MPEG-4');
mywriter.FrameRate = 20 ;
open(mywriter)
writeVideo(mywriter,MovieVector)
close(mywriter)

