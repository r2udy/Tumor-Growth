% A plot of the Kuramoto-Sivashinsky equation:
%
% u_t = -u*u_x - u_xx - u_xxxx, periodic boundary conditions on [0,32*pi]
%
% 
% solving = 0 ----> u_t = - u_xx - u_xxxx, periodic boundary conditions on [0,32*pi]
% 
% solving = 1 ----> u_t = - u_xx - u_xxxx - 1/2 * uu_x, periodic boundary conditions on [0,32*pi]
% 
% solving = 2 ----> u_t = - u_xx - u_xxxx + gaussian_noise , periodic boundary conditions on [0,32*pi]
%
% computation is based on v = fft(u), so linear term is diagonal
%
% Using this program:
% u is the initial condition
% dt is the time step
% N is the number of points calculated along x

%% Parameters 
nu = input('(u_x) -> nu = ')
lbda = input('(u*u_x) -> lbda = ')

%Choice of the equation
solving = 1;

% Initial condition and grid setup
N = 128
x = 32*pi*(1:N)'/N; % this is used in the plotting
u = cos(x/16).*(1+sin(x/16)); % initial condition
h_fft = fft(u);

% Precompute various scalars for ETDRK4
dt = 0.25 % time step
k = [0:N/2-1 0 -N/2+1:-1]'/16; % wave numbers ( 2*pi*k/L ) 
L = nu*k.^2 - lbda*k.^4; % Fourier multipliers
E = exp(dt*L); E_2 = exp(dt*L/2);
M = 16;
r = exp(1i*pi*((1:M)-0.5)/M);
Lr = dt*L(:, ones(M,1)) + r(ones(N,1),:);

Q     = dt*real(mean( (exp(Lr/2) - 1)./Lr, 2 ));
alpha = dt*real(mean( (-4 - Lr + exp(Lr).*(4-3*Lr+Lr.^2))./Lr.^3 , 2));
beta  = dt*real(mean( (2 + Lr + exp(Lr).*(-2 + Lr))./Lr.^3 , 2));
gamma = dt*real(mean( (-4 - 3*Lr - (Lr).^2 + exp(Lr).*(4 - Lr))./Lr.^3 , 2));

gaussian_noise = wgn(length(h_fft), 1, -10, 'complex');

%% Time-Stepping loop:

uu = u; tt = 0;
tmax = t_heights(end); nmax = round(tmax/dt); nplt = floor((tmax/100)/dt);
% g = -1i*lbda*k;
g = 1i*lbda*k;
g_prime = -1i*zeta*k.^2;
g_tild = -1i*kappa*k.^2;
for n = 1:nmax
    t = n*dt;
    if solving==0
        the_title = 'The Kuramoto Sivashinsky equation : ^{d h}/_{d t} = -\nu ^{d²h}/_{dx²} - \lambda ^{d^{4}h}/_{dx^{4}} ';
        h_fft = E.*h_fft;

    elseif solving==1
        the_title = 'The Kuramoto Sivashinsky equation : ^{d h}/_{d t} = -\nu ^{d²h}/_{dx²} - \lambda ^{d^{4}h}/_{dx^{4}} - ^{1}/_{2} (^{d h}/_{d x})^2';
        Nh = g.*fft(real(ifft(h_fft)).^2);
    
        a = E_2.*h_fft + Q.*Nh;
        Na = g.*fft(real(ifft(a)).^2);
        
        b = E_2.*h_fft + Q.*Na;
        Nb = g.*fft(real(ifft(b)).^2);
    
        c = E_2.*h_fft + Q.*(2*Nb-Nh);
        Nc = g.*fft(real(ifft(c)).^2);
    
        h_fft = E.*h_fft + alpha.*Nh + 2*beta.*(Na+Nb) + gamma.*Nc;

    else
        the_title = 'The Kuramoto Sivashinsky equation : ^{d h}/_{d t} = -\nu ^{d²h}/_{dx²} - \lambda ^{d^{4}h}/_{dx^{4}} + \eta';
        h_fft = E.*h_fft + gaussian_noise;

    end
    if mod(n,nplt)==0
        u = real(ifft(h_fft));
        uu = [uu,u]; tt = [tt,t];
    end
%     u = real(ifft(h_fft));
%     uu = [uu,u]; tt = [tt,t];
end

close all;

%% Plot
surf(tt,x,uu), shading interp, lighting phong, axis tight
colormap(summer); set(gca,'zlim')
light('color',[1 1 0],'position',[-1,2,2])
material([0.30 0.60 0.60 40.00 1.00]);
the_parameters = [ ' \nu=' num2str(nu) ', \lambda=' num2str(lbda)];
title(the_title, the_parameters)
xlabel('t'); ylabel('x')
colorbar

% figure(2)
% nf_visu = 200;
% title(the_title,'Seen in section for t \in [0,'+tt(nf_visu)+']')
% for i=1:20:nf_visu
%     txt = ['t=', num2str(i)];
%     plot(x, uu(:, i), 'DisplayName', txt);
%     hold on;
% end
% legend show
