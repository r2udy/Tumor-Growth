clear all
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

%% Data saving param
saving = 0;

%% Parameters
%Filling coefficients by the user
pde_number = 0;
% khi = input('(u) -> khi = ')
% mu = input('(u^2) -> mu = ')
% eta = input('(u^3) -> eta = ')
% kappa = input('(u^4) -> kappa = ')
nu = input('(u_x) -> nu = ')
epsilon = input('(u_xx) -> epsilon = ');
% omega = input('(v_xx) -> omega = ');
iota = input('(u_xxxx) -> iota = ')
% lbda = input('(u*u_x) -> lbda = ')
% psi = input('Div(v*h) -> psi = ')

if ~exist("khi", "var")
    khi = 0;
end
if ~exist("mu", "var")
    mu = 0;
end
if ~exist("eta", "var")
    eta = 0;
end
if ~exist("kappa", "var")
    kappa = 0;
end
if ~exist("nu", "var")
    nu = 0;
end
if ~exist("omega", "var")
    omega = 0;
end
if ~exist("epsilon", "var")
    epsilon = 0;
end
if ~exist("iota", "var")
    iota = 0;
end
if ~exist("lbda", "var")
    lbda = 0;
end
if ~exist("psi", "var")
    psi = 0;
end

Lib_params = {{'\khi', ' h', khi}, {'\mu', ' h^2', mu}, {'\eta', ' h^3', eta}, {'\kappa', ' h^4', kappa}, {'\nu', ' ^{dh}/_{dx}', nu}, {'\omega', ' ^{d²v}/_{dx²}', omega}, {'\epsilon', ' ^{d²h}/_{dx²}', epsilon}, {'\iota', ' ^{d^{4}h}/_{dx^{4}}', iota}, {'\lambda', ' h ^{dh}/_{dx}', lbda}, {'\psi', ' \nabla(v h)', psi}};


%% Parameters : Filling coefficients from the generated file
% cd D:\CODING\data\'IF Files'\PDElearn_heights_filter\density_equation_filter\
% fields = dir();
% list_filename_intermd = [fields.name];
% list_filename = regexp(list_filename_intermd, 'pde\_stability\_([a-z]+[0-9]?\_)?[0-9]+\.txt', 'match');
% 
% %Choice of the PDE
% pde_number = input('pde identifier (integer) ->  ');
% 
% %Reach the coefficients of the PDE
% pat = '0' + string(pde_number);
% for p=1:size(list_filename,2)
%     if size(strfind(list_filename{p}, pat),1)>0
%         file_pat = string(list_filename{p});
%         break
%     end
% end
% coeffs = fileread(file_pat)
% matching_col1 = regexp(coeffs, '[+-]?[0-9]+\.[0-9]*', 'match');
% n_coeffs = size(matching_col1,2);
% 
% %Some coeffficent aren't in the first file of coefficients
% %Build the library of terms with the coefficients
% Lib_params = {{'\khi', ' h', 0}, {'\mu', ' h^2', 0}, {'\eta', ' h^3', 0}, {'\kappa', ' h^4', 0}, {'\nu', ' ^{dh}/_{dx}', 0}, {'\epsilon', ' ^{d²h}/_{dx²}', 0}, {'\iota', ' ^{d^{4}h}/_{dx^{4}}', 0}, {'\lambda', ' h ^{dh}/_{dx}', 0}, {'\psi', ' \nabla(v h)', 0}};
% for l=1:n_coeffs
%     Lib_params{l}{3} = str2double(matching_col1{l});
% end
% 
% khi = Lib_params{1}{3};
% mu = Lib_params{2}{3};
% eta = Lib_params{3}{3};
% kappa = Lib_params{4}{3};
% nu = Lib_params{5}{3};
% omega = Lib_params{6}{3};
% epsilon = Lib_params{7}{3};
% iota = Lib_params{8}{3};
% lbda = Lib_params{9}{3};
% psi = Lib_params{10}{3};
% 
% %% Data
% load("D:\CODING\data\IF Files\Surface_growth_defectors06_07_2023 13_30_145.txt_(x,t).mat")

N = 256
x = 32*pi*(1:N)'/N; % this is used in the plotting
% s = RandStream('mlfg6331_64');
% u = randsample(s,heights(:,1),N);
% u = exp(-(x-x(end)/2).^2/(2*epsilon^2)) / sqrt(2*pi*epsilon^2);
u = 10*sin(pi*x/x(end));
v = 0;
h_fft = fft(u);
v_fft = fft(v);

%% Precompute various scalars for ETDRK4
% dt = diff(t_heights(1:2))
dt = 0.01
k = [0:N/2-1 0 -N/2+1:-1]'/16; % wave numbers ( 2*pi*k/Long ) 
L = 1j*nu*k - epsilon*k.^2 + iota*k.^4; % Fourier multipliers
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
uu = u; vv = v; tt = 0;
tmax = 10; nmax = round(tmax*1/dt); nplt = floor((tmax/100)/dt);
g = 1j*lbda*k;
for n = 1:nmax
    t = n*dt;
%     display(n);
    %Computing of the solution at time t
    Nh = g.*fft(real(ifft(h_fft)).^2) + khi*fft(real(ifft(h_fft))) + mu*fft(real(ifft(h_fft)).^2) + eta*fft(real(ifft(h_fft)).^3) + kappa*fft(real(ifft(h_fft)).^4);
    
    a = E_2.*h_fft + Q.*Nh;
    Na = g.*fft(real(ifft(a)).^2) + khi*fft(real(ifft(a))) + mu*fft(real(ifft(a)).^2) + eta*fft(real(ifft(a)).^3) + kappa*fft(real(ifft(a)).^4);
        
    b = E_2.*h_fft + Q.*Na;
    Nb = g.*fft(real(ifft(b)).^2) + khi*fft(real(ifft(b))) + mu*fft(real(ifft(b)).^2) + eta*fft(real(ifft(b)).^3) + kappa*fft(real(ifft(b)).^4);
    
    c = E_2.*h_fft + Q.*(2*Nb-Nh);
    Nc = g.*fft(real(ifft(c)).^2) + khi*fft(real(ifft(c))) + mu*fft(real(ifft(c)).^2) + eta*fft(real(ifft(c)).^3) + kappa*fft(real(ifft(c)).^4);
    
    h_fft = E.*h_fft + alpha.*Nh + 2*beta.*(Na+Nb) + gamma.*Nc; 
    
    %Updating of the solution array
    u = real(ifft(h_fft));
    uu = [uu,u]; tt = [tt,t];
end

close all;

%% Plot
%Plot for every time t
figure(3)
surf(tt,x,uu), shading interp, lighting phong, axis tight
colormap(summer); set(gca,'zlim')
light('color',[1 1 0],'position',[-1,2,2])
material([0.30 0.60 0.60 40.00 1.00]);

the_title = ['PDE ' num2str(pde_number) ' : ^{d h}/_{d t} ='];
the_parameters = '';
for key = Lib_params
    if key{1}{3} ~= 0.0
        the_parameters = [the_parameters, [key{1}{1}, '=', num2str(key{1}{3}), ', ']];
        the_title = [the_title, [' + ', key{1}{1}, key{1}{2}]];
    end
end
title(the_title, the_parameters)
xlabel('t'); ylabel('x')
colorbar

% Plot for t = t_{final}
figure(2)
plot(x, [uu(:,1), uu(:,end)])
title(the_title, [ the_parameters ',  for t=t_{final}'])
xlabel('t'); ylabel('x')
legend('t = t_{0}','t = t_{finale}')
legend show

%% Test Fourier derivative 
% Nt = size(heights,1);
% h_fft = fft(heights);
% kt = 2*pi*fftshift(-Nt/2:Nt/2 -1)/Nt;
% KT = kt(ones(size(heights,2), 1), :)';
% dh_fft = 1j * h_fft .* KT;
% dh_t = real(ifft(dh_fft));
% 
% Nx = size(heights,2);
% h_fft = fft(heights(50,:));
% kx = 2*pi*fftshift(-Nx/2:Nx/2 -1)/Nx;
% KX = kx(ones(size(heights,1), 1), :);
% dh_fft = 1j * h_fft .* kx;
% dh_x = real(ifft(dh_fft));



% %% Plot the movie
% figure(3);
% y = uu;
% x = x;
% j_suppr = find(max(y)> 2*mean(y));
% %y(:, j_suppr) = [];
% height_max = max(max(y)) + 0.1*max(max(y));
% height_min = min(min(y));
% frame = getframe(gcf);
% fps = 40;
% video_name = 'video_PDEL_x_test';
% v = VideoWriter(video_name, 'MPEG-4'); 
% v.FrameRate = fps;
% open(v);               
% for time=1:size(y,2)
%     plot(x, y(:,time))
%     ylim([height_min, height_max])
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end   
% close(v);

%% Data Saving
if saving
    mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_h_t.mat'), 'heights_t');
    mfw(diff(uu'));
    release(mfw)
    mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_h.mat'), 'heights');
    mfw(uu(:, 1:end-1)');
    release(mfw)
    mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_x.mat'), 'x_heights');
    mfw(x);
    release(mfw)
    mfw = dsp.MatFileWriter(sprintf("Verif_algo" + '_t.mat'), 't_heights');
    mfw(tt(1:end-1)');
    release(mfw)
%     heights = uu;
%     t_heights = tt;
%     x_heights = x;
%     save(sprintf('Surface_growth_' + "verifalgo" + '_(x,t).mat'), 'heights', 't_heights', 'x_heights');
end