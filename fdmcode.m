function temperature_plots_same_as_paper()  

clear; clc; close all; 

% --- Parameters (as used in the paper's numerical example)
rho = 1; Cp = 1; kappa0 = 0.1; 
chi_values = [-1.0 -0.5 0 0.5 1.0]; 

a = 1; b = 3; N = 41; h = (b-a)/(N-1); 

tEnd = 15; M = 31; tau = tEnd/(M-1); 

alpha = 2; % u(a,t)
beta  = 1; % u(b,t)

% --- Grid + initial condition (same form as paper) 
x = zeros(N,1); 
u0 = zeros(N,1); 
for i = 1:N
    x(i) = a + (i-1)*h; 
    u0(i) = 2 - (x(i)-1)/2 + (x(i)-1)*(x(i)-3); 
end 

t = zeros(M,1); 
for n = 1:M
    t(n) = (n-1)*tau; 
end 

% --- Newton settings (paper uses eps < 1e-4 stopping metric)
tol = 1e-4; 
maxIt = 50; % safety guard

A = rho*Cp/tau; 

for c = 1:length(chi_values) 
    chi = chi_values(c); 

    % Solution storage U(i,n) 
    U = zeros(N,M); 
    U(:,1) = u0; 

    % Arrays for Newton system 
    u     = zeros(N,1); 
    u_1   = zeros(N,1); 
    uNext = zeros(N,1); 
    G     = zeros(N,1); 
    L     = zeros(N,N); 

    % Boundary Jacobian rows 
    L(1,1) = 1; 
    L(N,N) = 1; 

    % ---- Time stepping ---- 
    for n = 2:M 
        % Initial guess u^(0) = u^{n-1} 
        u   = U(:,n-1); 
        u_1 = U(:,n-1); 

        eps = 1; 
        it = 0; 

        % ---- Newton iterations ---- 
        while (eps > tol) && (it < maxIt) 
            it = it + 1; 

            % Boundary residuals (Dirichlet) 
            G(1) = u(1) - alpha; 
            G(N) = u(N) - beta; 

            % Interior equations 
            for i = 2:N-1 
                % k(u) = k0 * exp(chi*u) 
                k  = kappa0 * exp(chi*u(i)); 
                Dk = chi*k; 
                D2k = chi*Dk;  % chi^2*k 

                % v = (u_{i+1} - u_{i-1})/(2h) 
                v = (u(i+1) - u(i-1)) / (2*h); 

                % phi and f (paper appendix) 
                phi = A*(u(i) - u_1(i)) - Dk*v*v; 
                f   = phi / k; 

                % p and q (paper appendix) 
                q = (-f*Dk + A - D2k*v*v) / k; 
                p = -2*Dk*v / k; 

                % Residual G_i (FDM form) 
                G(i) = u(i+1) - 2*u(i) + u(i-1) - h*h*f; 

                % Jacobian entries 
                L(i,i-1) = 1 + 0.5*h*p; 
                L(i,i)   = -2 - h*h*q; 
                L(i,i+1) = 1 - 0.5*h*p; 
            end 

            % Newton update: uNext = u - L^{-1} G 
            uNext = u - L\G; 

            % Paper stopping metric: 
            eps = sqrt(h * (uNext - u)' * (uNext - u)); 

            u = uNext; 
        end 

        if it == maxIt 
            warning('Newton hit maxIt for chi=%g at time step n=%d (t=%g). eps=%g', ...
                chi, n, t(n), eps);
        end 

        U(:,n) = u; 
    end 

    % ---- Plots ---- 
    figure; 
    mesh(x, t, U'); 
    xlabel('x'); ylabel('t'); zlabel('Temperature'); 
    title(['Temperature Evolution, \chi = ', num2str(chi)]); 

    figure; 
    plot(x, U(:,end), 'LineWidth', 2); 
    grid on; 
    xlabel('x'); ylabel('Temperature'); 
    title(['Final Temperature Distribution, \chi = ', num2str(chi)]); 
end 

end
