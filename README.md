function T = analytical_solution(nx, nt, nF, L, t_final, plotFlag)

% spatial grid
x = linspace(0, L, nx);

% time grid
t = linspace(0, t_final, nt);

% diffusivity
D = 1.6e-6;

% initialise solution
T = zeros(nx, nt);

for j = 1:nt
    for i = 1:nx
        
        sum_term = 0;
        
        for n = 1:nF
            
            lambda_n = n*pi/L;
            
            Bn = -84/(n*pi);
            
            sum_term = sum_term + Bn * sin(n*pi*x(i)/L) * exp(-(lambda_n^2)*D*t(j));
            
        end
        
        % full solution
        T(i,j) = sum_term - 42*x(i)/L + 2;
        
    end
    
    % optional plotting (animation)
    if plotFlag == 1
        plot(x, T(:,j), 'LineWidth', 2)
        xlabel('x (m)')
        ylabel('Temperature (°C)')
        title(['Time = ', num2str(t(j)/86400), ' days'])
        axis([0 256 -50 5])
        drawnow
    end
    
end

end
