% Branch from Hermit_coeff_nor.m

% DIFF;

% Integraged from Hermit_DIFF_GI479.m

% Options to run with noise or without noise

% Suppose introducing a shift would not affect the diff; so I introduced
% the shift b, to move all line profiles to the centre. 
% But it turns out that it matters!

%%%%%%%%%%
% Update %
%%%%%%%%%%
% Fix the noise calculation A + normrnd(0, A.^0.5/SN). @25/07/17
% Remove plotting the time series of Her_Spe. @25/07/17
% Introduce the "findpeaks" function -> find the highest few peaks in the periodogram @25/07/17

% Fix the noise calculation (A + normrnd(0, (1-A).^0.5/SN)). @26/07/17
% Change A to 1 - importdata(filename);

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 10000;
N_FILE          = 75;
ORDER           = 5;                                                        % Highest Hermite order 

array_order     = 0:ORDER;
idx_even        = mod(0:ORDER, 2) == 0;
order_even      = array_order(idx_even);
order_odd       = array_order(~idx_even);

SHIFT           = 11;
Her_Spe         = zeros(N_FILE, ORDER);
Her_Spe_rvc     = zeros(N_FILE, ORDER);
coeff           = zeros((ORDER+1), SHIFT);
coeff_rvc       = zeros((ORDER+1), SHIFT);

V_BEGIN         = -40;                                                      % m/s                                                     
V_END           = -V_BEGIN;
V_GRID          = (V_BEGIN : (V_END - V_BEGIN)/(SHIFT-1) : V_END) / 1000;   % km/s

grid_size       = 0.1;
v               = (-20 : grid_size : 20)';                                  % km/s
idx             = (v<10) & (v>-10);
v               = v(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
h               = waitbar(0,'Please wait...');
for n = 1:N_FILE
    
    v_planet    = 5 * sin(n/25/0.618*2*pi + 1) * 0.001; % km/s
    
    i           = n - 1;
    filename    = ['../CCF_dat/ccf', num2str(i), '.dat'];
    A           = 1 - importdata(filename);
    A           = A(idx);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    b           = f.b;  % shift    

    for order = 0:ORDER
        for shift = 1:SHIFT
            temp                        = hermite_nor(order, v - V_GRID(shift) + v_planet) * grid_size;
            temp_rvc                    = hermite_nor(order, v - b - V_GRID(shift)) * grid_size;
            coeff(order+1, shift)       = sum(temp .* (A + normrnd(0, (1-A).^0.5/SN)));     
            coeff_rvc(order+1, shift)   = sum(temp_rvc .* (A + normrnd(0, (1-A).^0.5/SN))); 
            % coeff(order+1, shift)   = sum(A .* temp);     
            % coeff_rvc(order+1, shift) = sum(A .* temp_rvc);             
        end
    end

    % fitting %
    for order = order_even
        p = polyfit(V_GRID * 1000, coeff(order+1, :), 2);
        Her_Spe(n, order+1) = p(1);

        p_rvc = polyfit(V_GRID * 1000, coeff_rvc(order+1, :), 2);
        Her_Spe_rvc(n, order+1) = p_rvc(1);          
    end

    for order = order_odd
        p = polyfit(V_GRID * 1000, coeff(order+1, :), 1);
        Her_Spe(n, order+1) = p(1);
        
        p_rvc = polyfit(V_GRID * 1000, coeff_rvc(order+1, :), 2);
        Her_Spe_rvc(n, order+1) = p_rvc(1);          
    end
    
    waitbar( n / N_FILE )
end
close(h)  

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
cd ../
for order = 0:ORDER
    
    [pxx,f] = plomb(Her_Spe(:, order+1), 1:N_FILE, 0.25, 100, 'normalized');
    [pmax,lmax] = max(pxx);
    f0 = f(lmax);
    disp(['T_planet: ', num2str(1/f0)]);

    [pxx_rvc,f_rvc] = plomb(Her_Spe_rvc(:, order+1), 1:N_FILE, 0.25, 100, 'normalized');
    [pmax_rvc,lmax_rvc] = max(pxx_rvc);
    f0_rvc = f_rvc(lmax_rvc);
    disp(['T_activity: ', num2str(1/f0_rvc)]);

    [pks,locs] = findpeaks(pxx, f);                 % find all the peaks in (pxx, f)
    [pks_maxs, idx_maxs] = sort(pks, 'descend');    % sort "pks" in descending order; mark the indexies 
    [pks_rvc,locs_rvc] = findpeaks(pxx_rvc, f_rvc);            
    [pks_maxs_rvc, idx_maxs_rvc] = sort(pks_rvc, 'descend'); 

    h = figure;
    
        hold on
        plot(f, pxx, 'r')
        plot(f_rvc, pxx_rvc, 'b')
        legend('Rest frame', 'Observed frame', 'Location', 'Best')
        
        for i = 1:5
            x = locs(idx_maxs(i));  % locations of the largest peaks -> harmonics
            y = pks_maxs(i);
            text(x, y, ['\leftarrow', num2str(1/x, '%3.2f')]);
            
            x_rvc = locs_rvc(idx_maxs_rvc(i));  % locations of the largest peaks -> harmonics
            y_rvc = pks_maxs_rvc(i);
            text(x_rvc, y_rvc, ['\leftarrow', num2str(1/x_rvc, '%3.2f')]);            
        end
        
        xlabel('Frequency')
        ylabel('Normalized Power')
        title_name = ['Order', num2str(order)];
        title(title_name);
    hold off

    out_eps = [title_name, '.eps'];
    print(out_eps, '-depsc')
    close(h);    
end