% Branch from Hermit_coeff_nor.m
% Correct the RV center before Hermite decomposition

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 10000;
N_CCF           = 74;                               % number of CCF files
N_hermite       = 5;                                % Highest Hermite order 
coeff           = zeros((N_hermite+1), (N_CCF+1));
coeff_noise     = zeros((N_hermite+1), (N_CCF+1));
coeff_rvc       = zeros((N_hermite+1), (N_CCF+1));
coeff_noise_rvc = zeros((N_hermite+1), (N_CCF+1));
grid_size       = 0.1;
v               = (-20 : grid_size : 20)';          % km/s
RV              = importdata('RV.dat') / 1000;      % activity induced RV [km/s]
RV_gauss        = zeros(size(0:N_CCF))';

idx             = (v > -10) & (v < 10);
v               = v(idx);

% template %
A_tpl           = 1 - importdata('CCF_tpl.dat');
A_tpl           = A_tpl(idx);
f_tpl           = fit( v, A_tpl, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
b_tpl           = f_tpl.b;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
for n_CCF = 0:N_CCF

    v_planet    = 5 * sin(n_CCF/25/0.618*2*pi + 1) * 0.001; % km/s
    
    filename    = ['CCF_dat/CCF', num2str(n_CCF), '.dat'];
    A           = 1 - importdata(filename);
    A           = A(idx);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    b           = f.b;  % shift
    RV_gauss(n_CCF+1) = b;

    disp([n_CCF, b*1000, (b-b_tpl)*1000])

    for n_hermite = 0:N_hermite
        temp                                = A .* hermite_nor(n_hermite, v - b_tpl + v_planet) * grid_size;
        coeff(n_hermite+1, n_CCF+1)         = sum(temp);                                % Without noise
        coeff_noise(n_hermite+1, n_CCF+1)   = sum(temp .* (1+normrnd(0, A.^0.5/SN)));   % With noise

        temp_rvc                                = A .* hermite_nor(n_hermite, v - b) * grid_size;
        coeff_rvc(n_hermite+1, n_CCF+1)         = sum(temp_rvc);                                % Without noise
        coeff_noise_rvc(n_hermite+1, n_CCF+1)   = sum(temp_rvc .* (1+normrnd(0, A.^0.5/SN)));   % With noise
    end

end     
RV_gauss = RV_gauss - b_tpl;


for n_hermite = 0:N_hermite
    
    [pxx_noise,f_noise] = plomb(coeff_noise(n_hermite+1, :), 0:N_CCF);
    [pmax_noise,lmax_noise] = max(pxx_noise);
    f0_noise = f_noise(lmax_noise);
    disp(['T_planet: ', num2str(1/f0_noise)]);

    [pxx_noise_rvc,f_noise_rvc] = plomb(coeff_noise_rvc(n_hermite+1, :), 0:N_CCF);
    [pmax_noise_rvc,lmax_noise_rvc] = max(pxx_noise_rvc);
    f0_noise_rvc = f_noise_rvc(lmax_noise_rvc);
    disp(['T_activity: ', num2str(1/f0_noise_rvc)]);

    h = figure;
        hold on
        plot(f_noise, pxx_noise / max(pxx_noise), 'r')
        plot(f_noise_rvc, pxx_noise_rvc / max(pxx_noise_rvc), 'b')
        plot(f0_noise, pmax_noise / max(pxx_noise), 'ro')
        
        text(f0_noise, pmax_noise / max(pxx_noise), ['\leftarrow', num2str(1/f0_noise)]);
        plot(f0_noise_rvc, pmax_noise_rvc / max(pxx_noise_rvc), 'bo')
        text(f0_noise_rvc, pmax_noise_rvc / max(pxx_noise_rvc), ['\leftarrow', num2str(1/f0_noise_rvc)])
        legend('Rest frame', 'Observed frame', 'Location', 'Best')
        xlabel('Frequency')
        ylabel('Normalized Power')
        title_name = ['Order', num2str(n_hermite)];
        title(title_name);
    hold off

    out_eps = [title_name, '.eps'];
    print(out_eps, '-depsc')
    close(h);
end

if 0
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficient Plotting %
    %%%%%%%%%%%%%%%%%%%%%%%%
    x = 0:N_CCF;

    for n_hermite = 0:N_hermite

        h = figure;
        hold on
            plot(x, coeff(n_hermite+1, :), 'r')
            plot(x, coeff_noise(n_hermite+1, :), 'rs--', 'MarkerSize',6, 'MarkerFaceColor', 'r')
            plot(x, coeff_rvc(n_hermite+1, :), 'b')
            plot(x, coeff_noise_rvc(n_hermite+1, :), 'bo--', 'MarkerSize',6, 'MarkerFaceColor', 'b')  
            legend('Rest frame', 'Rest frame w n', 'Shifted frame', 'Shifted frame w n', 'Location', 'Best')
            xlabel('Observation Number','Interpreter','latex')
            Y_label = ['a', num2str(n_hermite)];
            ylabel(Y_label,'Interpreter','latex')
            TITLE1 = ['SN', num2str(SN)];
            title_name = ['Order', num2str(n_hermite), '--', TITLE1];
            title(title_name)
        hold off

        out_eps = [title_name, '_nor.eps'];
        print(out_eps, '-depsc')
        close(h);
    end


    %%%%%%%%%%%%%%%%%%%%%
    % Correlation Plots %
    %%%%%%%%%%%%%%%%%%%%%

    for n_hermite = 0:N_hermite

        h = figure;
        hold on 
            plot(RV_gauss * 1000, coeff(n_hermite+1, :), 'r')
            plot(RV_gauss * 1000, coeff_rvc(n_hermite+1, :), 'b')
            plot(RV_gauss * 1000, coeff_noise(n_hermite+1, :), 'rs', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
            plot(RV_gauss * 1000, coeff_noise_rvc(n_hermite+1, :), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b')
            xlabel('RV_a (m/s)')
            ylabel(['a', num2str(n_hermite)])
            legend('Ideal', 'Observed Frame, S/N = 10000', 'Location', 'northeast')
            legend('Rest frame', 'Observed Frame', 'Location', 'southeast')
            title_name = ['Order', num2str(n_hermite), ' -- ', TITLE1];
            title(['a', num2str(n_hermite)])
        hold off

        out_eps = [title_name, '-nor.png'];
        print(out_eps, '-dpng')
        close(h);
    end   
end



