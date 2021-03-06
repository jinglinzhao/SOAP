% Branch from Hermit_coeff_nor.m

% Simulate with planets

%%%%%%%%%%
% Update %
%%%%%%%%%%
% Introduce the "findpeaks" function -> find the highest few peaks in the periodogram @25/07/17
% Fix the noise calculation (A + normrnd(0, (1-A).^0.5/SN)). @26/07/17
% Change indexing in order to implement parfor @12/08/17

% Based on Hermit_coeff_NOR_RV_correct_0720.m @07/12/17

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 10000;
N_FILE          = 75;                               % number of CCF files
grid_size       = 0.1;
v               = (-20 : grid_size : 20)';          % km/s
v0              = v;
RV              = importdata('../RV.dat') / 1000;      % activity induced RV [km/s]
RV_gauss        = zeros(N_FILE,1);

idx             = (v > -10) & (v < 10);
v               = v(idx);

% template %
A_tpl           = 1 - importdata('../CCF_tpl.dat');
A_tpl           = A_tpl(idx);
f_tpl           = fit( v, A_tpl, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
b_tpl           = f_tpl.b;
% idx_x           = (v < b_tpl + 10) & (v > b_tpl);
idx_x           = (v < b_tpl + 15) & (v > b_tpl - 15);
v_x             = v(idx_x);
y_fit           = zeros(N_FILE, length(v_x));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:N_FILE
    
    i           = n - 1;
    v_planet    = 20 * sin(i/25*0.38*2*pi + 1) * 0.001; % km/s
    
    filename    = ['../CCF_dat/CCF', num2str(i), '.dat'];
    A           = 1 - importdata(filename);
    A0          = A;
    A           = A(idx);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    b           = f.b;  % shift
    RV_gauss(n) = b;
    disp([i, b*1000, (b-b_tpl)*1000])

    y_fit(n, :) = spline(v0 + v_planet, A0 + normrnd(0, (1-A0).^0.5/SN), v_x);
end     

RV_gauss = RV_gauss - b_tpl;

cd ../

if 0
    
    figure;
    hold on 
    for i = 1:length(v_x)
        plot(1:N_FILE, y_fit(:, i), '.')
    end
    hold off
    
    
    figure;
    hold on 
    for i = 1:length(v_x)
        plot(1:N_FILE, y_fit(:, i) - y_fit(1,i), '.')
    end
    hold off

    figure;
    hold on 
    for i = 1:length(v_x)
        SUM = y_fit(:, i) - y_fit(1,i);
        plot(i, SUM, '.')
    end
    hold off    
    
    h = figure;
    hold on 
    for i = 1:length(v_x)
        plot(1:N_FILE, y_fit(:, i) / mean(y_fit(:, i)), '.')
    end
    hold off

    figure;
    hold on 
        plot(1:length(v_x), sum(y_fit), '.')
    hold off

    figure;
    hold on 
    for i = 1:length(v_x)
        plot(i, sum(y_fit(:, i) - y_fit(1,i)), '.')
    end
    hold off
end


if 0
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficient Plotting %
    %%%%%%%%%%%%%%%%%%%%%%%%


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



