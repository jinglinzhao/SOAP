[pxx,f] = plomb(coeff(1+1, :), 0:N_CCF);
[pmax,lmax] = max(pxx);
figure; plot(f,pxx)
f0 = f(lmax);
1/f0

% with noise
[pxx_noise,f_noise] = plomb(coeff_noise(1+1, :), 0:N_CCF);
[pmax_noise,lmax_noise] = max(pxx_noise);
figure; plot(f_noise,pxx_noise)
f0_noise = f_noise(lmax_noise);
1/f0_noise



