theta = -0.5*pi:0.01:0.5*pi;

a = 1; b =10; c=1.5;
cost = cos(theta);
sint = sin(theta);
f = sint;

psd = a ./ (1 + b*b * f.*f).^(0.5*(c+1)) ;

rho = cost.^2 .* psd ;

Frho = fftshift(fft(rho)) / max(size(rho));
Prho = sqrt(Frho .* conj(Frho));


figure(1); plot(rho);
figure(2); plot(Prho);
figure(3); plot(imag(Frho) ./ real(Frho));
