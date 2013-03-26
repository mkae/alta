function rho = ABS(theta_in, phi_in, theta_out, phi_out, a, b, c)

	fx = sin(theta_out)*cos(phi_out) - sin(theta_in);
	fy = sin(theta_out)*sin(phi_out);
	f  = sqrt(fx*fx + fy*fy);

	psd = a / (1+b*b * f*f)^(c);
	rho = 4*pi*pi * (cos(theta_in)+cos(theta_out))^2 * psd;

endfunction
