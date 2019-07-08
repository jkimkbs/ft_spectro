%% 10 Interferogram waves
% James Kim, Antoine wojdyla

% 7/8/2019

%% define a time scale

x_m = 0:0.01e-6:100e-6;

%% Defining Variables: wavelength

lambda0_m = 10e-6;
dlambda_m = 0.1e-6;

%% Interferences
    
% number of wavelength
Nl = 10;
E = zeros(1,length(x_m));

for i_l=1:Nl
    lambda_m = lambda0_m + (i_l-1)*dlambda_m;
    E = E+sin(2*pi/lambda_m*x_m);
end

plot(x_m*1e6,E)
xlabel('optical path [\mum]')
ylabel('electric field')
title('ten wavelength superimposed')

%% Now you can do interferences....

d_m = 0e-6;

Nl = 10;
E0 = zeros(1,length(x_m));
Ed = zeros(1,length(x_m));
for i_l=1:Nl
    lambda_m = lambda0_m + (i_l-1)*dlambda_m;
    E0 = E0+sin(2*pi/lambda_m*x_m);
    Ed = Ed+sin(2*pi/lambda_m*(x_m-d_m));
    I = abs(E0+Ed).^2;
end

plot(x_m*1e6,Ed)
xlabel('optical path [\mum]')
ylabel('electric field')
title('ten wavelength superimposed and shifted')

%%
% central wavelength
lambda0_m = 0.5e-6;
% wavelength increment
dlambda_m = 0.05e-6;
% number of wavelengths
Nl = 10;

% delay stage steps
dds_m = 100e-9;
% maximum delay
dsm_m = 1000e-6;

% delays
ds_m = (-dsm_m):dds_m:dsm_m;
Nd = length(ds_m);

E0 = zeros(1,length(x_m));
I = zeros(length(dds_m));

for i_d=1:Nd
    Ed = zeros(1,length(x_m));
    Id = 0;
    for i_l=1:Nl
        lambda_m = lambda0_m + (i_l-1)*dlambda_m;
        E0 = sin(2*pi/lambda_m*x_m);
        Ed = sin(2*pi/lambda_m*(x_m-ds_m(i_d)));
        Id = Id + sum(abs(E0+Ed).^2);
    end
    I(i_d) = Id;
end

%%
plot(ds_m*1e6,I)
xlabel('stage delay [\mum]')
ylabel('intensity [a.u.]')
title(sprintf('interferogram; lambda0 = %1.1fnm, dlambda = %1.2fnm, Nlambda = %i',...
    lambda0_m*1e9,dlambda_m*1e9,Nl))

%% Fourier Transform spectroscopy
FTXR = abs(fftshift(fft(I))).^2;
f_cpm = linspace(-1/(2*dds_m),1/(2*dds_m),length(ds_m));
plot(f_cpm*1e-6,FTXR)
xlabel('frequency [\mum^-^1]')
ylabel('PSD of interferogram')
title(sprintf('FTIR; l0 = %1.1fnm, dl=%1.2fnm, Nl = %i; Ds = %1.1fum, ds = %1.1fnm',...
    lambda0_m*1e9,dlambda_m*1e9,Nl,dsm_m*1e6, dds_m*1e9))


