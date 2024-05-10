function [x1,y1,h7,h81,h82,x83,y83,h83] = estimate(target,D,N,berth)
% assuming transmitter at the origin and lambda_0 = 1
% target: target cartesian coordinate vector [x,y]
% D: separation distance of the receivers
% wavelengths: vector of wavelengths to transmit

Tx = [0,0];
l_0 = 1;

if nargin == 3
    berth = 5*l_0;
end

spacing = l_0/8;
xlim = [target(1)-berth target(1)+berth];
ylim = [target(2)-berth target(2)+berth];

% create array of wavelengths
wavelengths = zeros([1 N]);
for n = 1:N
    wavelengths(n) = N/(n + N/2);
end

% gather data
R1 = [Tx(1) - D, 0];
R2 = [Tx(1) + D, 0];
R1m = [Tx(1) - D/2, 0];
R2m = [Tx(1) + D/2, 0];

r0 = sqrt((target(1)-Tx(1))^2+(target(2)-Tx(2))^2);
r1 = sqrt((target(1)-R1(1))^2+(target(2)-R1(2))^2);
r2 = sqrt((target(1)-R2(1))^2+(target(2)-R2(2))^2);

R1data = zeros([1 N]);
R2data = zeros([1 N]);
for i = 1:N
    lambda = wavelengths(i);
    % R1data(i) = (1j*lambda*r)^-0.5*exp(-1j*pi*r/lambda) % using magnitude
    R1data(i) = exp(1j*pi*(r0+r1)/lambda);
    R2data(i) = exp(1j*pi*(r0+r2)/lambda);
end

% bistatic model
[x1,y1,h1] = reconstruct(wavelengths,spacing,spacing,xlim,ylim,R1,Tx,R1data);
[x1,y1,h2] = reconstruct(wavelengths,spacing,spacing,xlim,ylim,R2,Tx,R2data);
h7 = h1 + h2;
%disp('Bistatic model estimation found')

% monostatic model
[x1,y1,h1] = reconstruct(wavelengths,spacing,spacing,xlim,ylim,R1m,R1m,R1data);
[x1,y1,h2] = reconstruct(wavelengths,spacing,spacing,xlim,ylim,R2m,R2m,R2data);
h81 = h1 + h2;
%disp('Monostatic model estimation found')

% range and bearing angle model
Rrdata = R1data.*R2data;
Rbdata = R1data.*conj(R2data);
[x1,y1,h1] = reconstruct(wavelengths,spacing,spacing,xlim,ylim,Tx,Tx,Rrdata,'range');
[x1,y1,h2] = reconstruct(wavelengths,spacing,spacing,xlim,ylim,R1,R2,Rbdata,'bearing angle');
h82 = h1 + h2;
%disp('Range and bearing angle model estimation found')

% fourier model
Nfft = N;
scale_factor = round(N/(2*D));
Nfft_scale = Nfft*scale_factor;
rr = abs(fft(Rrdata,Nfft));
bb = abs(fft(Rbdata,Nfft_scale));

domain = [1:D*scale_factor+1 Nfft_scale-D*scale_factor+1:Nfft_scale];

x3p = zeros(0);
y3p = zeros(0);
h83 = zeros(0);
for i = 1:length(domain)
    for j = 1:Nfft
        if domain(i) > Nfft_scale/2
            x3p = [x3p -(Nfft_scale-(domain(i)-1))*90/(scale_factor*D)];
        else
            x3p = [x3p (domain(i)-1)*90/(scale_factor*D)];
        end
        y3p = [y3p (j-1)/2];
        h83 = [h83 bb(domain(i))+rr(j)];
    end
end

% convert polar to cartesian
x3c = zeros(size(x3p));
y3c = zeros(size(y3p));
for i = 1:length(x3c)
    theta = x3p(i);
    r = y3p(i);
    x3c(i) = r*cos((90-theta)*pi/180);
    y3c(i) = r*sin((90-theta)*pi/180);
end
%disp('Fourier model estimation found')

x83 = [x3p; x3c];
y83 = [y3p; y3c];

end