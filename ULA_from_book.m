% Uniform Linear Array Simulation trying to follow the book

clc; clear all; close all;

N = 21;            % number of receivers

f = 1e9;           % Hertz
c = 299792458;     % speed of light
lambda = c/f;
d = lambda/2;
w = 2*pi*f;
deg2rad = pi/180;

beam_angle = 20;

% Placing the array along the z-axis as in the book (for displaying
% purpose)

Tx_pos = zeros(N,3);
Tx_pos(:,3) = ([0:N-1]-(N-1)/2)*d; % Eq. 2.53

% Array manifold vector vk
% theta is polar angle defined as coming from the z axis down

%theta = 20;       % degrees (need to be swept)
theta = -180:1:180;
freq_steps = length(theta);

kz = -2*pi/lambda * cos(theta*deg2rad);
vk = zeros(freq_steps,N);
for i = 1:freq_steps  % sweep through all the degrees
    kz = -2*pi/lambda * cos(theta(i)*deg2rad);
    vk(i,:) = exp(j*(((N-1)/2)-[0:N-1])*kz*d);  % Eq. 2.55
end

%%%% Weighting to change sidelobe level and mainbeam-width %%%%%

n_tilda = ([0:N-1]-(N-1)/2);

% cosine weighting Eq. 3.12
%weights = sin(pi/(2*N))*cos(pi*n_tilda/N); 

% raised cosine Eq. 3.16 
p = 0;
c_p = p/N+((1-p)/2)*sin(pi/(2*N));
weights = c_p*(p+(1-p)*cos(pi*n_tilda/N));

% Hamming weighting
weights = 0.54+0.46*cos(2*pi*n_tilda/N);
weights = weights/sum(weights);




phase_delay = exp((j*2*pi*d.*[0:N-1]./lambda)*sin(beam_angle*deg2rad));
w_h  = repmat(phase_delay.*ones(1,N)/N,freq_steps,1);
w_h2  = repmat(phase_delay.*weights,freq_steps,1);

% Sigma frequency-wavenumber function in sigma space

sigma = w_h.*vk; %freq_steps x N 
sigma = sum(sigma,2);

sigma2 = w_h2.*vk; %freq_steps x N 
sigma2 = sum(sigma2,2);

%%%% plotting %%%%
figure
subplot(2,1,1)
plot(theta,real(sigma))
title(['Uniform Linear Array Gain w/ N = ', num2str(N)])
xlabel('degrees')
ylabel('linear')
hold on
plot(theta,real(sigma2),'--')
subplot(2,1,2)
polar(theta*deg2rad,abs(sigma)')
title(['steering angle ', num2str(beam_angle) ' degrees'])
hold on
polar(theta*deg2rad,abs(sigma2)','--')

%%% theoretical %%%%
radians = 0:pi/100:2*pi;
sig = 2*pi*cos(radians)*d/lambda;
sigma2 = (1/N)*(sin(N.*sig/2))./sin(sig/2); %theoretical value

%figure
%plot(radians, sigma2)




