% Clean workspace
clear all; close all; clc

% Preamble from CP1_sample.m.  The autograder will take care of any files
% needed for this step.  Please only submit solution file.

load('Kraken.mat') %load data matrix
L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x; %create 3D axis arrays with 64 points
k = (2*pi/(2*L))*[0:(n/2 - 1), -n/2:-1]; %create frequency array and rescale them to be 2pi periodic 
ks = fftshift(k); %shift values to order them correctly

%Create 3D grids for both spatial domain and frequency domain 
[X,Y,Z] = meshgrid(x,y,z);
[Kx,Ky,Kz] = meshgrid(ks,ks,ks);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum all realizations in frequency space
% after your for loop ends, save the sum as variable A1

average = zeros(n, n, n);
for f = 1:49
    KrakenTensor(:, :, :) = reshape(Kraken(:, f), n, n, n);
    fKrakenTensor = fftn(KrakenTensor);
    average = fKrakenTensor + average;
end

A1 = average;

% Average the sum over the 49 realizations (i.e., A1/49) and save as A2

average = average / 49;

A2 = average;

% find the peak frequencies in x, y, and z directions; i.e., find the
% max in each direction of the normalized sum A2.
% save these variables as A3, A4, and A5

%find the indices of the peak frequencies
[maxFreq, maxIndex] = max(average(:));
[KxMaxIdx, KyMaxIdx, KzMaxIdx] = ind2sub(size(average), maxIndex);
%find the coordinates of the peak frequencies
KxMax = k(KxMaxIdx);
KyMax = k(KyMaxIdx);
KzMax = k(KzMaxIdx);

A3 = KxMax
A4 = KyMax
A5 = KzMax

%create an appropriate Gaussian filter and save it as A6

filter = g(k, [KxMax, KyMax, KzMax]);

A6 = filter;

% Using the peak frequencies for the filtered signal, estimate the x, y, and z coordinates of the Kraken over time and save as A7, A8, A9

threshold = 0.7;
xLocation = zeros(49, 1);
yLocation = zeros(49, 1);
zLocation = zeros(49, 1);
for j = 1:49
    Un(:, :, :) = reshape(Kraken(:, j), n, n, n); % We need to reshape our data into a tensor, which represents a cube of Fourier modes in x-y-z space
    fUn = fftn(Un);
    FilteredfUn = filter .* fUn;
    KrakenLocation = ifftn(FilteredfUn);
    [maxSpot, maxSpotIndex] = max(KrakenLocation(:));
    [xMaxIdx, yMaxIdx, zMaxIdx] = ind2sub(size(KrakenLocation),maxSpotIndex);
    xLocation(j) = x(xMaxIdx);
    yLocation(j) = y(yMaxIdx);
    zLocation(j) = z(zMaxIdx);
    %{
    %plotting stuff
    fM = max(abs(FilteredfUn), [], 'all');
    MLoc = max(abs(KrakenLocation), [], 'all');
    %close all, isosurface(X, Y, Z, abs(FilteredfUn)/fM, threshold)
    close all, isosurface(X, Y, Z, abs(KrakenLocation)/MLoc, threshold)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis([min(x) max(x) min(y) max(y) min(z) max(z)]), grid on, drawnow
    pause
    %} 
end

A7 = xLocation';
A8 = yLocation';
A9 = zLocation';

% Plot the location in x-y-z space over time for your report (not for the autograder)

%plot3(xLocation, yLocation, zLocation)

% Plot the projection onto the x-y plane for your reprot (not for the autograder)

%plot(xLocation, yLocation)

% Save a table of the x-y-z coordinates for your report (not for the autograder)

xyzCoords = [xLocation, yLocation, zLocation];
xyzTable = array2table(xyzCoords);
xyzTable.Properties.VariableNames = ["x", "y", "z"];
disp(xyzTable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include all helper functions below since the autograder can only accept
% one file.

function res = g(k, center)

    res = ones([length(k), length(k), length(k)]);

    for a = 1:length(k)
        for b = 1:length(k)
            for c = 1:length(k)
                x = k(a) - center(1);
                y = k(b) - center(2);
                z = k(c) - center(3);
                res(a, b, c) = exp((-1 / 10) * (x^2 + y^2 + z^2));
            end
        end
    end
end