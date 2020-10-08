function [ IFFT_Image,Fourier_Pattern,Filtered_Pattern,FFT_Filter] = Fourier_Filter_TEM(serfile,r,cali)
%Writen by Leixin Miao to do Fourier Filtering on TEM images
%   Detailed explanation goes here

%% Input Parameter


if nargin == 3
    RealSpaceCalibration = cali;
    RealSpaceCalibration = (10^9)*RealSpaceCalibration;
    original_image = serfile;
    ImageSize = size(original_image);
else
    original_image = serfile.image;
    RealSpaceCalibration = serfile.calibration(1);
    ImageSize = size(original_image);
    RealSpaceCalibration = (10^9)*RealSpaceCalibration;
end
    


if nargin == 1
prompt = 'What is the resolution value in pm? ';
resolution = input(prompt);
resolution = resolution/1000;
else
    resolution = r/1000;
end


%% Calculating the Fourier Space Pixel Size
Nx = ImageSize(1);
Lx = Nx*RealSpaceCalibration; %RealSpaceCalibration is in nm
MomCalib = ((-Nx/2):(Nx/2-1))/Lx;
MomCalib = mean(diff(MomCalib));% pixel size in Fourier Space in nm-1

%% Creating the filter
r = 1/resolution;
r = r/MomCalib;
radius = ceil(r);


if radius < Nx
qxx = circshift(((-Nx/2+1):(Nx/2)),[Nx/2 0]);
qyy = qxx;
[qyya, qxxa] = meshgrid(qyy,qxx);
aa(:,:,1) = qxxa;aa(:,:,2)=qyya;
Circle = aa(:,:,1).^2+aa(:,:,2).^2;
FFT_Filter = ones(ImageSize(1),ImageSize(2));
FFT_Filter(Circle>((radius/2)^2)) = 0;
FFT_Filter = imgaussfilt(FFT_Filter,round(radius*0.02));
FFT_Filter = imgaussfilt(FFT_Filter,1);
%% Filtering the FFT
ff = fftshift(fft2(fftshift(original_image)));
Filtered_ff = ff.*FFT_Filter;
IFFT_image = ifftshift(ifft2(ifftshift(Filtered_ff)));
IFFT_Image = ImNorm(abs(IFFT_image));
Fourier_Pattern = (abs(ff));
Filtered_Pattern = (abs(Filtered_ff));
Fourier_Difference = Fourier_Pattern - Filtered_Pattern;

if nargin ~= 2
figure;
subplot(1,3,1);imagesc((log(Fourier_Pattern)));title('Fourier Pattern');axis image off;colormap(jet);
subplot(1,3,2);imagesc((log(Filtered_Pattern)));axis image off;colormap(jet);title('Filtered Pattern');
subplot(1,3,3);imagesc(((Fourier_Difference)));axis image off;colormap(jet);title('Difference');

figure;
subplot(1,2,1);imshow(ImNorm(original_image));title('original Image');subplot(1,2,2);imshow(IFFT_Image);title('IFFT Image');
end


else
    IFFT_Image = original_image;
    Fourier_Pattern = [];
    Filtered_Pattern = [];
    FFT_Filter = [];
end

%% Plot

end

