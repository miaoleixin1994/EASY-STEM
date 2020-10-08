function sdata = MPFit(sdata,saveAll)
%% Explanation
% Uses an iterative multi-step fitting process to represent each atom
% position in the aberration corrected STEM image as an analytical sum of 
% 2D Gaussians. Stores all the Gaussians too - generate later as needed. 
% Code originally developed by Debangshu Mukherjee on September 18th, 2017 - 
% inspired by STEMlat01_gs and STEMlat01m_dm. sdata is the struct file. 
% sdata.image must be the image file, while sdata.pos lists the intensity 
% minima or the intensity maxima. Code modified by Debangshu Mukherjee on 
% March 25th, 2018. Email: debangshu24@gmail.com

%% User Defined Parameters
tic;
%% Start parallel pool
if isempty(gcp('nocreate'))
    clust = parcluster;
    parpool((clust.NumWorkers - 1));
end

if nargin < 2
    saveAll = 0;
end
tol = 1e-12; %residual error tolerance

%% Function Parameters
rsl_positions = sdata.pos; %Initial positions corresponding to brightest 
                       %or darkest pixels from RealspaceLattice01
NumberAtoms = length(rsl_positions);
Img = sdata.image - min(sdata.image(:)); %The STEM image
Img = Img/(max(Img(:)));
immedian = median(Img(:));
correction1 = ceil(mean(ceil(rsl_positions(:,1)) - floor(rsl_positions(:,1))));

%% Get nearest neighbors and calculation size
nearestNeighbors = 6;
iterCount = 2*(2+nearestNeighbors);
neighborPosList = zeros(NumberAtoms,(2*nearestNeighbors));
initial_positions = rsl_positions(:,1:2);
nearest_idx1 = knnsearch(initial_positions,initial_positions,'K',(nearestNeighbors+1));
for ii = 1:nearestNeighbors
    neighborPosList(:,(2*ii - 1):(2*ii)) = initial_positions(nearest_idx1(:,(ii+1)),1:2);
end
clearvars ii
neighbor_distance = zeros(NumberAtoms,nearestNeighbors);
for jj = 1:nearestNeighbors
    ydist = initial_positions(:,1) - neighborPosList(:,(2*jj - 1));
    xdist = initial_positions(:,2) - neighborPosList(:,(2*jj));
    neighbor_distance(:,jj) = ((ydist.^2) + (xdist.^2)).^0.5;
end
clearvars jj ydist xdist
minPeakDistances = min(neighbor_distance,[],2);
CutPeakSize = ceil(0.5*(median(neighbor_distance(:))));
InterPeakDistance = (2*CutPeakSize) - 1 + correction1;
padVal = ceil(1.75*InterPeakDistance);
PaddedImage = padarray(Img,[padVal padVal],immedian,'both');
ImageSize = size(PaddedImage);

%% Output Functions
GaussList1 = zeros(NumberAtoms,iterCount,7);
[xV,yV] = meshgrid(1:ImageSize(2),1:ImageSize(1));
LengthSub1 = InterPeakDistance^2;
ListOfFits1 = zeros(NumberAtoms,LengthSub1); %Fitted Peaks
GaussEvol1 = zeros(NumberAtoms,LengthSub1,iterCount);
ListOfZ1 = ListOfFits1; %Orignal Z values
ListOfTails1 = ListOfFits1; %Tail cutoffs
fitsize1 = zeros(NumberAtoms,1);
DistList1 = zeros(NumberAtoms, iterCount);
xPosStep1 = initial_positions(:,2);
yPosStep1 = initial_positions(:,1);

%% Iterative Peak Fit Step 1
dq1 = parallel.pool.DataQueue;
N = iterCount*NumberAtoms;
waitbar1 = waitbar(0, 'Fitting Multiple Gaussians...');
% Use the waitbar's UserData to track progress
waitbar1.UserData = [0 N];
afterEach(dq1, @(varargin) iIncrementWaitbar(waitbar1));
parfor PeakPos = 1:NumberAtoms
    yPos = initial_positions(PeakPos,1) + padVal;
    xPos = initial_positions(PeakPos,2) + padVal;
    minNeighbor = minPeakDistances(PeakPos);
    sub = (abs(yV - yPos) < CutPeakSize) & (abs(xV - xPos) < CutPeakSize);
    xVals = xV(sub);
	yVals = yV(sub);
	zVals1 = PaddedImage(sub);
    ListOfZ1(PeakPos,:) = zVals1;
    fitr = zeros(iterCount,7);
    realfits = zeros(iterCount,7);
    fitpeak = [0 0];
    fitp = zeros(iterCount,2);
    FittedZ = zeros(LengthSub1,1);
    FittedRest = zeros(LengthSub1,1);
    fitsize1(PeakPos) = length(zVals1);
    for iter_number=1:iterCount
        if (median(abs(zVals1)) > tol)
            [fitresult, zfit] = Custom_Gaussian(xVals,yVals,zVals1,(0.25*InterPeakDistance));
            fitr(iter_number,:) = fitresult;
            distRSL = (((fitresult(6) - yPos)^2) + ((fitresult(5) - xPos)^2))^0.5;
            GaussList1(PeakPos,iter_number,:) = fitresult;
            zVals1 = zVals1 - zfit;
            DistList1(PeakPos,iter_number) = distRSL;
            if (distRSL > (0.25*minNeighbor))
                FittedRest = FittedRest + zfit;
                realfits(iter_number,:) = zeros(1,7);
            else
                FittedZ = FittedZ + zfit;
                realfits(iter_number,:) = fitresult;
            end
            GaussEvol1(PeakPos,:,iter_number) = FittedZ + FittedRest;
        else
            fitr(iter_number,:) = 0;
        end
        send(dq1,iter_number);
    end
    ListOfFits1(PeakPos,:) = FittedZ;
    ListOfTails1(PeakPos,:) = FittedRest;
    
    fitp(:,1) = realfits(:,7).*realfits(:,5);
    fitp(:,2) = realfits(:,7).*realfits(:,6);
    fitpeak(1,1) = (sum(fitp(:,1)))/(sum(realfits(:,7))) - padVal;
    fitpeak(1,2) = (sum(fitp(:,2)))/(sum(realfits(:,7))) - padVal;
    xPosStep1(PeakPos,1) = fitpeak(1,1);
    yPosStep1(PeakPos,1) = fitpeak(1,2);
end
xcorrect1 = (isnan(xPosStep1)).*initial_positions(:,2);
ycorrect1 = (isnan(yPosStep1)).*initial_positions(:,1);
xPosStep1(isnan(xPosStep1)) = 0.001;
yPosStep1(isnan(yPosStep1)) = 0.001;
xPosStep1 = xPosStep1 + xcorrect1;
yPosStep1 = yPosStep1 + ycorrect1;
close(waitbar1);

%% Save the calculated variables
nearest_idx1(:,1) = [];
posRefine(:,3) = yPosStep1;
posRefine(:,4) = xPosStep1;
posRefine(:,1:2) = rsl_positions(:,1:2);
sdata.posRefine = posRefine;
sdata.nearest_neighbors1 = nearest_idx1;
sdata.neighbor_distance1 = neighbor_distance;
sdata.rFit = InterPeakDistance;
sdata.Img = Img;
sdata.iterCount = iterCount;
sdata.padVal = padVal;
if saveAll > 0
    sdata.GaussList1 = GaussList1;
    sdata.ListOfFits1 = ListOfFits1;
    sdata.fitsize1 = fitsize1;
    sdata.ListOfZ1 = ListOfZ1;
    sdata.GaussEvol1 = GaussEvol1;
    sdata.ListOfTails1 = ListOfTails1;
    sdata.DistList1 = DistList1;
end
% delete(gcp('nocreate')); %Shut down parallel pool
toc;
end

%% Gaussian Peak Refiner
function [fitresult, zfit] = Custom_Gaussian(xx,yy,zz,boundary_checker,x0,y0)
[xCleaned, yCleaned, zCleaned] = prepareSurfaceData(xx,yy,zz);
[amp, ind] = max(abs(zCleaned));
if nargin < 5
    x0 = xCleaned(ind); 
    y0 = yCleaned(ind);
end
%% Set up boundary conditions
sigmaMax = boundary_checker; 
sigmaMin = 0.05*boundary_checker;
sigma0 = 0.5*boundary_checker;
peak_range = boundary_checker;

xy_values = {xCleaned,yCleaned};

%% Set up the startpoint

ang = 45; % angle in degrees.
Sigma_Y = sigma0; %initial sigma value of Y
Sigma_X = sigma0; %initial sigma value of X
z0 = median(zCleaned(:));

%% Force the fitted gaussians to stay within a range of the original RealspaceLattice peak
xmax = x0 + peak_range;
ymax = y0 + peak_range;
xmin = x0 - peak_range;
ymin = y0 - peak_range;

%% Give the initial fitting options.
LowerBound = [-Inf, 0, sigmaMin, sigmaMin, xmin, ymin, -Inf];
UpperBound = [Inf, 180, sigmaMax, sigmaMax, xmax, ymax, Inf]; % angles greater than 90 are redundant
StartPoint = [amp, ang, Sigma_X, Sigma_Y, x0, y0, z0];%[amp, sx, sxy, sy, xo, yo, zo];

tols = 1e-8;
options = optimset('Algorithm','trust-region-reflective',...
    'Display','off',...
    'MaxFunEvals',5e2,...
    'MaxIter',5e2,...
    'TolX',tols,...
    'TolFun',tols,...
    'TolCon',tols); 

%% Perform Gaussian Peak Fitting
fitresult = lsqcurvefit(@gaussian2D,StartPoint,xy_values,zCleaned,LowerBound,UpperBound,options);

%% Calculate Peak Fit Errors
zfit = gaussian2D(fitresult,xy_values);
end

%% Gaussian Peak Function Used
function z = gaussian2D(par,xy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
end

%% Waitbar Function
function iIncrementWaitbar(wb)
ud = wb.UserData;
ud(1) = ud(1) + 1;
waitbar(ud(1) / ud(2), wb);
wb.UserData = ud;
end