function [s] = STEMlat01_LM(s,r)
% Refine lattice positions of the A and B sites in the HAADF image obtained
% using RealspaceLattice01. Intensity distrubution fit with a single 2D
% elliptical Gaussian.

% Code written by Greg Stone (gregastone@gmail.com).
% Code last revised by Greg Stone on August 23, 2015.
% Program based on upon code provided by and discussions with Colin Ophus 
% (clophus@lbl.gov) at NCEM.

% posRefine output: [A0 I0 x0 y0 sigma_x sigma_y theta]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fit_acc = 0; % If fit_acc = 1 the iteration differences for random peaks will be displayed in the command window.
imsize = size(s.image);
% Fitting parameters
iterTol = 0.005;%iteration tolerate
conTol = 0.1;
NsubPxIter = 32;%number of iterations
rad = 1;
%rad = 0.825; %by Greg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup image area to extract for fitting
if nargin ==1
rCut = round(rad*min([norm(s.lat(2,:))/s.latSiteFrac(1) ...
    norm(s.lat(3,:))/s.latSiteFrac(2)]));
else
rCut = r;
end
%s.lat(2,:) is the u vector, (3,:) isthe v vector; 
%latSiteFrac is the number of vectors in the unit cell.
rFit = rCut/3;
%calculate the radius of the atom for gaussian fitting
xyTol = 0.25*rCut;
sigmaMax = 0.75*rCut; sigmaMin = 0.1*rCut;
sigma0 = 0.25*rCut;
dxy = 0.5;

% Image size and (x,y) coordinates
img = s.image;
imgN = size(img);
[ya,xa] = meshgrid(1:imgN(2),1:imgN(1));

%% delete points on the edge of the image
id1 = s.pos(:,1) <= rCut/2;
id2 = s.pos(:,2) <= rCut/2;
id3 = s.pos(:,1) >= imgN(1) - rCut/2;
id4 = s.pos(:,2) >= imgN(2) - rCut/2;
s.pos(id1,1) = NaN;
s.pos(id2,1) = NaN;
s.pos(id3,1) = NaN;
s.pos(id4,1) = NaN;
s.pos(isnan(s.pos(:,1)),:)=[];
s.rCut = rCut;


% creates the coordinates for each pixel on the image. 
% Get information regarding the A and B site coordinates determine by
% RealspaceLattice01
xp = s.pos(:,1);
yp = s.pos(:,2);

Np = length(xp);
con = false(Np,1); %generates zeros of size Np*1
% Define the fitting function and options

%TolX is the lower bound of the size of a step of x direction. MaxFunEvals is the maximum
%value of function value allowed,TolFun is the minimum size of a step of y
%direction

% Creates matrix that will store the peak fitting values - GS

% posRefineM = zeros(Np,7);

% Step through each peak and use fitting to refine its position
hbar = parfor_progressbar(Np,'Please wait...');
parfor ii = 1:Np
    
%     options = optimset('TolX',1e-5,'MaxFunEvals',1e4,...
%     'TolFun',10^-5);
  options = optimset('TolX',1e-10,'MaxFunEvals',1e7,...
  'TolFun',10^-10,'Display','off');

    x0 = xp(ii);%Gets coordinates of x and y postion 
    y0 = yp(ii);
    
    x = x0; y = y0;
    sigmax = sigma0;%fitting sigmax and sigmay for ellipse shape
    sigmay = sigma0;
    
    xc = round(x);%get integer
    yc = round(y);
    
    imgN = size(img);
    % Cut out image subsection around the peak
    xv = max((xc-rCut),1):min((xc+rCut),imgN(1));
    yv = max((yc-rCut),1):min((yc+rCut),imgN(2));
    xcut = xa(xv,yv);
    ycut = ya(xv,yv);    
    % A or B sites to do the gaussian fit
    
    Icut = double(img(xv,yv));%get intensity on image with 2 precision
    k = double(min(Icut(:)));%minimum value of Icut
    I0 = double(max(Icut(:))-k);%get relative value of intensity
    jj = 1;%while loop parameter
    
    % Movement of the fitted position, if less than fititerdif, fitting
    % will stop
    fitIterDif = 0.5;
    
    % While loops cuts down on the number of fittings performed if it
    % reaches the defined tolerance threshold between susccesive fits.
    while (jj <= NsubPxIter && fitIterDif >= iterTol)
%     while (jj <= NsubPxIter)
        fun = @(c,x) c(1) + c(2)*exp(-((cos(c(7)).^2/c(5).^2+sin(c(7)).^2/c(6).^2)*(x(:,1)-c(3)).^2/2 + (-sin(2*c(7))/c(5).^2+sin(2*c(7))/c(6).^2)*(x(:,1)-c(3)).*(x(:,2)-c(4))/4 + (sin(c(7)).^2/c(5).^2+cos(c(7)).^2/c(6).^2)*(x(:,2)-c(4)).^2/2));


        xi = x; yi = y;
        sub = (xcut-x).^2 + (ycut-y).^2  < rFit^2;
        %if xcut and ycut is inside the radius of rFit. it generates the
        %matrix of 0 or 1 as variable sub.
        
        xfit = xcut(sub);%x coordinate of fitting
        yfit = ycut(sub);%y coordinate of fitting
        zfit = Icut(sub);%intensity of that position
        
        % Initial guesses
        c0 = [k I0 x y sigmax sigmay 0];
        dI = abs(I0)*.2;%I0 is the difference of min and max in Intensity
        
        lb = [k-dI I0-dI x0-xyTol y0-xyTol max(.8*sigmax,sigmaMin) max(.8*sigmay,sigmaMin) -pi/4];% this two paramters I don't understand
        ub = [k+dI I0+dI x0+xyTol y0+xyTol min(1.2*sigmax,sigmaMax) min(1.2*sigmay,sigmaMax) pi/4.01];
        
        % Perform fitting
        [cc,~,~,~] = lsqcurvefit(fun,c0,[xfit yfit],zfit,lb,ub,options);%Apply fitting function 
        %the output parameter
        k = cc(1,1);
        I0 = cc(1,2);
        x = cc(1,3);
        y = cc(1,4);
        sigmax = cc(1,5);
        sigmay = cc(1,6);
        theta = cc(1,7);
        
        fitIterDif = norm([x y] - [xi yi]);
        kk = [k I0 x y sigmax sigmay theta];
        jj = jj + 1;
        
    end
    
    

    posRefineM(ii,:) = kk;
    %it generates a matrix of from ii=1:Np, each fitted results are
    %included. 
    startP = [xp(ii) yp(ii)];
    priorP = [xi yi];
    endP = [x y];
    
    
    if fit_acc == 1 && rand(1)*50 < 1
        fitText = ['a0: ', num2str(ii), ', Iter. #: ', num2str(jj), ...
            ', Iter. dif: ', num2str(fitIterDif)];
        disp(fitText);
    end
    
    if norm(endP - priorP) < conTol && ...
            norm(startP - endP) < xyTol
        con(ii) = true;
    else
        con(ii) = false;
    end
    hbar.iterate(1); % update progress by one iteration 

end
close(hbar); % close the progress bar





%This part generates the absolute position of atoms
abPos(:,1:2) = s.pos(:,3:4);%1 and 2 column is x y coordinates
abPos(:,3) = posRefineM(:,3) - s.lat(1,1);%s.lat(1,1) is point chosen to be x coordinate of origin
abPos(:,4) = posRefineM(:,4) - s.lat(1,2);%same as (1,2), it is the y coordinate of origin, substract to get relative position
abPos(:,3:4) = abPos(:,3:4)/s.lat(2:3,:);%2:3 is the vector length of u and v
abPos(:,5) = abPos(:,3) - abPos(:,1);%difference of original position and relative position
abPos(:,6) = abPos(:,4) - abPos(:,2);
[abPos(:,7),abPos(:,8)] = cart2pol(abPos(:,5),abPos(:,6));%transfer for cartesian to polar coordinates





s.posRefineM = [posRefineM s.pos(:,3:4)];
s.abPos = abPos;% define output variables
s.fit_ID = con;
poorFit = length(con) - sum(con);%define poor fits
ij = 1:length(con);
poorFitn = ij(~con);
outText = ['Number of nonconverged fits: ', num2str(poorFit)];
disp(outText)

if poorFit > 0
   disp(num2str(poorFitn)); 
end
toc;
beep;
end
