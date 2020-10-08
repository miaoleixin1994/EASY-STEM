function [s] = RealSpaceStrain(s,strain_range)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    idxy = mod(s.pos(:,3),1) == 0 & mod(s.pos(:,4),1) == 0; 
    
    xp0 = s.posRefineM(idxy,3);
    yp0 = s.posRefineM(idxy,4);
  
    up = s.pos(idxy,3);
    vp = s.pos(idxy,4);
    
    ui = zeros(length(xp0),2);
    vi = ui;
    uc = zeros(length(xp0),2);
    
    nn = 1;
    for ii = 1:length(up)
        
        try
        id1 = s.pos(:,3) == (up(ii,1) + nn) & s.pos(:,4) == (vp(ii,1));
        xp1 = s.posRefineM(id1,3);
        yp1 = s.posRefineM(id1,4);
        
        id2 = s.pos(:,3) == (up(ii,1)) & s.pos(:,4) == (vp(ii,1) + nn);
        xp2 = s.posRefineM(id2,3);
        yp2 = s.posRefineM(id2,4);        
        
        id3 = s.pos(:,3) == (up(ii,1) + nn) & s.pos(:,4) == (vp(ii,1)+nn);
        xp3 = s.posRefineM(id3,3);
        yp3 = s.posRefineM(id3,4);         
        
        
        ui(ii,1:2) = - [xp0(ii) yp0(ii)] + [xp1 yp1];
        vi(ii,1:2) = - [xp0(ii) yp0(ii)] + [xp2 yp2];
        
        uc(ii,1) = (xp0(ii) + xp1 + xp2 + xp3)/4;             
        uc(ii,2) = (yp0(ii) + yp1 + yp2 + yp3)/4;
        
        
        catch
%             ui(ii,:) = NaN;
%             vi(ii,:) = NaN;
%             uc(ii,:) = NaN;
        end
        
    end
    
        ifnonzero1 = all(ui,2);
        ifnonzero2 = all(vi,2);
        ifnonzero = logical(ifnonzero1.*ifnonzero2);
        ui(~ifnonzero,1) = NaN;
        ui(isnan(ui(:,1)),:)=[];
        vi(~ifnonzero,1) = NaN;
        vi(isnan(vi(:,1)),:)=[];
        uc(~ifnonzero,1) = NaN;
        uc(isnan(uc(:,1)),:)=[];        
        Imat = eye(2,2);
        epsilonMat = zeros(2,2,length(ui));
        omegaMat = epsilonMat;
        uvmat = zeros(2,2,length(ui)); 
        
        
        for ii = 1:length(ui)
            ref_uv = nn*s.lat(2:3,:);
            uvmat(:,:,ii) = [ui(ii,:);vi(ii,:)];
            T = uvmat(:,:,ii)*inv(ref_uv);
            D = T - Imat;
            epsilonMat(:,:,ii) = (D+transpose(D))/2;
            omegaMat(:,:,ii) = (D-transpose(D))/2;        
        end
        s.uc = uc;
        s.ui = ui;
        s.vi = vi;
        s.epsilonMat = epsilonMat;
        s.omegaMat = omegaMat;
        s.uvmat = uvmat;
        
        
        imsize = size(s.image);
        [xi,yi] = meshgrid(1:1:imsize(1,2),1:1:imsize(1,1));
        YList(:) = uc(:,1);
        XList(:) = uc(:,2);
        espxx_list = squeeze(epsilonMat(1,1,:));
        espyy_list = squeeze(epsilonMat(2,2,:));        
        s.epsxx = griddata(YList,XList,espxx_list,xi,yi,'v4');
        s.epsyy = griddata(YList,XList,espyy_list,xi,yi,'v4');
        newmap = bluewhitered_v3(256);
        up_lim = strain_range;low_lim = -strain_range;
        eps_xx_colorscale = ((s.epsxx - low_lim)/(up_lim - low_lim));
        eps_xx_colorscale = ceil(eps_xx_colorscale*256);
        eps_xx_colorscale = 1+uint8(eps_xx_colorscale);
        s.epsxx_color = zeros(imsize(1),imsize(2),3);
        
        for ii = 1:imsize(1)
            for jj = 1:imsize(2)
        %         u_color(ii,jj,1:3) = newmap(Quantized_U(ii,jj),1:3);
                s.epsxx_color(ii,jj,1:3) = newmap(eps_xx_colorscale(ii,jj),1:3);
            end
        end
        
        eps_yy_colorscale = ((s.epsyy - low_lim)/(up_lim - low_lim));
        eps_yy_colorscale = ceil(eps_yy_colorscale*256);
        eps_yy_colorscale = 1+uint8(eps_yy_colorscale);
        s.epsyy_color = zeros(imsize(1),imsize(2),3);
        for ii = 1:imsize(1)
            for jj = 1:imsize(2)
        %         u_color(ii,jj,1:3) = newmap(Quantized_U(ii,jj),1:3);
                s.epsyy_color(ii,jj,1:3) = newmap(eps_yy_colorscale(ii,jj),1:3);
            end
        end        
        
        figure;imagesc(s.epsxx);axis image off;
        figure;imshow(s.epsxx_color);
        figure;imshow(s.epsyy_color);
end

