function [ Y_color ] = ColorMap_ScaleBar(min_ex,max_ex,mean_ex,cmap )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%     bottom = [0 0.1 1];
%     mid_bot = [0 0.5 1];
%     middle = [1 1 1];
%     mid_top = [1 0.5 0];
%     top = [1 0 0];
    scale = 256;
    aa = 1:scale;
    a = meshgrid(aa,1:30);
    Quantized_y = uint8(a);

    if nargin == 0
    bottom = [1 0 0];
    mid_bot = [];
    middle = [1 1 1];
    mid_top = [];
    top = [0 0 1];        
    
    colorscale = scale;
    new = [bottom;mid_bot;middle;mid_top;top];
    ss = size(new);
    len = ss(1);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, colorscale);
    newmap = zeros(colorscale, 3);
    
    
    
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end


    Y_color = zeros(30,scale,3);

    for ii = 1:30
        for jj = 1:scale
            Y_color(ii,jj,1:3) = newmap(Quantized_y(ii,jj),1:3);
        end
    end
    
    
    figure; imagesc(Y_color);axis image off;
    
    
    else
        scale = 256;
        width = 32;
        aa = 1:scale;
        a = meshgrid(aa,1:width);
        Quantized_y = uint8(a);
        
        newmap = cmap;
        Y_color = zeros(width,scale,3); 
    
        for ii = 1:width
            for jj = 1:scale
                Y_color(ii,jj,1:3) = newmap(Quantized_y(ii,jj),1:3);
            end
        end
    
    Y_color = permute(Y_color,[2 1 3]);
    Y_color = flipud(Y_color);
    figure; imagesc((Y_color));axis image off;hold on;...
%     figure; imagesc((Y_color));axis image off;hold on;...

    text(15,128,[num2str(mean_ex,2) 'pm'],'Color','black','FontWeight','bold','HorizontalAlignment','center','FontSize',15);    
    text(15,250,[num2str(min_ex,2) 'pm'],'Color','black','FontWeight','bold','HorizontalAlignment','center','FontSize',15);
    text(15,6,['+' num2str(max_ex,2) 'pm'],'Color','black','FontWeight','bold','HorizontalAlignment','center','FontSize',15);
    hold off;
    
    figure; imagesc((Y_color));axis image off;
%     
%     text(15,250,[num2str(min_ex - mean_ex,3) 'pm'],'Color','black','FontWeight','bold','HorizontalAlignment','center','FontSize',15);
%     text(15,6,['+' num2str(max_ex - mean_ex,3) 'pm'],'Color','black','FontWeight','bold','HorizontalAlignment','center','FontSize',15);
    end

    
end

