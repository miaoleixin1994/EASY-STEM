function Point_moving_origin(src,evt)
evname = evt.EventName;
switch(evname)
    case{'MovingROI'}
        refreshdata;
    case{'ROIMoved'}
        assignin('base','Origin_Pos',(evt.CurrentPosition));
        atom_pos = evalin('base','atom_pos');
        Origin_Pos = evalin('base','Origin_Pos'); 

        cursor_point = [Origin_Pos(2) Origin_Pos(1)];  %Extract Coordinates of the point ROI
        Idx = knnsearch(atom_pos.pos(:,1:2),cursor_point);
        Origin_Pos = fliplr(atom_pos.pos(Idx,1:2));       
        
        lat = evalin('base','lat');
        LatNumPlot = evalin('base','LatNumPlot');
        hold on;
        imshow(atom_pos.image);scatter(atom_pos.pos(:,2),atom_pos.pos(:,1),'y.');
        q1 = quiver(Origin_Pos(1),Origin_Pos(2),lat(2,2)*LatNumPlot,lat(2,1)*LatNumPlot,'LineWidth',4,'Color',[1 0.2 0],'MaxHeadSize',0.3,'AutoScale',"off");
        q2 = quiver(Origin_Pos(1),Origin_Pos(2),lat(3,2)*LatNumPlot,lat(3,1)*LatNumPlot,'LineWidth',4,'Color',[0 0.2 1],'MaxHeadSize',0.3,'AutoScale',"off");
        
        q1.XDataSource = 'Origin_Pos(1)';
        q1.YDataSource = 'Origin_Pos(2)';
        q2.XDataSource = 'Origin_Pos(1)';
        q2.YDataSource = 'Origin_Pos(2)';
        assignin('base','if_moved',1);
        
        lat(1,2) = Origin_Pos(1);
        lat(1,1) = Origin_Pos(2);
        assignin('base','lat',lat);
        linkdata on
end

end
