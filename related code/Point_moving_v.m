function Point_moving_v(src,evt)
evname = evt.EventName;
switch(evname)
    case{'MovingROI'}
        refreshdata;
    case{'ROIMoved'}
        assignin('base','V_Pos',(evt.CurrentPosition));
        V_Pos = evalin('base','V_Pos');        
        lat = evalin('base','lat');
        atom_pos = evalin('base','atom_pos');
        LatNumPlot = evalin('base','LatNumPlot');
        cursor_point = [V_Pos(2) V_Pos(1)];  %Extract Coordinates of the point ROI
        Idx = knnsearch(atom_pos.pos(:,1:2),cursor_point);
        V_Pos = fliplr(atom_pos.pos(Idx,1:2));               
        
        
        hold on;
        u1 = V_Pos(1)-lat(1,2);
        v1 = V_Pos(2)-lat(1,1);
        assignin('base','u1',(u1));
        assignin('base','v1',(v1));
        imshow(atom_pos.image);scatter(atom_pos.pos(:,2),atom_pos.pos(:,1),'y.');
%         q1 = quiver(lat(1,2),lat(1,1),u1,v1,'LineWidth',4,'Color',[1 0.2 0],'MaxHeadSize',0.3,'AutoScale',"off");
        q1 = quiver(lat(1,2),lat(1,1),lat(2,2)*LatNumPlot,lat(2,1)*LatNumPlot,'LineWidth',4,'Color',[1 0.2 0],'MaxHeadSize',0.3,'AutoScale',"off");
        q2 = quiver(lat(1,2),lat(1,1),u1,v1,'LineWidth',4,'Color',[0 0.2 1],'MaxHeadSize',0.3,'AutoScale',"off");
        linkdata on
        
        q2.UDataSource = 'u1';
        q2.VDataSource = 'v1';
%         q2.XDataSource = 'U_Pos(1)';
%         q2.YDataSource = 'U_Pos(2)';
        assignin('base','if_moved',1);
        

        lat(3,2) = (V_Pos(1)-lat(1,2))/LatNumPlot;
        lat(3,1) = (V_Pos(2)-lat(1,1))/LatNumPlot;
        assignin('base','lat',lat);
        
end

end
