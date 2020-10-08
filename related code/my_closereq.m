function my_closereq(app,src,callbackdata)
% Close request function 
% to display a question dialog box 
   selection = questdlg('Stop the Finding Origin and Close This Figure?',...
      'Close Request Function',...
      'Yes','No','Yes'); 
   switch selection 
      case 'Yes'
%          set(gcf,'UserData',0)
         delete(gcf)
         
         
      case 'No'
      return 
   end
end
