function Out=AskUser(text)
    
% %     function Yes(btn,ax)
% %         Out=1;
% %         closereq
% %         close all
% %     end
% %     function No(btn,ax)
% %         Out=0;
% %         closereq
% %         close all
% %     end
% % 
% %     fig = uifigure('Position',[300 300 400 200], 'Color', 'r');
% %     bg = uibuttongroup(fig,'Position',[100 50 200 50]);
% %     tb1 = uibutton(bg,'Position',[50 25 100 25],'ButtonPushedFcn', @(btn,event) Yes(btn));
% %     tb2 = uibutton(bg,'Position',[50 0 100 25],'ButtonPushedFcn', @(btn,event) No(btn));
% %     
% %     
% %     label1 = uilabel(fig,...
% %     'Position',[100 164 200 15],...
% %     'Text','Do you have an transformation file?');
% %     
% %     tb1.Text = 'Yes';
% %     tb2.Text = 'No';
% %     uiwait(fig)

    
    prompt =text;
    x = input(prompt);
    
    if x == 1
        Out=1;
    elseif x==0
        Out = 0;
    else
       disp(' Enter 0 or 1 ') ;
       return
    end

        
end