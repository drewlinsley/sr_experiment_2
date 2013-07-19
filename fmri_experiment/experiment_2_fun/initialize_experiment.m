function status = fmri_initialize_experiment(params,instructions)


%% Set up our displays and i/o
try
    
    %% Instructions
    %%
    
    switch params.instruction
        case 'on'
            
            
            for i = 1:numel(instructions),
                
                Screen('TextFont',params.w,'Helvetica');
                Screen('TextSize',params.w,29);
                
                %Read instructions:
                fd = fopen(instructions{i},'rt');
                if fd == -1,
                    error('Problem opening instructions*.txt')
                end
                mytext = '';
                tl = fgets(fd);
                lcount = 0;
                while (tl~=-1) & (lcount < 48)
                    mytext = [mytext tl];
                    tl = fgets(fd);
                    lcount = lcount + 1;
                end
                fclose(fd);
                mytext = [mytext char(10)];
                
                [~, ~, bbox] = DrawFormattedText(params.w,mytext,'center','center',[params.white params.white params.white]);
                
                %Show computed bounding box
                Screen('FrameRect',params.w, params.white, bbox);
                Screen('Flip',params.w);
                
                while 1
                    pressed = KbCheck(-1);
                    if pressed
                            break;
                    end
                end
                
                
            end
            
            
            
            %nothing....
            
    end
    
catch status
    Screen('CloseAll');
    save('initialize_error','status')
    throw(status)
end

if ~exist('status','var'),
    status = 1;
end