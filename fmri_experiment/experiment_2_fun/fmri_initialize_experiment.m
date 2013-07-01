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
                    [ pressed, firstPress]=KbQueueCheck(params.keyboard_device);
                    if pressed
                        if firstPress(params.trigger),
                            fprintf('\rStarting the Experiment')
                            break
                        elseif firstPress(params.bathroom_response) || firstPress(params.kitchen_response),
                            if i == numel(instructions),
                            Screen('DrawLine', params.w, params.white, params.centx-params.CrossWidth, params.centy, params.centx+params.CrossWidth, params.centy, 10);
                            Screen('DrawLine', params.w, params.white, params.centx, params.centy-params.CrossWidth, params.centx, params.centy+params.CrossWidth, 10);
                            Screen(params.w, 'Flip');
                            else
                               break;
                            end
                        elseif firstPress(params.escapeKey),
                            Screen('CloseAll');
                            ListenChar(0);
                            KbQueueRelease(params.keyboard_device);
                            break;
                        end
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
