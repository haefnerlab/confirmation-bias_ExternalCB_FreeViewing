function genCube(is_template,contrast,wPtr,w,xc,yc,category)
diff = 100;
shft = w/3;
gray = [127 127 127];
white = [255 255 255];
diff_color = 255-127;

opp_cat = 383 - 256 * contrast;
same_cat = 255;
% opp_cat = 256*contrast-1;
% same_cat = 255;
ln_width = 7;
if is_template==1
    Screen('DrawLine', wPtr, opp_cat, xc-2*w/3-diff + shft, yc, xc-2*w/3-diff + w/2 + shft, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, opp_cat, xc-2*w/3-diff + shft, yc, xc-2*w/3-diff + w/2 + shft, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, opp_cat, xc-2*w/3-diff , yc, xc-2*w/3-diff + shft, yc,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3-diff, yc, xc-2*w/3-diff + w/2, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-diff, yc, xc+w/3-w/2-diff, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-diff, yc, xc+w/3-w/2-diff, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3-diff, yc, xc-2*w/3-diff + w/2, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-diff + shft, yc, xc+w/3-w/2-diff + shft, yc+w/2 ,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-diff + shft, yc, xc+w/3-w/2-diff + shft, yc-w/2 ,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3-diff + w/2, yc+w/2, xc-2*w/3-diff + w/2 + shft, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-diff, yc, xc+w/3-diff + shft, yc,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-w/2-diff, yc-w/2, xc+w/3-w/2-diff + shft, yc-w/2,ln_width);
    
    
    Screen('DrawLine', wPtr, opp_cat, xc+w/3+diff, yc, xc+w/3-w/2+diff, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, opp_cat, xc+w/3+diff, yc, xc+w/3-w/2+diff, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, opp_cat, xc+w/3+diff, yc, xc+w/3+diff + shft, yc,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3+diff, yc, xc-2*w/3+diff + w/2, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3+diff, yc, xc-2*w/3+diff + w/2, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3+diff + shft, yc, xc-2*w/3+diff + w/2 + shft, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3+diff + shft, yc, xc+w/3-w/2+diff + shft, yc+w/2 ,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3+diff + shft, yc, xc+w/3-w/2+diff + shft, yc-w/2 ,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3+diff + shft, yc, xc-2*w/3+diff + w/2 + shft, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3+diff + w/2, yc+w/2, xc-2*w/3+diff + w/2 + shft, yc+w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc+w/3-w/2+diff, yc-w/2, xc+w/3-w/2+diff + shft, yc-w/2,ln_width);
    Screen('DrawLine', wPtr, same_cat, xc-2*w/3+diff , yc, xc-2*w/3+diff + shft, yc,ln_width);
    
    
else
    if category==-1
        Screen('DrawLine', wPtr, opp_cat, xc-2*w/3 + shft, yc, xc-2*w/3 + w/2 + shft, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, opp_cat, xc-2*w/3 + shft, yc, xc-2*w/3 + w/2 + shft, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, opp_cat, xc-2*w/3 , yc, xc-2*w/3 + shft, yc,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3, yc, xc-2*w/3 + w/2, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3, yc, xc+w/3-w/2, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3, yc, xc+w/3-w/2, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3, yc, xc-2*w/3 + w/2, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3 + shft, yc, xc+w/3-w/2 + shft, yc+w/2 ,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3 + shft, yc, xc+w/3-w/2 + shft, yc-w/2 ,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3 + w/2, yc+w/2, xc-2*w/3 + w/2 + shft, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3, yc, xc+w/3 + shft, yc,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3-w/2, yc-w/2, xc+w/3-w/2 + shft, yc-w/2,ln_width);
    else
        Screen('DrawLine', wPtr, opp_cat, xc+w/3, yc, xc+w/3-w/2, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, opp_cat, xc+w/3, yc, xc+w/3-w/2, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, opp_cat, xc+w/3, yc, xc+w/3 + shft, yc,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3, yc, xc-2*w/3 + w/2, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3, yc, xc-2*w/3 + w/2, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3 + shft, yc, xc-2*w/3 + w/2 + shft, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3 + shft, yc, xc+w/3-w/2 + shft, yc+w/2 ,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3 + shft, yc, xc+w/3-w/2 + shft, yc-w/2 ,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3 + shft, yc, xc-2*w/3 + w/2 + shft, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3 + w/2, yc+w/2, xc-2*w/3 + w/2 + shft, yc+w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc+w/3-w/2, yc-w/2, xc+w/3-w/2 + shft, yc-w/2,ln_width);
        Screen('DrawLine', wPtr, same_cat, xc-2*w/3 , yc, xc-2*w/3 + shft, yc,ln_width);
    end
end

end