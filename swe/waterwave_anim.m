% animate result from generate_waterwave
%
% copied from waterwave_orig.m
%
%   input:
%       Y   :   4d array from waterwave2
%       dt  :   time diference between samples
%       fps :   frames pre seccond 
function waterwave_anim(Y,dt,fps,file)
    [surfplot,top] = waterwave_initgraphics(Y);
    rep = size(Y,4);
    set(gcf,'color','w'); % set figure background to white
    drawnow;
        
    for nstep = 1:rep
        C = squeeze(abs(Y(:,:,2,nstep)) + abs(Y(:,:,3,nstep)));  % Color shows momemtum
        t = nstep*dt;
        set(surfplot,'zdata',squeeze(Y(:,:,1,nstep)),'cdata',C);
        set(top,'string',sprintf('t = %6.2f',t))
        drawnow;

        if(1)
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
    
            if nstep == 1 
                % On the first loop, create the file. In subsequent loops, append
                imwrite(imind,cm,file,'gif','DelayTime',0,'loopcount',inf);
            else
                imwrite(imind,cm,file,'gif','DelayTime',0,'writemode','append');
            end
        pause(1/fps);
        end
    end
end
      