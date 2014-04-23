% animate result from generate_waterwave
%
% copied from waterwave_orig.m
function waterwave_anim(Y,dt)
    [surfplot,top,start,stop] = waterwave_initgraphics(Y);
    rep = size(Y,4);
    nstep = 1;
    while get(start,'value')==0 && get(stop,'value')==0 && nstep <= rep
        C = squeeze(abs(Y(:,:,2,nstep)) + abs(Y(:,:,3,nstep)));  % Color shows momemtum
        t = nstep*dt;
        set(surfplot,'zdata',squeeze(Y(:,:,1,nstep)),'cdata',C);
        set(top,'string',sprintf('t = %6.2f',t))
        drawnow
        nstep = nstep+1;
        pause(dt);
    end
end
      