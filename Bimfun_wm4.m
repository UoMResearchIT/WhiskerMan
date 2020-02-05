function [E,gradE] = Bimfun_wm4(z, r, np0, np1, imh, imv, dt, sigma1, sigma2, calib)
    % sigma1 keeps solution close to the previous frame's solution
    % sigma2 keeps P1 near the middle
    %
    % see notebook 190216 for explicit equations
    
    % extension to 3d nec before debug can be used...
    % debug = 0;
    
    % column vectors that span plane normal to bezier at t==0 - r(:,1):
    na = np0(:,1);
    nb = np0(:,2);
    % column vectors that span plane normal to bezier at t==1 - r(:,3):
    nc = np1(:,1);
    nd = np1(:,2);
    
    mv = calib.matrix(1,:);
    mw = calib.matrix(2,:);
    ov = calib.vector(1);
    ow = calib.vector(2);
    
    % Bezier control points:
    P = zeros(3,3);
    P(:,1) = r(:,1) + z(1)*na + z(2)*nb;
    P(:,2) = z(3:5)';
    P(:,3) = r(:,3) + z(6)*nc + z(7)*nd;
    
    % temporal contiguity constraint:
    Ereg1 = 0.5 * sigma1 * sum(sum((P-r).^2));
    gradEreg1 = zeros(7,1);
    gradEreg1(1) = sigma1 * (P(:,1)-r(:,1))'*na;
    gradEreg1(2) = sigma1 * (P(:,1)-r(:,1))'*nb;
    gradEreg1(3:5) = sigma1 * (P(:,2)-r(:,2));
    gradEreg1(6) = sigma1 * (P(:,3)-r(:,3))'*nc;
    gradEreg1(7) = sigma1 * (P(:,3)-r(:,3))'*nd;
    
    % Constraint to centre the middle control point:
    p = P(:,3)-P(:,1);
    q = P(:,2)-P(:,1);
    Ereg2 = 0.5 * sigma2 * (q'*p/norm(p) - .5*norm(p)).^2;
    gradEreg2 = zeros(7,1);
    gradEreg2(3:5) = sigma2 * (q'*p/norm(p) - .5*norm(p)) * p;
    
    % following will need extension to 3d...
    % if debug
    %     figure, subplot 121, hold on
    %     plot(P(1,:),P(2,:),'r.')
    %     subplot 122, hold on
    %     plot(p(1,:),p(2,:),'o',q(1,:),q(2,:),'sq')
    %     legend('p','q')
    %     pause
    %     close
    % end
    clear p q
    
    clear r
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the 'posterior' energy term - ie line integral of bezier over the image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = 0:dt:1;
    b = bezierval(P,t);
    x = b(1,:);
    y = b(2,:);
    v = mv*b + ov;
    w = mw*b + ow;
    
    % find the points on the projected bezier curves that are internal to the images:
    wid = size(imh.s,2);
    hgt = size(imh.s,1);
    gdpt.h = find((x>=1)&(x<=wid)&(y>=1)&(y<=hgt));
    Eh = dt*sum(interp2(imh.s,x(gdpt.h),y(gdpt.h)));
    
    if isfield(imv,'s')
        wid = size(imv.s,2);
        hgt = size(imv.s,1);
        gdpt.v = find((v>=1)&(v<=wid)&(w>=1)&(w<=hgt));
        Ev = dt*sum(interp2(imv.s,v(gdpt.v),w(gdpt.v)));
    else
        Ev=0;
    end
    clear wid hgt
    E = Eh + Ev + Ereg1 + Ereg2;
    
    
    
    % compute the gradients of Eh and Ev wrt the 'z' pmtrs:
    % gradients of images:
    id.x = interp2(imh.dx,x(gdpt.h),y(gdpt.h));
    id.y = interp2(imh.dy,x(gdpt.h),y(gdpt.h));
    
    % d Eh/dzi:
    th = t(gdpt.h);
    gradEh = zeros(7,1);
    gradEh(1) = dt*sum(id.x.*((1-th).^2)*na(1)+id.y.*((1-th).^2)*na(2));
    gradEh(2) = dt*sum(id.x.*((1-th).^2)*nb(1)+id.y.*((1-th).^2)*nb(2));
    gradEh(3) = dt*sum(id.x.*(2*(1-th).*th));
    gradEh(4) = dt*sum(id.y.*(2*(1-th).*th));
    gradEh(5) = 0;
    gradEh(6) = dt*sum(id.x.*(th.^2)*nc(1)+id.y.*(th.^2)*nc(2));
    gradEh(7) = dt*sum(id.x.*(th.^2)*nd(1)+id.y.*(th.^2)*nd(2));
    clear th
    % d Ev/dzi:
    if isfield(imv,'dx')
        id.v = interp2(imv.dx,v(gdpt.v),w(gdpt.v));
        id.w = interp2(imv.dy,v(gdpt.v),w(gdpt.v));
        
        tv = t(gdpt.v);
        gradEv = zeros(7,1);
        gradEv(1) = dt*sum(id.v.*((1-tv).^2)*(mv*na)+id.w.*((1-tv).^2)*(mw*na));
        gradEv(2) = dt*sum(id.v.*((1-tv).^2)*(mv*nb)+id.w.*((1-tv).^2)*(mw*nb));
        gradEv(3) = dt*sum(id.v.*(2*(1-tv).*tv)*mv(1)+id.w.*(2*(1-tv).*tv)*mw(1));
        gradEv(4) = dt*sum(id.v.*(2*(1-tv).*tv)*mv(2)+id.w.*(2*(1-tv).*tv)*mw(2));
        gradEv(5) = dt*sum(id.v.*(2*(1-tv).*tv)*mv(3)+id.w.*(2*(1-tv).*tv)*mw(3));
        gradEv(6) = dt*sum(id.v.*(tv.^2)*(mv*nc)+id.w.*(tv.^2)*(mw*nc));
        gradEv(7) = dt*sum(id.v.*(tv.^2)*(mv*nd)+id.w.*(tv.^2)*(mw*nd));
        clear tv
    else
        gradEv = zeros(7,1);
    end
    clear gdpt id
    
    % put it all together:
    gradE = gradEh + gradEv + gradEreg1 + gradEreg2;