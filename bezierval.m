function b = bezierval(bez,t)
    % Evaluate bezier curve with parameters bez(:,1),bez(:,2),... at points t
    % size(t) = [1,N], where N is number of points at which to evaluate
    % function.  Typically, t = [0,1]
    % size(bez) = [2,order+1]
    order = size(bez,2)-1;
    switch order
        case 1
            p0 = bez(:,1);
            p1 = bez(:,2);
            b = p0 + (p1-p0)*t;
        case 2
            p0 = bez(:,1);
            p1 = bez(:,2);
            p2 = bez(:,3);
            b = p0*(1-t).^2 + 2*p1*((1-t).*t) + p2*t.^2;
        otherwise
            error('Write more code!')
    end