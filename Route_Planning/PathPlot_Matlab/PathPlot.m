function PathPlot()

    type(1)=1;
    type(2)=3;
    type(3)=1;
    type(4)=0;
    type(5)=0;
    %type = path.type;
%     x = [];
%     y = [];
    path.t=1.21158;
    path.u=-0.9397;
    path.v=0.9903;
    path.w=0;
    path.x=0;
    seg = [path.t,path.u,path.v,path.w,path.x];
    pvec = [0,0,0];
    rmin = 5;
    for i = 1:5        
        if type(i) == 2
            theta = pvec(3);
            dl = rmin*seg(i);
            dvec = [dl*cos(theta), dl*sin(theta), 0];
            dx = pvec(1)+linspace(0,dvec(1));
            dy = pvec(2)+linspace(0,dvec(2));
%             x = [x,dx];
%             y = [y,dy];
            pvec = pvec+dvec;
        elseif type(i) == 1
            theta = pvec(3);
            dtheta = seg(i);
            cenx = pvec(1)-rmin*sin(theta);
            ceny = pvec(2)+rmin*cos(theta);
            t = theta-pi/2+linspace(0,dtheta);
            dx = cenx+rmin*cos(t);
            dy = ceny+rmin*sin(t);
%             x = [x,dx];
%             y = [y,dy];
            theta = theta+dtheta;
            pvec = [dx(end),dy(end),theta];
            dl = dtheta;
        elseif type(i) == 3
            theta = pvec(3);
            dtheta = -seg(i);
            cenx = pvec(1)+rmin*sin(theta);
            ceny = pvec(2)-rmin*cos(theta);
            t = theta+pi/2+linspace(0,dtheta);
            dx = cenx+rmin*cos(t);
            dy = ceny+rmin*sin(t);
%             x = [x,dx];
%             y = [y,dy];
            theta = theta+dtheta;
            pvec = [dx(end),dy(end),theta];
            dl = -dtheta;
        else
            % do nothing
        end
        if dl > 0
            plot(dx,dy,'b');
        else
            plot(dx,dy,'r');
        end
        hold on
    end
    hold off
    axis equal
end