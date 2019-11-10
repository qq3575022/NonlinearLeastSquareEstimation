function [xy_meas, A] = rtoxy(r_sim1, r_sim2, r_sim3, r1dot_e, r2dot_e, r3dot_e, AXY, T)

    x0 = [1.95, 1.51];
    
    x1 = [-0.05, 1.5];
    x2 = [2, 3.0];
    x3 = [2.7, 0.05];

    xy_meas = NaN(6,length(r_sim1));
    A = NaN(2,length(r_sim1));
    
    xy_meas(1,1) = 0.5*((x3(2)-x2(2))*(r_sim1(1)^2-r_sim2(1)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(2)-x1(2))*(r_sim2(1)^2-r_sim3(1)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
    xy_meas(1,1) = xy_meas(1,1)/((x3(2)-x2(2))*(x2(1)-x1(1))-(x2(2)-x1(2))*(x3(1)-x2(1)));
    xy_meas(2,1) = 0;xy_meas(3,1) = 0;
    
    xy_meas(4,1) = 0.5*((x3(1)-x2(1))*(r_sim1(1)^2-r_sim2(1)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(1)-x1(1))*(r_sim2(1)^2-r_sim3(1)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
    xy_meas(4,1) = xy_meas(4,1)/((x3(1)-x2(1))*(x2(2)-x1(2))-(x2(1)-x1(1))*(x3(2)-x2(2)));
    xy_meas(5,1) = 0;xy_meas(6,1) = 0;
    
    A(1,1) = 0;A(2,1) = 0;
    
    for j = 2:1:length(r_sim1)
        xy_meas(1,j) = 0.5*((x3(2)-x2(2))*(r_sim1(j)^2-r_sim2(j)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(2)-x1(2))*(r_sim2(j)^2-r_sim3(j)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
        xy_meas(1,j) = xy_meas(1,j)/((x3(2)-x2(2))*(x2(1)-x1(1))-(x2(2)-x1(2))*(x3(1)-x2(1)));
        xy_meas(2,j) = 0.5*((x3(2)-x2(2))*(r1dot_e(j)^2-r2dot_e(j)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(2)-x1(2))*(r2dot_e(j)^2-r3dot_e(j)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));%(xy_meas(1,j) - xy_meas(1,j-1))/T;
        xy_meas(2,j) = xy_meas(2,j)/((x3(2)-x2(2))*(x2(1)-x1(1))-(x2(2)-x1(2))*(x3(1)-x2(1)));
        xy_meas(3,j) = (xy_meas(2,j) - xy_meas(2,j-1))/T;
        
        xy_meas(4,j) = 0.5*((x3(1)-x2(1))*(r_sim1(j)^2-r_sim2(j)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(1)-x1(1))*(r_sim2(j)^2-r_sim3(j)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
        xy_meas(4,j) = xy_meas(4,j)/((x3(1)-x2(1))*(x2(2)-x1(2))-(x2(1)-x1(1))*(x3(2)-x2(2)));
        xy_meas(5,j) = 0.5*((x3(1)-x2(1))*(r1dot_e(j)^2-r2dot_e(j)^2-x1(1)^2+x2(1)^2-x1(2)^2+x2(2)^2)-(x2(1)-x1(1))*(r2dot_e(j)^2-r_sim3(j)^2-x2(1)^2+x3(1)^2-x2(2)^2+x3(2)^2));
        xy_meas(5,j) = xy_meas(5,j)/((x3(1)-x2(1))*(x2(2)-x1(2))-(x2(1)-x1(1))*(x3(2)-x2(2)));
        
        %xy_meas(5,j) = %(xy_meas(4,j) - xy_meas(4,j-1))/T;
        xy_meas(6,j) = (xy_meas(5,j) - xy_meas(5,j-1))/T;
        
        sin = (xy_meas(1,j) - x0(1))/sqrt((xy_meas(1,j) - x0(1))^2+(xy_meas(4,j) - x0(2))^2);
        cos = (xy_meas(4,j) - x0(2))/sqrt((xy_meas(1,j) - x0(1))^2+(xy_meas(4,j) - x0(2))^2);
        
        A(1,j) =  AXY(2,j)*cos + AXY(1,j)*sin;
        A(2,j) = -AXY(2,j)*sin + AXY(1,j)*cos;
    end

end