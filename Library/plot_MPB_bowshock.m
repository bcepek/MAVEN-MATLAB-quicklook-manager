function [plot1,plot2,plot_3,plot4] = plot_MPB_bowshock
%***************  Bow shock ************
yb=0:0.01:4.5;
xb=(-exp(yb).^0.2)*0.2+0.04*yb-0.2.*yb.^2.15+1.6;
plot1 = plot(xb,yb,'k',xb,-yb,'k'); %******** plot bow shock ******
hold on

%************   Magnetosphere boundary ******
ymb=0:0.01:2.5;
xmb=(-exp(ymb).^2)*0.05+0.1*ymb-0.2.*ymb.^2+1.25;
rmb=sqrt(xmb.^2 + ymb.^2);
plot2 = plot(xmb,ymb,'k',xmb,-ymb,'k'); %***** plot magnetosphere boundary ***
hold on
axis equal tight

t=(3*pi/2:0.001:5/2*pi);

%*************** plot Mars ****
plot_3 = plot3(cos(t),sin(t),zeros(size(t)),'k');  %*************** plot Mars
hold on 
t1=(pi/2:0.001:3*pi/2);
plot4 = area(cos(t1),sin(t1));
end