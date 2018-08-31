function [ ] = plot_slip(opt_results, record_video)
%plot_slip gives visualization for slip results
%   Feed it an opt_results struct and watch the pretty output
     if record_video==1
        v=VideoWriter('vel_4','MPEG-4');
        v.FrameRate=10;
        open(v);
     end    

    s = opt_results;
    
    Fs = s.param.k * (s.r0 - s.r);
    xcop = -s.Tankle * s.param.transmission_ankle .* s.r ./ (Fs .* s.y);
    leg_angle = atan2(s.y,s.x) - pi/2;
    
    %X Y Animation
    slipbodyanimate=figure('Name','Slip body animation','pos',[1921,0,1921,1000]);
    figure(slipbodyanimate);
    ax1 = subplot(4,4,[1,5,9]);
    ax1.XLim = [-.7 .7];
    ax1.YLim = [-.1 1.2];
    axis equal
    %COM circle
    theta = linspace(0, 2 * pi, 100);
    x_c = 0 + .1 * cos(theta);
    y_c = 0 + .1 * sin(theta);
    h(1) = patch(x_c, y_c, 'red');
    %Spring Line
    x_l = [-.01, .01, .01, -.01];
    y_l = [0, 0, -s.r(1), -s.r(1)];
    h(2) = patch(x_l, y_l, 'blue');
    tf = hgtransform('Parent', ax1);
    set(h, 'Parent', tf)
    rline = refline([0 0]);
    rline.Color = 'r';
    rline.LineStyle = '--';
    title('X Y SLIP')
    %Trajectory line
    hold on;
    plot(s.x, s.y, 'b--')

    %COP Animation
    figure(slipbodyanimate)
    ax2 = subplot(4,4,2);
    ax2.XLim = [-.2 .2];
    ax2.YLim = [-.1 .2];
    x_foot = [-s.param.lf/2, s.param.lf/2, s.param.lf/2, -s.param.lf/2];
    y_foot = [0, 0, .05, .05];
    patch(x_foot, y_foot, 'blue') %foot
    arrowx = .25*[0, .1, 0, -.1];
    arrowy = .25*[0, -.2, -.1, -.2];
    arrow = patch(arrowx, arrowy, 'red');
    arrow_move = @(x_pos) arrowx + x_pos;
    title('COP')

    %COP plot
    figure(slipbodyanimate)
    ax3 = subplot(4,4,6);
    cop_plot = plot(0,0);
    ax3.XLim = [min(s.t), max(s.t)];
    ax3.YLim = [-s.param.lf/2 - .01, s.param.lf/2 + .01];
    cop_ref_upper = refline(ax3, [0, s.param.lf/2]); cop_ref_upper.Color = 'r'; cop_ref_upper.LineStyle = '--';
    cop_ref_lower = refline(ax3, [0, -s.param.lf/2]); cop_ref_lower.Color = 'r'; cop_ref_lower.LineStyle = '--';
    title('COP')

    %Ankle torque plot
    figure(slipbodyanimate)
    ax4 = subplot(4,4,[3,4,7,8]);
    ankle_torque_plot = plot(0,0);
    ax4.XLim = [min(s.t), max(s.t)];
    ax4.YLim = [min(s.Tankle), max(s.Tankle)];
    title('Ankle torque trajectory')

    %Leg animation
    figure(slipbodyanimate)
    ax5 = subplot(4,4,[10,14]);
    leg_anim = plot(0,0);
    ax5.XLim = [-.2, .2];
    ax5.YLim = [0, max(s.r0) + .1];
    %COM circle
    theta = linspace(0, 2 * pi, 100);
    x_c = 0 + .1 * cos(theta);
    y_c = 0 + .1 * sin(theta);
    circle = patch(x_c, y_c, 'r');
    new_circle_y = @(y) .1 * sin(theta) + y;
    %Spring
    x_l = [-.01, .01, .01, -.01];
    y_l = [0, 0, 1.4, 1.4];
    patch(x_l, y_l,'b')
    title('Spring length animation')

    %Leg torque plot
    figure(slipbodyanimate)
    ax6 = subplot(4,4, [11,12,15,16]);
    leg_torque_plot = plot(0,0);
    ax6.XLim = [min(s.t), max(s.t)];
    ax6.YLim = [min(s.Tleg), max(s.Tleg)];
    title('Leg torque trajectory')


    for q = 1:2 %Number of times to play
        for i = 1:length(s.x)
            %Body animation stuff
            %Repatch the spring line
            h(2).Vertices(3:4,2) = -(s.r(i) + 0);
            %HG transform stuff
            Tx = makehgtform('translate',[s.x(i),s.y(i),0],'zrotate',leg_angle(i));
            set(tf, 'Matrix', Tx)


            %COP animation
            %Repatch the arrow                       
            arrow.Vertices(:,1) = arrow_move(xcop(i));

            %COP Plot
            cop_plot.XData = s.t(1:i);
            cop_plot.YData = xcop(1:i);

            %Ankle torque plot
            ankle_torque_plot.XData = s.t(1:i);
            ankle_torque_plot.YData = s.Tankle(1:i);

            %Leg length animation
            circle.Vertices(:,2) = new_circle_y(s.r0(i));

            %Leg torque plot
            leg_torque_plot.XData = s.t(1:i);
            leg_torque_plot.YData = s.Tleg(1:i);

            drawnow
            pause(.2)

                    if record_video==1
                        F=getframe(gcf);
                        writeVideo(v,F);
                    end
        end
    end            
            if record_video == 1
                close(v)
            end

end

