function [h] = drawEllipse(ellipseParams, color, style)
	t = 0:0.01:2*pi;
	ex=ellipseParams(3)+ellipseParams(1)*cos(t)*cos(ellipseParams(5)) - ellipseParams(2)*sin(t)*sin(ellipseParams(5));
	ey=ellipseParams(4)+ellipseParams(2)*sin(t)*cos(ellipseParams(5)) + ellipseParams(1)*cos(t)*sin(ellipseParams(5));
	h = plot(ex,ey,[color style],'LineWidth',2);
end
