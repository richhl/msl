function drawCircle()
	theta = linspace(0,2*pi,500);
	fill(cos(theta), sin(theta),[.7 .7 .7], 'LineWidth', 2);
	axis equal;
end
