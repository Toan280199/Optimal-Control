function rad_out = ToTrigonometricCircle(rad_in)

% input: any angle
% output: [-pi:pi]
% example: 3.4 pi => -0.6 pi

rad_out = rad_in;

while(rad_out > 1.01*pi )
    rad_out = rad_out - 2*pi;
end
while(rad_out < -1.01*pi)
	rad_out = rad_out + 2*pi;
end