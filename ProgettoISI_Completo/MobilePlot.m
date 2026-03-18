function  Mobile_Graph=MobilePlot(dx,dy,angz)
global Mobile;
dz=0;

% Matriz de rotación z
Rz=[ cos(angz), -sin(angz) 0 0; sin(angz) cos(angz) 0 0; 0 0 1 0;0 0 0 1];
Rot_TheWholeArm=Rz;
Rot_TheWholeArm(1:3,4) = [dx dy dz]';

WheelColor='g';
BodyColor=[0 0 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUERPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tam=0;
for ii = 1:(length(Mobile.Base))
    robotPatch = Rot_TheWholeArm*Mobile.Base{ii};
    Mobile_Graph(tam+ii) = patch(robotPatch(1,:),robotPatch(2,:),robotPatch(3,:),BodyColor,'LineWidth',1.8);
end

Rot_TheWholeArm(1:3,4) = [dx-0.38*sin(angz) dy+0.38*cos(angz) dz]';

tam=tam+ii;

for ii =1:(length(Mobile.Wheel))
    robotPatch = Rot_TheWholeArm*Mobile.Wheel{ii};
    Mobile_Graph(tam+ii) = patch(robotPatch(1,:),robotPatch(2,:),robotPatch(3,:),WheelColor,'LineWidth',1.8);
end

Rot_TheWholeArm(1:3,4) = [dx+((0.38+0.1)*sin(angz)) dy-((0.38+0.1)*cos(angz)) dz]';

tam=tam+ii;
for ii =1:(length(Mobile.Wheel))
    robotPatch = Rot_TheWholeArm*Mobile.Wheel{ii};
    Mobile_Graph(tam+ii) = patch(robotPatch(1,:),robotPatch(2,:),robotPatch(3,:),WheelColor,'LineWidth',1.8);
end





