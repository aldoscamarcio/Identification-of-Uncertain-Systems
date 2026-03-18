function MobileRobot
width=1.2;
lth=1;
high=0.3;
radio=0.25;

global Mobile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cuerpo %%%%%%%%%%%%%%%%%%%%%%%%%

base=[19.35 -12.39 0,
18.366 -14.764 0,
16.802 -16.802 0,
14.764 -18.366 0,
12.39 -19.35 0,
9.843 -19.685 0,
9.843 -14.764 0,
-9.843 -14.764 0,
-9.843 -19.685 0,
-12.39 -19.35 0,
-14.764 -18.366 0,
-16.802 -16.802 0,
-18.366 -14.764 0,
-19.35 -12.39 0,
-19.685 -9.843 0,
-19.685 9.843 0,
-19.35 12.39 0,
-18.366 14.764 0,
-16.802 16.802 0,
-14.764 18.366 0,
-12.39 19.35 0,
-9.843 19.685 0,
-9.843 14.764 0,
9.843 14.764 0,
9.843 19.685 0,
12.39 19.35 0,
14.764 18.366 0,
16.802 16.802 0,
18.366 14.764 0,
19.35 12.39 0,
19.685 9.843 0,
18.973 8.297 0,
18.384 6.699 0,
17.922 5.061 0,
17.59 3.391 0,
17.39 1.701 0,
17.323 -0 0,
17.39 -1.701 0,
17.59 -3.391 0,
17.922 -5.061 0,
18.384 -6.699 0,
18.973 -8.297 0,
19.685 -9.843 0]/39.37;

base(:,1)=width*base(:,1);
base(:,2)=lth*base(:,2);

zs=high;

for k=1:length(base)-1
    Mobile.Base{k}=[base(k,1) base(k,2) base(k,3) 1;base(k+1,1) base(k+1,2) base(k,3) 1;base(k+1,1)  base(k+1,2) zs 1 ;base(k,1) base(k,2) zs  1]';
if k==length(base)-1
    Mobile.Base{k+1}=[base(length(base),1) base(length(base),2) base(1,3) 1; base(1,1) base(1,2) base(1,3) 1;base(1,1)  base(1,2) zs 1;base(length(base),1) base(length(base),2) zs 1]';
end
end

Mobile.Base{k+2}=[base'; ones(1,length(base))];
Mobile.Base{k+3}=[[base(:,1) base(:,2) base(:,3)+zs]'; ones(1,length(base))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% LLANTAS %%%%%%%%%%%%%%%%%%%%%%%%%
th=0:0.3:2*pi;
x=radio*cos(th);
z=radio*sin(th);
y=zeros(1,length(th));

wheel=[x' y' z'+high/2];

zs=0.1;

for k=1:length(wheel)-1
    Mobile.Wheel{k}=[wheel(k,1) wheel(k,2) wheel(k,3) 1;wheel(k+1,1) wheel(k,2) wheel(k+1,3) 1;wheel(k+1,1) zs wheel(k+1,3) 1 ;wheel(k,1)  zs wheel(k,3)  1]';
if k==length(wheel)-1
    Mobile.Wheel{k+1}=[wheel(length(wheel),1) wheel(length(wheel),2) wheel(length(wheel),3) 1; wheel(1,1) wheel(length(wheel),2) wheel(1,3) 1;wheel(1,1)  zs wheel(1,3)  1;wheel(length(wheel),1) zs wheel(length(wheel),3)  1]';
end
end

Mobile.Wheel{k+2}=[wheel'; ones(1,length(wheel))];
Mobile.Wheel{k+3}=[[wheel(:,1) wheel(:,2)+zs wheel(:,3)]'; ones(1,length(wheel))];






