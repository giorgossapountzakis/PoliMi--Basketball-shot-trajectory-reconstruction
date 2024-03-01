%% Load video
videoObject = VideoReader('videoplayback.mp4');
numberOfFrames = videoObject.NumberOfFrame;
%%
%This is where we save ball positions
centerFinal=[];
radiiFinal=[];
frame=[];
%Run video for 30 frames, starting from when ball leaves hands
for k = 23:1:53	

% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

thisFrame=read(videoObject,k);
imshow(thisFrame);
axis on;
hold on;
drawnow;
[centers, radii, metric] = imfindcircles(thisFrame,[20 50],'ObjectPolarity','dark','Sensitivity',0.92);
if(length(radii)==1)
    if(centers(2)<500)
    radiiFinal=[radiiFinal ; radii];
    centerFinal=[centerFinal ; centers];
    frame=[frame ; k];
    end
end
if(length(radii)==2)
    centers
    if(centers(1,2)<500)
    radiiFinal=[radiiFinal ; radii(1)];
    centerFinal=[centerFinal ; centers(1,:)];
    frame=[frame ; k];
    else if(centers(2,2)<500)
    radiiFinal=[radiiFinal ; radii(2)];
    centerFinal=[centerFinal ; centers(2,:)];
    frame=[frame ; k];
    end
end

end
end
viscircles(centerFinal, radiiFinal,'EdgeColor','b');
pause(2)
close all;

%%
Img=read(videoObject,30);
%%
%Coordinates of chosen points on the basketball court
ptsEuc=[660 120 300;
    660 120 405;
    840 120 300;
    840 120 405;
    570 579 0;
    930 579 0];
%bottom right backboard,top right backboard,bottom left backboard, top left backboard, far away ft line, close ft
%line,
x=ptsEuc(:,1);
y=ptsEuc(:,2);
z=ptsEuc(:,3);


%%
%imshow(Img)   This is for selecting points manually which we have done
%[xi,yi] = getpts
%axis on
xi=[480.771493212670 ; 490.545248868778; 356.970588235294 ; 363.486425339367; 1376.699095022624; 1625.929864253394];
yi=[345.024886877828 ; 136.518099547511 ; 279.866515837104 ; 17.604072398190 ; 890.726244343891 ; 1068.282805429864];
%% AC=B 6 points
%This is where we compute the camera projection matrix
A=[x(1) y(1) z(1) 1 0 0 0 0 -xi(1)*x(1) -xi(1)*y(1) -xi(1)*z(1);
    0 0 0 0 x(1) y(1) z(1) 1 -yi(1)*x(1) -yi(1)*y(1) -yi(1)*z(1);
    x(2) y(2) z(2) 1 0 0 0 0 -xi(2)*x(2) -xi(2)*y(2) -xi(2)*z(2);
    0 0 0 0 x(2) y(2) z(2) 1 -yi(2)*x(2) -yi(2)*y(2) -yi(2)*z(2);
    x(3) y(3) z(3) 1 0 0 0 0 -xi(3)*x(3) -xi(3)*y(3) -xi(3)*z(3);
    0 0 0 0 x(3) y(3) z(3) 1 -yi(3)*x(3) -yi(3)*y(3) -yi(3)*z(3);
    x(4) y(4) z(4) 1 0 0 0 0 -xi(4)*x(4) -xi(4)*y(4) -xi(4)*z(4);
    0 0 0 0 x(4) y(4) z(4) 1 -yi(4)*x(4) -yi(4)*y(4) -yi(4)*z(4);
    x(5) y(5) z(5) 1 0 0 0 0 -xi(5)*x(5) -xi(5)*y(5) -xi(5)*z(5);
    0 0 0 0 x(5) y(5) z(5) 1 -yi(5)*x(5) -yi(5)*y(5) -yi(5)*z(5);
    x(6) y(6) z(6) 1 0 0 0 0 -xi(6)*x(6) -xi(6)*y(6) -xi(6)*z(6);
    0 0 0 0 x(6) y(6) z(6) 1 -yi(6)*x(6) -yi(6)*y(6) -yi(6)*z(6)];

B=[xi(1) yi(1) xi(2) yi(2) xi(3) yi(3) xi(4) yi(4) xi(5) yi(5) xi(6) yi(6)]';

%C=inv(A'*A)*A'*B;
%C=(A'*A)\(A'*B)
C=A\B;
P=[C(1) C(2) C(3) C(4); C(5) C(6) C(7) C(8); C(9) C(10) C(11) 1];

%Enforce this so that C31x+C32y+C33z + 1 = 1
C(9)=0;
C(10)=0;
C(11)=0;

%% Reconstruction DE=F
fps=videoObject.FrameRate;
t=(frame-23)./fps; %Select one timepoint per frame.
xBall=centerFinal(:,1);
yBall=centerFinal(:,2);

%Here we solve the system to get parameters for the parabolic motion
%parameters
D=[C(1)-xBall(1).*C(9) C(1).*t(1)-xBall(1).*C(9).*t(1)   C(2)-xBall(1).*C(10)   C(2).*t(1)-xBall(1).*C(10).*t(1)   C(3)-xBall(1).*C(11)  C(3).*t(1)-xBall(1).*C(11).*t(1);
 C(5)-yBall(1).*C(9)   C(5).*t(1)-yBall(1).*C(9).*t(1)   C(6)-yBall(1).*C(10)   C(6).*t(1)-yBall(1).*C(10).*t(1)   C(7)-yBall(1).*C(11)  C(7).*t(1)-yBall(1).*C(11).*t(1);
 C(1)-xBall(2).*C(9)   C(1).*t(2)-xBall(1).*C(9).*t(2)   C(2)-xBall(2).*C(10)   C(2).*t(2)-xBall(2).*C(10).*t(2)   C(3)-xBall(2).*C(11)  C(3).*t(2)-xBall(1).*C(11).*t(2);
 C(5)-yBall(2).*C(9)   C(5).*t(2)-yBall(1).*C(9).*t(2)   C(6)-yBall(2).*C(10)   C(6).*t(2)-yBall(2).*C(10).*t(2)   C(7)-yBall(2).*C(11)  C(7).*t(2)-yBall(1).*C(11).*t(2);
 C(1)-xBall(3).*C(9)   C(1).*t(3)-xBall(1).*C(9).*t(3)   C(2)-xBall(3).*C(10)   C(2).*t(3)-xBall(3).*C(10).*t(3)   C(3)-xBall(3).*C(11)  C(3).*t(3)-xBall(1).*C(11).*t(3);
 C(5)-yBall(3).*C(9)   C(5).*t(3)-yBall(1).*C(9).*t(3)   C(6)-yBall(3).*C(10)   C(6).*t(3)-yBall(3).*C(10).*t(3)   C(7)-yBall(3).*C(11)  C(7).*t(3)-yBall(1).*C(11).*t(3);
 C(1)-xBall(4).*C(9)   C(1).*t(4)-xBall(4).*C(9).*t(4)   C(2)-xBall(4).*C(10)   C(2).*t(4)-xBall(4).*C(10).*t(4)   C(3)-xBall(4).*C(11)  C(3).*t(4)-xBall(4).*C(11).*t(4);
 C(5)-yBall(4).*C(9)   C(5).*t(4)-yBall(4).*C(9).*t(4)   C(6)-yBall(4).*C(10)   C(6).*t(4)-yBall(4).*C(10).*t(4)   C(7)-yBall(4).*C(11)  C(7).*t(4)-yBall(4).*C(11).*t(4);
 C(1)-xBall(7).*C(9)   C(1).*t(7)-xBall(7).*C(9).*t(7)   C(2)-xBall(7).*C(10)   C(2).*t(7)-xBall(7).*C(10).*t(7)   C(3)-xBall(7).*C(11)  C(3).*t(7)-xBall(7).*C(11).*t(7);
 C(5)-yBall(7).*C(9)   C(5).*t(7)-yBall(7).*C(9).*t(7)   C(6)-yBall(7).*C(10)   C(6).*t(7)-yBall(7).*C(10).*t(7)   C(7)-yBall(7).*C(11)  C(7).*t(7)-yBall(7).*C(11).*t(7);
 C(1)-xBall(10).*C(9) C(1).*t(10)-xBall(10).*C(9).*t(10) C(2)-xBall(10).*C(10) C(2).*t(10)-xBall(10).*C(10).*t(10) C(3)-xBall(10).*C(11) C(3).*t(10)-xBall(10).*C(11).*t(10);
 C(5)-yBall(10).*C(9) C(5).*t(10)-yBall(10).*C(9).*t(10) C(6)-yBall(10).*C(10) C(6).*t(10)-yBall(10).*C(10).*t(10) C(7)-yBall(10).*C(11) C(7).*t(10)-yBall(10).*C(11).*t(10);
 C(1)-xBall(13).*C(9) C(1).*t(13)-xBall(13).*C(9).*t(13) C(2)-xBall(13).*C(10) C(2).*t(13)-xBall(13).*C(10).*t(13) C(3)-xBall(13).*C(11) C(3).*t(13)-xBall(13).*C(11).*t(13);
 C(5)-yBall(13).*C(9) C(5).*t(13)-yBall(13).*C(9).*t(13) C(6)-yBall(13).*C(10) C(6).*t(13)-yBall(13).*C(10).*t(13) C(7)-yBall(13).*C(11) C(7).*t(13)-yBall(13).*C(11).*t(13);
 C(1)-xBall(16).*C(9) C(1).*t(16)-xBall(16).*C(9).*t(16) C(2)-xBall(16).*C(10) C(2).*t(16)-xBall(16).*C(10).*t(16) C(3)-xBall(16).*C(11) C(3).*t(16)-xBall(16).*C(11).*t(16);
 C(5)-yBall(16).*C(9) C(5).*t(16)-yBall(16).*C(9).*t(16) C(6)-yBall(16).*C(10) C(6).*t(16)-yBall(16).*C(10).*t(16) C(7)-yBall(16).*C(11) C(7).*t(16)-yBall(16).*C(11).*t(16)];
g=-980.665;%gravity
F=[xBall(1).*(0.5.*C(11).*g.*t(1).^2+1)-(0.5.*C(3).*g.*t(1).^2+C(4));
    yBall(1).*(0.5.*C(11).*g.*t(1).^2+1)-(0.5.*C(3).*g.*t(1).^2+C(8));
    xBall(2).*(0.5.*C(11).*g.*t(2).^2+1)-(0.5.*C(3).*g.*t(2).^2+C(4));
    yBall(2).*(0.5.*C(11).*g.*t(2).^2+1)-(0.5.*C(3).*g.*t(2).^2+C(8));
    xBall(3).*(0.5.*C(11).*g.*t(3).^2+1)-(0.5.*C(3).*g.*t(3).^2+C(4));
    yBall(3).*(0.5.*C(11).*g.*t(3).^2+1)-(0.5.*C(3).*g.*t(3).^2+C(8));
    xBall(4).*(0.5.*C(11).*g.*t(4).^2+1)-(0.5.*C(3).*g.*t(4).^2+C(4));
    yBall(4).*(0.5.*C(11).*g.*t(4).^2+1)-(0.5.*C(3).*g.*t(4).^2+C(8));
    xBall(7).*(0.5.*C(11).*g.*t(7).^2+1)-(0.5.*C(3).*g.*t(7).^2+C(4));
    yBall(7).*(0.5.*C(11).*g.*t(7).^2+1)-(0.5.*C(3).*g.*t(7).^2+C(8));
    xBall(10).*(0.5.*C(11).*g.*t(10).^2+1)-(0.5.*C(3).*g.*t(10).^2+C(4));
    yBall(10).*(0.5.*C(11).*g.*t(10).^2+1)-(0.5.*C(3).*g.*t(10).^2+C(8));
    xBall(13).*(0.5.*C(11).*g.*t(13).^2+1)-(0.5.*C(3).*g.*t(13).^2+C(4));
    yBall(13).*(0.5.*C(11).*g.*t(13).^2+1)-(0.5.*C(3).*g.*t(13).^2+C(8));
    xBall(16).*(0.5.*C(11).*g.*t(16).^2+1)-(0.5.*C(3).*g.*t(16).^2+C(4));
    yBall(16).*(0.5.*C(11).*g.*t(16).^2+1)-(0.5.*C(3).*g.*t(16).^2+C(8))];
E=D\F;


%% Plot Basketball court
c1 = [0,0,0];
c2 = [1500,0,0];
c3 = [1500,2800,0];
c4 = [0,2800,0];
plot3( [0 1500 1500 0 0], [0 0 2800 2800 0], [0 0 0 0 0] )
hold on
plot3( [570 930 930 570 570], [0 0 579 579 0], [0 0 0 0 0] )
axis equal
plot3( [570 930 930 570 570], [2800 2800 2800-579 2800-579 2800], [0 0 0 0 0] )
th = linspace( 0, pi, 100);
R = 360/2;  
xC = R*cos(th) + 750;
yC = R*sin(th) + 579;
zC=th.*0;
plot3(xC,yC,zC);
th = linspace( 0, -pi, 100);
xC = R*cos(th) + 750;
yC = R*sin(th) + 2800-579;
zC=th.*0;
plot3(xC,yC,zC);
th = linspace( -pi, pi, 100);
xC = R*cos(th) + 750;
yC = R*sin(th) + 1400;
zC=th.*0;
plot3(xC,yC,zC);

plot3( [0 1500], [1400 1400], [0 0] )
plot3( [660 840 840 660 660], [120 120 120 120 120], [300 300 405 405 300] )
plot3( [660 840 840 660 660], [2680 2680 2680 2680 2680], [300 300 405 405 300])

%% Plot ball trajectory
X=E(1)+E(2)*t;
Y=E(3)+E(4)*t;
Z=E(5)+E(6)*t+0.5*g*t.^2;
plot3(X,Y,Z,'linewidth',4)
