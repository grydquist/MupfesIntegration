fileid=fopen('../pos.txt');
cell=fscanf(fileid,'%f');
x1=cell(1:6:end);
y1=cell(2:6:end);
z1=cell(3:6:end);
x2=cell(4:6:end);
y2=cell(5:6:end);
z2=cell(6:6:end);
plot3(cell(1:6:end),cell(2:6:end),cell(3:6:end))
hold on
%plot(cell(3:6:end))
plot3(cell(4:6:end),cell(5:6:end),cell(6:6:end))
pbaspect([1,1,1])
%axis([-1,1,-.2,1.2,10,11.5])

[x,y,z]=sphere(20);
x=x/2;
y=y/2;
z=z/2;

hold off
for i=1:length(cell)/6
  axis([-2,2,-2,2,0,30])
  pbaspect([4,4,30])
  %view(90,0)
  hold on;
  
  surf(x+cell(i*6-5),y+cell(i*6-4),z+cell(i*6-3),'facecolor', 'r', 'edgealpha', 0)
  surf(x+cell(i*6-2),y+cell(i*6-1),z+cell(i*6),'facecolor', 'b', 'edgealpha', 0)
  light;
    lighting gouraud;
  
  plot3(cell(i*6-5),cell(i*6-4),cell(i*6-3),'or');
  plot3(cell(i*6-2),cell(i*6-1),cell(i*6),'or');
  
  pause(0.1);

  clf;
end
% [xc,yc,zc] = cylinder(2);
% zc(2,:) = 30;
% surf(xc,yc,zc)

fclose all;