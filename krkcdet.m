
while 1>0
fclose all;
fileid=fopen('../pos.txt');
cell=fscanf(fileid,'%f');
x000 = cell(1:8:end);
x100 = cell(2:8:end);
x101 = cell(3:8:end);
x001 = cell(4:8:end);
x010 = cell(5:8:end);
x110 = cell(6:8:end);
x111 = cell(7:8:end);
x011 = cell(8:8:end);



clf;
%plot(x000./x000);
hold on
plot(x100./x000);
plot(x010./x000);
plot(x001./x000);
plot(x110./x000);
plot(x011./x000);
plot(x101./x000);
plot(x111./x000);

pause(0.3)

end
fclose all;