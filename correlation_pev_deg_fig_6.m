%%
clear all 
%clc
# Please first run the code: "wheel_random_regular_model_fig_4.py" 
# contact information regarding codes
# "priodyuti pradhan" <priodyutipradhan@gmail.com>
# Complex Systems Lab, Indian Institute of Technology Indore 


path1 = 'wheel_RRd_deg_26408.txt';
path2 = 'wheel_RR_graph_pev_dloc_26408.txt';

formatSpec = '%d %f';
sizeA = [2 Inf];

fd1 = fopen(path1,'rt');
fd2 = fopen(path2,'rt');

Y1 = fscanf(fd1,formatSpec,sizeA);
Y1 = Y1';
fclose(fd1);

Y2 = fscanf(fd2,formatSpec,sizeA);
Y2 = Y2';
fclose(fd2);

d = Y1(:,2);

deg = d./sqrt(d'*d);
sqrt(deg'*deg)
pev = Y2(:,2);
corr(deg,pev)


