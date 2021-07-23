%{ 
FEPIC version 01

Taha Y. Posos*, Oksana Chubenko+, and Sergey V. Baryshev* 
*Materials Accelerators and Microscopy Labarotory (MAMLab) - Michigan State University
+Department of Physics - Arizona State University

FEPic_v01 is based on the method given in https://arxiv.org/abs/2012.03578
Equations and terminology referred through the code are in the paper.
The software is developed to compute apparent emission area from given emission pattern.
It assumes discrete emission spots to resolve area. Watch the tutorial to learn how to use this code.

This code can work with any size of images. However, we strongly recommend
using 450x450 pixel image because default parameters are set by assuming
it and because of speed considerations. If you want to use higher resolution
images you need to experiment with many paramenters. To reseize image to 
450x450 pixel simply uncomment line 83.

By default, the code assumes RGB images. If you would like to use grayscale
images, simply comment line 78.

Run the code. It will open a file selection window. Choose a picture file.
Then, the software will compute emission area and visualize the results.
It supports many common file formats (png,jpg,bmp,tif...)
%}

%%%%%%%%%%%%%%%%%%%%%%% Clear %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; % Clear all variables
clc; % Clear command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Fixed Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
99.9% of the time you do not need to change below parameters if images
are 450x450 pixels.
%}
seed = 0;
rng(seed); % fix random seed
r_c = 10; % radius of the search region
A = 300; % amplitude for the decision boundary (Eq.5)
d_mu = 0; % mean for the decision boundary (Eq.5)
radfit = 7; % half-size of the fit region for the Gaussian surface fit
sc = 1; % sigma multiplier for the Gaussian surface fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Tuanable Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
You may need to change these parameters. If the spots are very close to each
other, decrease d_c. If there is glowing background,increase d_s. If the spots
are faint, decrease A3d. Most of the time you do not need to change std3d,
but if spots are small and there is cloudy-like backgound, you may consider
decreasing it to prevent false positives.
%}
d_s = 10; % standard deviation for the decision boundary (Eq.5). Default is 10.
d_c = 5; % offset for the decision boundary (Eq.5). Default is 5.
A3d = 10; % lower limit of amplitude for the Gaussian surface fit (Eq.6). Default is 10.
std3d = 7; % upper limit of standard deviation for the Gaussian surface fit (Eq.6). Default is 7.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% File Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Opens a window for file selection.
[file,path] = uigetfile('*.png');
if isequal(file,0)
    disp('User selected Cancel');
    return
else
    disp(['User selected ', fullfile(path,file)]);
end
data=imread(fullfile(path,file));
[~,file,~] = fileparts(file);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Data Conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Coverts RGB to grayscale. If the input image is already grayscale,
uncomment the following line.
%}
data=rgb2gray(data);
%{
If the picture is not 450x450 pixel, uncomment the following line to reseize image to recomended 
size. This is optional.
%}
% data=imresize(data,[450 450]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Noise Addition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Add a small and recoverable random noise to prevent equalivance probelem.
This will be removed later on.
%}
sz=size(data);
data_noise=double(data) + 0.1*rand(sz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Compute Distance-Value Pairs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part decides MISR and GM pixels, and compute distances for all pixels.
[X,Y]=ndgrid(1:sz(1),1:sz(2));
L=zeros(sz,'uint8'); % matrix to store labels
D=zeros(sz); % matrix to store distances
mymap=data_noise ~= max(data_noise,[],'all'); % decide GM

L(~mymap)=2; % label GM

for index=[X(mymap) Y(mymap)]'
    d=find_distance(data_noise,X,Y,index,sz,r_c); % compute distances for all pixels
    if isempty(d) || (d>r_c)
        L(index(1),index(2))=1; % decide MISR
    else
        D(index(1),index(2))=d; % store distances for ordinary pixels
    end
end

mymap = L==1;
mymap2 = (L==1)|(L==2);

for index=[X(mymap) Y(mymap)]'
    % compute distances for MISR pixels
    D(index(1),index(2))=find_distance2(data_noise(mymap2),X(mymap2),Y(mymap2),data_noise(index(1),index(2)),index);
end

mymap = L==1;
mymap2 = L==2;
D(mymap2)=sqrt(min((X(mymap)-X(mymap2)).^2 + (Y(mymap)-Y(mymap2)).^2)); % compute distance for GM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Apply Decison Boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gauss = @(x) A*exp(-(x-d_mu).^2./(2*d_s^2)) + d_c; % (Eq.5)
i_c = gauss(double(data));
LM = zeros(sz,'uint8'); % matrix to store LM labels
LM(D > i_c) = 1; % apply the decision boundary
LM_coor = [X(LM==1) Y(LM==1)]; % store LM coordinates
LM_count = sum(LM==1,'all'); % count number of possible LMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Emission Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part filters false positives in LM list and computes emission area (unit is px^2)
e_spot = []; % to store coordinates of emission pixels (pixels inside a spot)
discard = []; % to store coordinates of false positives
for index = 1:LM_count
    pxydata = fitdata(data,X,Y,LM_coor(index,:),sz,radfit); % get fit region for LM candidates
    mu = LM_coor(index,:);
    fitform = @(p,x) p(1).*exp(-0.5*((x(:,1)-mu(1))./p(2)).^2 -0.5*((x(:,2)-mu(2))./p(2)).^2) + p(3); % (Eq.6)
    par = lsqcurvefit(fitform,[10 10 10],pxydata(:,2:3),pxydata(:,1),[0 0 0]); % fit a Gaussian surface
    if (par(2)>std3d) || (par(1)<A3d) % find false positives
        discard = [discard; LM_coor(index,:)]; % add them to discard list and not make area computation for them
        continue;
    end
    % get pixels one std away from the center
    mymap = (pxydata(:,2)-mu(1)).^2 + (pxydata(:,3)-mu(2)).^2 < (par(2)*sc).^2;
    e_spot = [e_spot; pxydata(mymap,2:3)]; % add them to emission spot list
end

e_spot = unique(e_spot,'rows'); % if any spots overlap count pixels only once

if ~isempty(discard)
    lindex = sub2ind(sz,discard(:,1),discard(:,2));
    LM(lindex) = 0; % discard false positives from LM list
    LM_count = sum(LM==1,'all'); % count final number of LMs
    LM_coor = [X(LM==1) Y(LM==1)]; % coordinates of filtered LMs
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Visualization - Decision Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
scatter(data(:),D(:),18,'blue','filled');
hold;
scatter(data(LM==1),D(LM==1),18,'red','filled');
plot(0:255,gauss(0:255),'k','LineWidth',1);
hold;
title(sprintf('count = %d',LM_count));
xlim([0 256]);
ylim([0 300]);
xlabel('Pixel Value');
ylabel('Distance');
box on;
% saveas(gcf,strcat(file,'_DP.png')); % uncomment to save the plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Finalize if there is no detected emission spot %%%%%%%%%%%%%%%%
if isempty(e_spot)
    disp('!!!!!!! No emission area is detected !!!!!!!');
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Visualize - Local Maximas (LM) %%%%%%%%%%%%%%%%%%%%%%%%%
RGB = insertMarker(data,[Y(LM==1) X(LM==1)],'color','blue');
figure(2);
imshow(RGB);
% imwrite(RGB,strcat(file,'_LMs.png')); % uncomment to save it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Visualize - Emission Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lindex = sub2ind(sz,e_spot(:,1),e_spot(:,2));
EP_count = size(lindex,1);
R = data;
R(lindex) = 0;
G = data;
G(lindex) = 0;
B = data;
B(lindex) = 255;
RGB = cat(3,R,G,B);
figure(3);
imshow(RGB);
imwrite(RGB,strcat(file,'_earea.png')); % uncomment to save it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Write results and used paramters to an Excel file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file2 = 'Results';
out = {file seed r_c A d_mu radfit sc d_s d_c A3d std3d LM_count EP_count};
c_names = {'file','seed','rc','A','d_mu','radfit','sc','d_s','d_c','A3d','std3d','LM_count','EP_count'};
if isfile(strcat(file2,'.csv'))
    writecell(out,strcat(file2,'.csv'),'WriteMode','append');
else
    writecell(c_names,strcat(file2,'.csv'));
    writecell(out,strcat(file2,'.csv'),'WriteMode','append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Function to get surface fit region for each LM %%%%%%%%%%%
function pxy = fitdata(data,X,Y,index,sz,radfit)
u=index(1)-radfit;
d=index(1)+radfit;
l=index(2)-radfit;
r=index(2)+radfit;
if index(1)-radfit<1
    u=1;
end
if index(1)+radfit>sz(1)
    d=sz(1);
end
if index(2)-radfit<1
    l=1;
end
if index(2)+radfit>sz(2)
    r=sz(2);
end
p = data(u:d,l:r);
x = X(u:d,l:r);
y = Y(u:d,l:r);
p = double(p(:));
x = x(:);
y = y(:);
pxy = [p x y];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Function to calculate distance for ordinary pixels %%%%%%%%%%%%%%
function d=find_distance(data,X,Y,index,sz,r_c)
u=index(1)-r_c;
d=index(1)+r_c;
l=index(2)-r_c;
r=index(2)+r_c;
if index(1)-r_c<1
    u=1;
end
if index(1)+r_c>sz(1)
    d=sz(1);
end
if index(2)-r_c<1
    l=1;
end
if index(2)+r_c>sz(2)
    r=sz(2);
end
mymap=data(u:d,l:r)>data(index(1),index(2));
if any(mymap,'all')
    X_temp=X(u:d,l:r);
    Y_temp=Y(u:d,l:r);
    d=sqrt(min((X_temp(mymap)-index(1)).^2 + (Y_temp(mymap)-index(2)).^2));
else
    d=[];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Function to calculate distance for MISR pixels %%%%%%%%%%%%%%%%%%%%
function d=find_distance2(data,X,Y,value,index)
mymap = data > value;
d=sqrt(min((X(mymap)-index(1)).^2 + (Y(mymap)-index(2)).^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%