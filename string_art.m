clear all; close all; clc;

% settings
frame_radius = 250; % mm
thread_width = .15; % mm
num_of_nails = 316; %200;
num_of_lines = 4000;
grid_size = 5; % mm, for computing the density
pxl_size = thread_width;

grid_num = round(2*frame_radius/grid_size);

% define canvas
aux = ceil((2*frame_radius)/pxl_size);
cnvs = zeros(aux);

% load image
im = imread('besancon.jpg');
im = rgb2gray(im);

% crop image to be square (take the biggest square region from the center of the image)
side = min(size(im));
im = im(floor((size(im,1)-side)/2)+1:floor((size(im,1)-side)/2)+1+side-1,floor((size(im,2)-side)/2)+1:floor((size(im,2)-side)/2)+1+side-1);
im = imresize(im,size(cnvs));
im = 1-double(im)/255;

% nail coordinates
xc = frame_radius*cos(2*pi/num_of_nails*(1:num_of_nails));
yc = frame_radius*sin(2*pi/num_of_nails*(1:num_of_nails));

% prepare cost function
[~,scaled_im] = cost_fcn(cnvs,im,grid_num);

% prepare canvas
xcnvs = linspace(-frame_radius,frame_radius,size(cnvs,1));
ycnvs = xcnvs;
[X,Y] = meshgrid(1:size(im,2),1:size(im,1));
xc_pxl = interp1(xcnvs,1:size(im,2),xc,'nearest');
yc_pxl = interp1(xcnvs,1:size(im,1),yc,'nearest');
cnvs = zeros(size(im));
lines = floor(rand(num_of_lines,1)*num_of_nails)+1;
objective = [];
best_f = inf;
l = 0;
%%
figure(2); 
set(gcf,'visible','off');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1);
h2 = imagesc(-scaled_im); title('original density');
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); c = caxis;
subplot(1,2,2);
h3 = imagesc(zeros(size(-scaled_im))); t2 = title(['canvas density (0 succesful mutations)']);
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); caxis(c);

% figure(3); 
% % set(gcf,'visible','off');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% h = imagesc(1-zeros(size(1-cnvs))); t = title(['canvas (0 succesful mutations)']); set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray);

figure(4);
subplot(1,2,1);
ho = plot([0 0],[0 0]);
grid;
xlabel('iteration');
ylabel('objective');
subplot(1,2,2);
ho2 = plot([0 0],[0 0]);
grid;
xlabel('iteration');
ylabel('objective');
pause(eps);

%% optimize
recording = NaN*ones(num_of_lines,10000);
while 1 % number of iterations
%     profile on
    % perform random mutation
    node = floor(rand()*numel(lines)+1);
    current_nail = lines(node);
    candidate_nail = floor(rand()*num_of_nails+1);
    next_node = node+1; if next_node>numel(lines), next_node = 1; end
    prev_node = node-1; if prev_node<1, prev_node = numel(lines); end
    if candidate_nail==current_nail || candidate_nail==lines(next_node) || candidate_nail==lines(prev_node), continue; end
    lines(node) = candidate_nail;
        
	% plot into canvas
    cnvs = zeros(size(cnvs));
    for k=1:numel(lines)-1
        [aux_x, aux_y]=bresenham(xc_pxl(lines(k)),yc_pxl(lines(k)),xc_pxl(lines(k+1)),yc_pxl(lines(k+1)));
        ind = sub2ind(size(cnvs), aux_y, aux_x);
        cnvs(ind) = 1;
    end
        
    % compare
    [f,~,scaled_cnvs] = cost_fcn(cnvs);
    if f<best_f
        best_f = f;
        objective = [objective; best_f];
        l = l+1; % count successful mutations
        if l<=10000
            recording(:,l) = lines;
        end
        
        % plot
    %     h.CData = -cnvs;
        h2.CData = -scaled_im;
        h3.CData = -scaled_cnvs;
%         t.String = ['canvas (' num2str(l) ' succesful mutations)'];
        t2.String = ['canvas density (' num2str(l) ' succesful mutations)'];
        ho.XData = 1:numel(objective);
        ho.YData = objective;
        ho2.XData = 1:min(numel(objective),10); % last ten values
        ho2.YData = objective(end-(min(numel(objective),10)-1):end);
        pause(eps);
        datetime
    else
        lines(node) = current_nail;
    end
    
    % get rid of too short segments
    % TBD: optionally get rid of two short segments

%     profile viewer
end
%% clean after finding the image
close(fig); close all;
%% show the result
figure(2); 
% set(gcf,'visible','off');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1);
h2 = imagesc(-scaled_im); title('original density');
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); c = caxis;
subplot(1,2,2);
h3 = imagesc(zeros(size(-scaled_im))); t2 = title(['canvas density (0 succesful mutations)']);
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); caxis(c);
h2.CData = -scaled_im;
h3.CData = -scaled_cnvs;
%         t.String = ['canvas (' num2str(l) ' succesful mutations)'];
t2.String = ['canvas density (' num2str(l) ' succesful mutations)'];
%% display stats
lines = [lines(1:end-1), lines(2:end)];
string_lenght = sum(sqrt((xc(lines(:,1))-xc(lines(:,2))).^2+(yc(lines(:,1))-yc(lines(:,2))).^2));
display(['total string length: ' num2str(string_lenght/1000) ' m']);
%% make video
%% optimization process
file = 'optimization_process.mat';

load(file);
[~,scaled_im] = cost_fcn(cnvs,im,grid_num*2); % set the resolution of images

% setup the figure dimensions and appeareance
fig = figure(777); 
% set(gcf,'visible','off');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(1,2,[0.01],0,[0.02]);
axes(ha(1));
h2 = imagesc(-scaled_im); title('original');
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); c = caxis;
set(gca,'visible','off')
axes(ha(2));
h3 = imagesc(zeros(size(-scaled_im))); t2 = title(['canvas (0 succesful mutations)']);
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); caxis(c);
set(fig,'color','white')
set(findall(fig,'-property','FontSize'),'FontSize',10)
fig.Units = 'inch';
aux = get(fig, 'Position'); 
% aux(3) = 10;
aux(4) = 3.2;
set(fig, 'Position', aux);
set(gca,'visible','off')

% initialize video file
% v = VideoWriter('optimization_process.mp4','MPEG-4');
% v.FrameRate = 1;
v = VideoWriter('optimization_process');
v.FrameRate = 1;
open(v);

num_of_iterations = find(isnan(recording(1,:)),1,'first')-1;
for l = 1:num_of_iterations
    l
    % plot into canvas
    cnvs = zeros(size(cnvs));
    for k=1:numel(recording(:,l))-1
        [aux_x, aux_y]=bresenham(xc_pxl(recording(k,l)),yc_pxl(recording(k,l)),xc_pxl(recording(k+1,l)),yc_pxl(recording(k+1,l)));
        ind = sub2ind(size(cnvs), aux_y, aux_x);
        cnvs(ind) = 1;
    end
    % compare
    [f,~,scaled_cnvs] = cost_fcn(cnvs);
    % update plot
    h3.CData = -scaled_cnvs;
    t2.String = ['canvas density (' num2str(l) ' succesful mutations)'];
    if l==1, pause(10); end
    pause(eps);
    % write to video file
    frame = getframe(gcf);
    writeVideo(v,frame); 
end
close(v);
%% knitting process
file = 'knitting_process.mat';

load(file);
grid_num = 2*grid_num;
[~,scaled_im] = cost_fcn(cnvs,im,grid_num); % set the resolution of images

% setup the figure dimensions and appeareance
fig = figure(777); 
% set(gcf,'visible','off');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(1,2,[0.01],0,[0.02]);
axes(ha(1));
h2 = imagesc(-scaled_im); title('original');
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); c = caxis;
set(gca,'visible','off')
axes(ha(2));
h3 = imagesc(zeros(size(-scaled_im))); t2 = title(['canvas (0 succesful mutations)']);
set(gca,'xtick',[]); set(gca,'ytick',[]); axis equal tight; colormap(gray); caxis(c);
set(fig,'color','white')
set(findall(fig,'-property','FontSize'),'FontSize',10)
fig.Units = 'inch';
aux = get(fig, 'Position'); 
% aux(3) = 10;
aux(4) = 3.2;
set(fig, 'Position', aux);
set(gca,'visible','off')

% initialize video file
v = VideoWriter('knitting_process');
v.FrameRate = 0.5;
open(v);

grid_pxl_size = round(size(im,1)/grid_num);
R = round(size(scaled_im,1)/2);
[X,Y] = meshgrid(1:size(scaled_im,1));
mask = logical(((X-R).^2+(Y-R).^2)>R^2);
cnvs = zeros(size(cnvs));
for k=1:numel(lines)-1
    [aux_x, aux_y]=bresenham(xc_pxl(lines(k)),yc_pxl(lines(k)),xc_pxl(lines(k+1)),yc_pxl(lines(k+1)));
    ind = sub2ind(size(cnvs), aux_y, aux_x);
    cnvs(ind) = 1;
end
scaled_cnvs = blkproc(cnvs, [grid_pxl_size grid_pxl_size], 'mean2');
scaled_cnvs(mask) = 0;
final_scaling = sum(scaled_cnvs(:));
for l = 1:numel(lines)-1
    % plot into canvas
    cnvs = zeros(size(cnvs));
    for k=1:l
        [aux_x, aux_y]=bresenham(xc_pxl(lines(k)),yc_pxl(lines(k)),xc_pxl(lines(k+1)),yc_pxl(lines(k+1)));
        ind = sub2ind(size(cnvs), aux_y, aux_x);
        cnvs(ind) = 1;
    end
    % compare
%     [f,~,scaled_cnvs] = cost_fcn(cnvs);
    scaled_cnvs = blkproc(cnvs, [grid_pxl_size grid_pxl_size], 'mean2');
    scaled_cnvs(mask) = 0;
    scaled_cnvs = scaled_cnvs./final_scaling;

    % update plot
    h3.CData = -scaled_cnvs;
    t2.String = ['canvas density (' num2str(l) ' succesful mutations)'];
    pause(eps);
    % write to video file
    frame = getframe(fig);
    writeVideo(v,frame);
end
close(v);
%% generate auxiliary svg for number labels around the circular frame
% generating SVG
fid = fopen('numbers2.svg','w'); % save plot as SVG
fprintf(fid,'<?xml version=\"1.0\" standalone=\"no\"?><!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\"><svg width=\"10mm\" height=\"10mm\" viewBox="0 0 100 100" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"><desc>numbers</desc>');
for k=1:num_of_nails 
    fprintf(fid,'<g>');
    fprintf(fid,['<text x=\"0\" y=\"' num2str(k*1e-3) '\" text-anchor=\"middle\" font-size=\"12\">' num2str(k) '</text>']); % face number
    fprintf(fid,'</g>');
end
fprintf(fid,'</svg>');
fclose(fid);
%% plot the "knitting" guide
global k play tau num_of_lines
load('knitting_process.mat');
% lines = [1:316].';
lines = [lines(1:end-1), lines(2:end)];

offset = 316/2-7;

figure(1);
plot(xc,yc,'-','Color',.2*ones(3,1));
axis equal tight; set(gca,'xtick',[]); set(gca,'ytick',[]);
hold on;
hl1 = plot([0 0],[0 0],'r.-','MarkerSize',20,'LineWidth',5);
% hl2 = plot([0 0],[0 0],'.-','MarkerSize',20,'Color',[255 204 221]./255);
tn = text(0,0,num2str(0),'FontSize',150,'HorizontalAlignment','center','VerticalAlignment','middle'); % print the goal nail
hold off;
% set(gcf,'units','centimeters','outerposition',[0 0 2*frame_radius/10 2*frame_radius/10]);
set(gcf,'KeyPressFcn',@myfun);

k=1961;
play = 1;
tau = 4;
num_of_lines = size(lines,1);
disp(['start at ' num2str(lines(1,1))]);
%%
k_last = NaN;
while 1
    if k~=k_last
        k_last = k;
        hl1.XData = [xc(mod(lines(k,1)+offset,num_of_nails)+1) xc(mod(lines(k,2)+offset,num_of_nails)+1)];% xc(mod(lines(k+1,2)+offset,num_of_nails)+1)];
        hl1.YData = [yc(mod(lines(k,1)+offset,num_of_nails)+1) yc(mod(lines(k,2)+offset,num_of_nails)+1)];% yc(mod(lines(k+1,2)+offset,num_of_nails)+1)];   
    %     hl2.XData = xc(mod(lines(max(1,k-5):k,1)+offset,num_of_nails)+1);
    %     hl2.YData = yc(mod(lines(max(1,k-5):k,1)+offset,num_of_nails)+1);
        tn.String = num2str(lines(k,2));
        pause(eps);
        tts(num2str(lines(k,2)),'Microsoft Hortense Desktop - French'); % say the goal nail
    end
    if play
        pause(tau);
        if play % checking play for the second time since it could have changed
            k=min(k+1,num_of_lines-1);
        end
    end
    pause(eps);
end