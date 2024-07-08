close all
clear 

% Plot parameters
height = 0.7
bottom = 0.15
nrplots = 20 
limits = [5100,5700]

datapath = "X:\04_PROJECTS\2022_Lake_response_warming_SNF\Scripts\Finalized_working_scripts\Masterfile_Processed_Soppensee.xlsx"
% sheetnames = ["NGRIP","GERZ","HSIRED","XRFRED","LOWRED","ALPS","Stats"]
sheetnames = ["Carotenoids","HSI","Lowstats"]

for m = 1:length(sheetnames)
    % Data 1 read table Low resolution usually XRF
    data_raw = readcell(datapath,Sheet = sheetnames{m});                       % Loading the dataframe
    datain.headers{m} = data_raw(3:5,:);                                       % Loading the headers of the data (plot titles)
    datain.info{m} = data_raw(6:8,:);                                          % Loading the info of the data (plot type, color, in-or out criterion)
    datain.data{m} = readmatrix(datapath,Sheet = sheetnames{m});               % Loading the data
    datain.data{m} = datain.data{m}(4:size(datain.data{m},1),:);               % Cutting off the front columns containing the ages
    clear data_raw
end

for m = 1:length(datain.data)
    sel.data{m} = datain.data{m}(3:size(datain.data{m},1),cell2mat(datain.info{m}(3,:)) == 1);
    sel.headers{m} = datain.headers{m}(:,cell2mat(datain.info{m}(3,:)) == 1);
    sel.info{m} = datain.info{m}(:,cell2mat(datain.info{m}(3,:)) == 1);
end

ageL = {'11300 +/- 300','11950 +/- 150','13050 +/- 100','13870 +/- 220','14810 +/-300','15780 +/- 800','16530 +/- 1000'}

%% dont run this if cars are sorted

wetn = sel.headers{1}(1,2:end)
% OM = readtable("OM_tab.csv")
wet_ma = sel.data{1}
wet_mat = sel.data{1}(:,2:end)

% OM_mat = readmatrix("OM_tab.csv")
indi = sum(wet_mat == 0) < 20
wet_mat = wet_mat(:,indi)
wetn = wetn(1,indi)
wet_mat = wet_mat(sum(isnan(wet_mat) == 0,2) == max(sum(isnan(wet_mat) == 0,2)),:)
sz = size(wet_mat)

x = 1:sz(2)
y = wet_ma(:,1)
z1 = wet_mat
z = normalize(z1,'zscore')
% z = zscore(z)
r = corr(z);
labels = wetn
% rtab = array2table(r,"VariableNames",labels,"RowNames",labels)
% sortedonce = sorty(rtab)
% sortedtwice = sortx(sortedonce)
% figure, h = heatmap(sortedtwice,'MissingDataColor','w','Colormap',parula);

[xq,yq] = meshgrid(1:0.2:sz(2), 5189.5:1:5704.5);
vq = griddata(x,y,z,xq,yq,"nearest");
figure,
mesh(xq,yq,vq)
xticks(x)
xticklabels(labels)
xlim([1,sz(2)-2])
ylim(limits)
ylabel("mm")

% -- Running on reversed matrix ---

figure,subplot(2,1,1)
dendrogram(linkage(z','ward'),0,'ColorThreshold',7.5);set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)
M = cluster(linkage(z','ward'),'cutoff',7.5,'Criterion','distance');

zflip = z'
zflip(:,size(z,1)+1) = M
temp = sortrows(zflip,size(z,1)+1,"ascend")
temp2 = temp(:,1:size(z,1)+1)'

% temp2(:,size(z,2)+1) = M
% double_sort = sortrows(temp2,size(z,2)+1,"ascend")
% double_sort = double_sort(:,1:size(z,2))
[mat, idx] = sort(M,"ascend")
sorted_labels = labels(idx)

figure,
[xq,yq] = meshgrid(1:0.1:50, 5007:1:5797);
temp2(end,:) = []
vq = griddata(x,y,temp2,xq,yq,"nearest");

labels = sorted_labels
figure,
mesh(xq,yq,vq)
colormap('turbo')
xticks(x)
xticklabels(labels)
xlim([1,size(wet_mat,2)-2])
ylim([5007,5797])
set(gca, 'YDir', 'reverse')
ylabel("mm")

% --- running clusters on corr matrix ---

figure,subplot(2,1,1)
dendrogram(linkage(r,'ward'),0,'ColorThreshold',1.5);set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)
M = cluster(linkage(r,'ward'),'cutoff',1.5,'Criterion','distance');

r(size(z,2)+1,:) = M
temp = sortrows(r',size(z,2)+1,"ascend")
temp2 = temp(:,1:size(z,2))'
temp2(:,size(z,2)+1) = M
double_sort = sortrows(temp2,size(z,2)+1,"ascend")
double_sort = double_sort(:,1:size(z,2))

[mat, idx] = sort(M,"ascend")
sorted_labels = labels(idx)

figure,
% 
% isupper = logical(triu(ones(size(r)),1));
% r(isupper) = NaN;

h = heatmap(double_sort,'MissingDataColor','w','Colormap',parula);
h.XDisplayLabels = sorted_labels;
h.YDisplayLabels = sorted_labels;

title('Correlation matrix');xlabel('Proxies'); ylabel('Proxies');

% ---

z(44,:) = M/2
tempz = sortrows(z',44,"ascend")
z = tempz(:,1:43)'
figure,
[xq,yq] = meshgrid(1:0.1:length(M), 5007:1:5797);
vq = griddata(x,y,z,xq,yq,"nearest");

labels = sorted_labels
figure,
surf(xq,yq,vq,'EdgeColor','none')
colormap('turbo')
c = colorbar()
c.Ticks = [-1, 1, 4, 7]; c.TickLabels =  {'absent','present','abundant','peak'}
xticks(x)
xticklabels(labels)
xlim([1,size(wet_mat,2)-2])
ylim([5007,5797])
set(gca, 'YDir', 'reverse')
ylabel("mm")

figure, heatmap(x,y,z)
colormap('turbo')
xticks(x)
xticklabels(labels)
xlim([1,size(wet_mat,2)-2])
ylim([5007,5797])
set(gca, 'YDir', 'reverse')
ylabel("mm")

array2table(temp2,"VariableNames",labels)

%% continue here

z = sel.data{1}(:,[8:size(sel.data{1},2)]);
labels = sel.headers{1}(1,[8:size(sel.data{1},2)]);
y = sel.data{1}(:,1);
x = 1:1:size(sel.data{1},2)-7;
z = zscore(z);

figure, subplot(1,11,[5:10])

[xq,yq] = meshgrid(0:1:max(x)+1, 5100:10:5700);
vq = griddata(x,y,z,xq,yq,"nearest");
surf(xq,yq,vq,'EdgeColor',[0.1 0.1 0.1],'LineStyle',':','LineWidth',0.3)
colormap('turbo')
c = colorbar();
c.Ticks = [-1, 1, 4, 7]; c.TickLabels =  {'absent','present','abundant','peak'};
xticks(x)
xticklabels(labels)
xlim([1,max(x)])
ylim(limits)
xlim([1 27])
set(gca, 'YDir', 'reverse')
ylabel("mm")
yticks([5000:100:5700])
yticklabels(AgeL)
set(gca,"YAxisLocation","right","YColor","none","FontName","Gill Sans MT", "FontSize",12,"XTickLabelRotation",90,"View",[0.42,90],"CameraPosition",[23,5351,88.58545955248599],"CameraTarget",[23,5351,2.276277865933995],'CameraViewAngle',6.608610360311924)
title('z-scored carotenoids')

% Phosphorous
clear x y z vq xq yq

colorarrayP = ([11 102 35;112 130 56;157 193 131;208 240 192]/255)

x = sel.data{1}(:,1);
subplot(1,11,2),
b = bar(x,sel.data{1}(:,2:5),1,'stacked','EdgeColor',[0 0 0],"Horizontal","on") % Edgecolor is set to black but can be changed
hold on, stairs(sum(sel.data{1}(:,2:5),2),x+5,"LineWidth",0.5,"Color",'k') 
hold on, plot(sum(sel.data{1}(:,2:5),2),x,"LineStyle",'none',"Marker",".","MarkerEdgeColor",colorarrayP(1,:),'MarkerSize',10)
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','left');ylabel('depth (mm)');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:4
    st(h) = string(sel.headers{1}(1,1+h))
    b(h).FaceColor = colorarrayP(h,:);
end
legend(st,"EdgeColor",'none',"Color","none")
xlabel("Green pigments nmol/g_{wet}")
ylim(limits); % Depths
xticks([5:10:45])

sp = subplot(1,11,1)
x = sel.data{3}(:,1); y = sel.data{3}(:,2)
xst = x+((x(5)-x(4))/2)
x = xst
nrofswitch = sum(abs(diff(y))>0)
fY = y(abs(diff(y))>0)
fY(2:nrofswitch+1) = y(find(abs(diff(y))>0)+1)
rectXarr(1) = min(x)
rectXarr(2:nrofswitch+1) = x(abs(diff(y))>0)
rectXarr(nrofswitch+2) = max(x)
Clr = bone(max(y))
for z = 1:length(rectXarr)-1
    rectX = [rectXarr(z),rectXarr(z+1)]
    rectY = xlim([sp]);
    pch(z) = patch(sp, rectY([1 1 2 2]), rectX([1,2,2,1]),Clr(fY(z),:),'EdgeColor', 'k', 'FaceAlpha', 0.3); % FaceAlpha controls transparency
    hold on,
end
clear rectXarr fY rectY rectX 
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
box off;
ylim(limits); % Depths
% hold on, scatter(sel.data{3}(:,2),sel.data{3}(:,1),20,sel.data{3}(:,2)) 

subplot(1,11,3)
plot(sel.data{2}(:,3),sel.data{2}(:,1),"LineWidth",2,"LineStyle","-","Color",sel.info{2}{2,2}); hold on,
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right','YColor','none');ylabel('depth (mm)');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");ylim(limits)
xlabel("Total chl-a derivatives \mug/g_{wet}");box off; 
xticks([25:25:175])

% Purple pigments
subplot(1,11,4),
plot(sel.data{2}(:,5),sel.data{2}(:,1),"LineWidth",2,"LineStyle","-","Color",sel.info{2}{2,3}); hold on,
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right','YColor','none');ylabel('depth (mm)');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
xlabel("Bphe-a sum \mug/g_{wet}");ylim(limits); % Depths
xticks([5:10:40])
box off; 

subplot(1,11,11),
x = sel.data{1}(:,1);
bar(sel.data{1}(:,1),sel.data{1}(:,7),1,'EdgeColor','none',"Horizontal","on","FaceColor",[230, 100, 56]/255,"FaceAlpha",0.6); hold on,
hold on, stairs(sel.data{1}(:,7),x+5,"LineWidth",0.5,"Color",'k') 
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel('depth (mm)');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
xlabel("Total carotenoids nmol/g_{wet}");ylim(limits); % Depths
yticks([5000:100:5700])
yticklabels(AgeL)

box off; 
