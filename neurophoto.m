%% neurophoto.m 
% 
% Description
% ---------------------------------------------------------------
%
%     neurophoto provides a standard (but customizable) 
%     analysis pipeline for neural image sets. Specifically
%     this app is optimized for timeseries tiff stacks 
%     from 2-photon microscopy experiments using GCaMP
%     fluorescent indicators or other photoactivatable
%     bio molecules.
%
%
%     This file has a companion script called...
%
%         neurophoto_tutorial.m
%
%     ...which provides an annotated step-by-step walkthrough
%     of the code, and is intended to help you make modifications
%     to this standard pipeline based on your experiment particulars
%     and analysis needs.
%     
% ---------------------------------------------------------------


%#######################################################################
%% NEUROPHOTO STARTUP - CLEAR WORKSPACE & ADD NEUROPHOTO DIR TO PATH
%#######################################################################


clc; close all; clear;
P.home = pwd; P.home = fileparts(which('neurophoto.m')); cd(P.home)
if ~any(regexp(P.home,'neurophoto') > 0)
    disp(['Run this code from the neurophoto directory; '...
      'your current working directory is currently:'])
    disp(P.home);
end
P.funs  = [P.home filesep 'neurophoto_functions'];
P.imgs  = [P.home filesep 'neurophoto_images'];
P.data  = [P.home filesep 'neurophoto_outputdata'];
addpath(join(string(struct2cell(P)),':',1))





%###############################################################
%% SELECT A TIFF STACK AND GET INFO
%###############################################################

[PIX] = getIMGpaths();

clearvars -except PIX





%###############################################################
%% IMPORT TIFF STACK
%###############################################################


IMG = IMPORTimages(PIX);

clearvars -except PIX IMG



%###############################################################
%% TRIM FIRST FEW FRAMES OF IMG STACK
%###############################################################


IMG(:,:,1:10) = [];

clearvars -except PIX IMG




%###############################################################
%% PREPROCESS TIFF STACK
%###############################################################



[IMG, BND] = PREPROCESSimages(IMG);


% IMG = uint8(rescale(IMG).*255);

IMG = rescale(IMG);

viewstack(IMG,.05)


clearvars -except PIX IMG BND




%###############################################################
%%               ADJUST IMAGE CONTRAST
%###############################################################



% IMG = adjustContrast(IMG);
% 
% clearvars -except PIX IMG
% 
% imstats(IM_A)






%###############################################################
%%              SMOOTH IMAGE
%###############################################################


% IM = smoothIMG(IMG);

SMIM = imgaussfilt3(IMG, 2);

close all
imagesc(SMIM(:,:,1))
colormap hot
title('FIRST FRAME OF GAUSSIAN SMOOTHED IMAGE STACK')
pause(2)



clearvars -except PIX IMG SMIM





%###############################################################
%%                        RUN PCA
%###############################################################
clc;

SZ.r = size(IMG,1);
SZ.c = size(IMG,2);
SZ.z = size(IMG,3);

% FIRST RESHAPE IMAGE STACK INTO A SINGLE MATRIX
IM = reshape(IMG,  SZ.r*SZ.c,[]  );


% RUN PRINCIPAL COMPONENTS ANALYSIS
IM = single(IM);

[PC.coef,PC.score,~] = pca(IM');

[PC(2).coef,PC(2).score,~] = pca(IM);


% RESHAPE COEF BACK INTO A STACK
PC(1).imc = reshape( PC(1).coef , SZ.r , SZ.c , [] );
PC(1).imc = PC(1).imc(:,:,1:25); % GET THE FIRST 25 COMPONENTS

PC(2).ims = reshape( PC(2).score , SZ.r , SZ.c , [] );
PC(2).ims = PC(2).ims(:,1:25); % GET THE FIRST 25 COMPONENTS



plotIMG(IMG, PC(1).imc ,'PCA4')

clearvars -except PIX IMG SMIM PC




%###############################################################
%% VIEW FIRST 4 COMPONENTS AFTER GETTING ABSOLUTE VALUE OF MEAN DEVIATION
%###############################################################
clc

I = PC(1).imc;

clear ABIM
for j=1:size(I,3)

    k=I(:,:,j);
    ABIM(:,:,j) = abs(I(:,:,j) - mean(k(:)));

end


plotIMG(IMG, ABIM ,'PCA4');
title('ABSOLUTE VALUE OF PCA MEAN DEVIATION')
cmappy(colormap('winter'));



pause(3);
clearvars -except PIX IMG SMIM PC ABIM




%###############################################################
%% PREVIEW THE FIRST 16 PRINCIPAL COMPONENTS (COEFFICIENT MATRIX)
%###############################################################



% plotIMG(IMG, PC(1).imc ,'PCA16');
% % cmappy(colormap('winter'));
% 
% pause(3);
% clearvars -except PIX IMG SMIM ROX PC ABIM




%###############################################################
%%      CHOOSE FROM FIRST 16 PRINCIPAL COMPONENTS
%###############################################################



clc; close all

% pickPCs(PC(1).imc)

[AXE] = pickPCs(PC(1).imc);
AXE(1) = [];

disp('CHOSEN PCs:'); disp(AXE)

PCI = abs(PC(1).imc(:,:,AXE));

viewstack(rescale(PCI),1)


clearvars -except PIX IMG SMIM PC ABIM PCI




%###############################################################
%%      DISPLAY HISTOGRAM AND BG CUTOFF
%###############################################################


[THRESH] = imhist(IMG, PCI);

clearvars -except PIX IMG SMIM PC ABIM PCI





%###############################################################
%%      GET MEAN MAX-MIN PROJECTION OF RAW IMAGE
%###############################################################
%{

% GET MEAN_MAX INTENSITY PROJECTION IMG
%------------------------------------------
maxI = zeros(size(IMG,1),size(IMG,2),2);

j=1;
for i = 1:10:(size(IMG,3)-10)

    maxI(:,:,j) = (max(IMG(:,:,i:(i+9)),[],3));

j=j+1;
end


MAXI = maxI;



close all; imagesc(MAXI(:,:,1)); colormap hot; 
title('AVERAGE MAX PIXEL INTENSITY OF RAW IMAGE STACK')
pause(2)






% GET MEAN_MIN INTENSITY PROJECTION IMG
%------------------------------------------
minI = zeros(size(IMG,1),size(IMG,2),2);

j=1;
for i = 1:10:(size(IMG,3)-10)

    minI(:,:,j) = (min(IMG(:,:,i:(i+9)),[],3));

j=j+1;
end


MINI = minI;



close all; imagesc(MINI(:,:,1)); colormap hot; 
title('AVERAGE MIN PIXEL INTENSITY OF RAW IMAGE STACK')
pause(2)



% GET MAX - MIN INTENSITY PROJECTION IMG
%------------------------------------------

IMAX = max(MAXI - MINI,[],3);



close all; imagesc(IMAX(:,:,1)); colormap hot; 
title('MAX - MIN PIXEL INTENSITY OF RAW IMAGE STACK')
pause(2)



clearvars -except PIX IMG SMIM PC ABIM PCI IMAX
%}



% GET MAX DIFF OF RAW IMAGE STACK
%------------------------------------------


% IMmin = min(double(IMG),[],3);
% IMmax = max(double(IMG),[],3);
% 
% IMAX = IMmax - IMmin;

nbins = 20;

sz = size(IMG);
t = round(linspace(1,sz(3),nbins));
IM = double(IMG);

IMAX = [];
for i = 1:numel(t)-1
    IMmin = min(IM(:,:, t(i):t(i+1) ) ,[],3);
    IMmax = max(IM(:,:, t(i):t(i+1) ) ,[],3);
    IMAX(:,:,i) = IMmax - IMmin;
end
IMAX = mean(IMAX,3);



close all; imagesc(IMAX); colormap hot
title('AVERAGE PIXEL VARIANCE OF RAW IMAGE STACK')
pause(2)


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX





%###############################################################
%%      GET PIXEL VARIANCE OF RAW IMAGE STACK
%###############################################################


nbins = 20;

sz = size(IMG);
t = round(linspace(1,sz(3),nbins));
IM = double(IMG);

IMV = [];
for i = 1:numel(t)-1
    IMV(:,:,i) = std(IM(:,:, t(i):t(i+1) ) ,[],3);
end
IMV = mean(IMV,3);



close all; imagesc(IMV); colormap hot
title('AVERAGE PIXEL VARIANCE OF RAW IMAGE STACK')
pause(2)


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV







%###############################################################
%%      CREATE COMPOSITE IMAGE USING COMBINATION OF ABOVE
%###############################################################


disp('IMG');  imstats(IMG);     % RAW IMAGE STACK
disp('SMIM'); imstats(SMIM);    % GAUSSIAN SMOOTHED VERSION OF IMG
disp('IMAX'); imstats(IMAX);    % MEAN MAX PIXEL INTENSITY OF IMG
disp('ABIM'); imstats(ABIM);    % ABSOLUTE MEAN DEVIATION OF ALL PCs
disp('PCI');  imstats(PCI);     % CHOSEN PRINCIPAL COMPONENTS
disp('IMV');  imstats(IMV);     % STDEV OF EACH IMG PIXEL ALONG 3RD DIM


%############   PLOT ALL 6 COMPOSITE OPTIONS   ################
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.02 .56 .3 .4],'Color','none'); axis off; hold on;
ax02 = axes('Position',[.35 .56 .3 .4],'Color','none'); axis off; hold on;
ax03 = axes('Position',[.67 .56 .3 .4],'Color','none'); axis off; hold on;
ax04 = axes('Position',[.02 .06 .3 .4],'Color','none'); axis off; hold on;
ax05 = axes('Position',[.35 .06 .3 .4],'Color','none'); axis off; hold on;
ax06 = axes('Position',[.67 .06 .3 .4],'Color','none'); axis off; hold on;

axes(ax01); imagesc(mean(IMG,3));  title('RAW IMAGE STACK');
axes(ax02); imagesc(mean(SMIM,3)); title('GAUSSIAN SMOOTHED VERSION OF IMG');
axes(ax03); imagesc(mean(IMAX,3)); title('MEAN MAX PIXEL INTENSITY OF IMG');
axes(ax04); imagesc(mean(ABIM,3)); title('ABSOLUTE MEAN DEVIATION OF ALL PCs');
axes(ax05); imagesc(mean(PCI,3));  title('CHOSEN PRINCIPAL COMPONENTS');
axes(ax06); imagesc(mean(IMV,3));  title('STDEV OF EACH IMG PIXEL ALONG 3RD DIM');

colormap hot
pause(2)





NIM.IMG  = rescale(IMG);
NIM.SMIM = rescale(SMIM);
NIM.IMAX = rescale(IMAX);
NIM.ABIM = rescale(ABIM);
NIM.PCI  = rescale(PCI);
NIM.IMV  = rescale(IMV);


disp('IMG');  imstats(NIM.IMG);     % RAW IMAGE STACK
disp('SMIM'); imstats(NIM.SMIM);    % GAUSSIAN SMOOTHED VERSION OF IMG
disp('IMAX'); imstats(NIM.IMAX);    % MEAN MAX PIXEL INTENSITY OF IMG
disp('ABIM'); imstats(NIM.ABIM);    % ABSOLUTE MEAN DEVIATION OF ALL PCs
disp('PCI');  imstats(NIM.PCI);     % CHOSEN PRINCIPAL COMPONENTS
disp('IMV');  imstats(NIM.IMV);     % STDEV OF EACH IMG PIXEL ALONG 3RD DIM


%############   PLOT ALL 6 COMPOSITE OPTIONS   ################
fh02 = figure('Units','normalized','OuterPosition',[.03 .07 .95 .90],...
              'Color','w','MenuBar','none');
ax11 = axes('Position',[.02 .56 .3 .4],'Color','none'); axis off; hold on;
ax12 = axes('Position',[.35 .56 .3 .4],'Color','none'); axis off; hold on;
ax13 = axes('Position',[.67 .56 .3 .4],'Color','none'); axis off; hold on;
ax14 = axes('Position',[.02 .06 .3 .4],'Color','none'); axis off; hold on;
ax15 = axes('Position',[.35 .06 .3 .4],'Color','none'); axis off; hold on;
ax16 = axes('Position',[.67 .06 .3 .4],'Color','none'); axis off; hold on;

axes(ax11); imagesc(mean(NIM.IMG,3));  title('RAW IMAGE STACK');
axes(ax12); imagesc(mean(NIM.SMIM,3)); title('GAUSSIAN SMOOTHED VERSION OF IMG');
axes(ax13); imagesc(mean(NIM.IMAX,3)); title('MEAN MAX PIXEL INTENSITY OF IMG');
axes(ax14); imagesc(mean(NIM.ABIM,3)); title('ABSOLUTE MEAN DEVIATION OF ALL PCs');
axes(ax15); imagesc(mean(NIM.PCI,3));  title('CHOSEN PRINCIPAL COMPONENTS');
axes(ax16); imagesc(mean(NIM.IMV,3));  title('STDEV OF EACH IMG PIXEL ALONG 3RD DIM');

colormap hot
pause(2)





MAGE = NIM.IMG(:,:,1);
MAGE(:,:,1) = rescale(mean(NIM.IMG,3));
MAGE(:,:,2) = rescale(mean(NIM.SMIM,3));
MAGE(:,:,3) = rescale(mean(NIM.IMAX,3));
MAGE(:,:,4) = rescale(mean(NIM.ABIM,3));
MAGE(:,:,5) = rescale(mean(NIM.PCI,3));
MAGE(:,:,6) = rescale(mean(NIM.IMV,3));



MAGE(:,:,1) = rescale(MAGE(:,:,1).^2)./20;
MAGE(:,:,2) = rescale(MAGE(:,:,2).^2)./20;
MAGE(:,:,3) = rescale(MAGE(:,:,3).^2)./20;

MAGE(:,:,4) = rescale(MAGE(:,:,4).^2);
MAGE(:,:,5) = rescale(MAGE(:,:,5).^2);
MAGE(:,:,6) = rescale(MAGE(:,:,6).^2);


PIC = rescale(mean(MAGE(:,:,[1 2 3 4 5 6]),3));


% ################   SINGLE AXIS ABSOLUTE LOCATION   ################

fh03 = figure('Units','pixels','Position',[100 35 800 750],...
    'Color','w','MenuBar','none');
ax31 = axes('Position',[.06 .06 .9 .9],'Color','none');


imagesc(PIC);
colormap hot
title('COMPOSITE IMAGE FOR AUTOMATED IMAGE SEGENTATION');



clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC














%########################################################################
%%              PERFORM IMAGE SEGMENTATION
%########################################################################

clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC



% SEGMENT IMAGE
%--------------------------------------------------------
% [BWMASK,BWRAW] = segIM2(PIC);
[BWMASK,BWRAW] = segIM3(PIC);


close all; p=imagesc(BWMASK); 
p.CData=BWMASK; pause(.2)



% GET REGION PROPERTIES & STATISTICS
%--------------------------------------------------------

IMFO.stats = regionprops(BWMASK);

[IMFO.bi,IMFO.labs,IMFO.n,IMFO.a] = bwboundaries(BWMASK,'noholes');





% PLOT BOUNDING COORDINATES AROUND ROIs
%--------------------------------------------------------

clc; close all
fh1 = figure('Units','normalized','OuterPosition',[.03 .06 .6 .8],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none',...
            'YDir','reverse','XColor','none','YColor','none'); hold on
ax2 = axes('Position',[.06 .06 .9 .9],'Color','none'); hold on

axis(ax1)
p=imagesc(BWMASK);p.CData=BWMASK; 
axis tight; colormap bone; 

axis(ax2)
for i = 1:numel(IMFO.stats)
    scatter(IMFO.bi{i}(:,2),IMFO.bi{i}(:,1),'.')
    hold on
end

pause(2)









% PULL OUT SOME STATS SO WE CAN DETERMINE IF ANYTHING
% NEEDS TO BE REMOVED OR RESHAPED
%--------------------------------------------------------
Area = [IMFO.stats.Area]';

c=[IMFO.stats.Centroid]; 
cx = c(1:2:end);
cy = c(2:2:end);
Centroid = [cx',cy'];

b = [IMFO.stats.BoundingBox];
br = b(1:4:end);
bc = b(2:4:end);
bw = b(3:4:end);
bh = b(4:4:end);
BBox = [bw',bh' (bw./bh)' (bh./bw)'];
%--------------





% ESTABLISH FILTERING PARAMETERS
%--------------------------------------------------------
AREA_FILTER = [12 , 400];      % <<<<<<<<< USER SHOULD ENTER THIS <<<<<<<<<<

% BOX_FILTER  = [5 5 1/4 1/4]; % <<<<<<<<< USER SHOULD ENTER THIS <<<<<<<<<<





% DETERMINE IF MORPHOLOGY STATS MEET FILTER CRITERIA
%--------------------------------------------------------

TOO.SMALL = AREA_FILTER(1) > Area;

TOO.BIG   = AREA_FILTER(2) < Area;

% TOO.BOXY  = sum((BBox < BOX_FILTER),2);


% FAIL = TOO.SMALL | TOO.BIG | TOO.BOXY;
FAIL = TOO.SMALL | TOO.BIG;
F = find(FAIL);



nMinExpectedROIs = 5;  % <<<<<<<<<<<< USER SHOULD ENTER THIS <<<<<<<<<<<<

nMaxExpectedROIs = 40; % <<<<<<<<<<<< USER SHOULD ENTER THIS <<<<<<<<<<<<


nROIs = numel(Area);   
fprintf('Total ROI count (first-pass): %0.f \n\n',nROIs)


nSmall = sum(TOO.SMALL);
fprintf('Number of ROIs below threshold: %0.f \n\n',nSmall)


nBig   = sum(TOO.BIG);
fprintf('Number of ROIs above threshold: %0.f \n\n',nBig)






% REMOVE ROIS THAT DID NOT PASS ALL THRESH TESTS
%--------------------------------------------------------

BW = IMFO.labs;
for i = 1:numel(F)
    BW(BW == F(i)) = 0;
end

BWMASK = BW > 0;





% AGAIN GET REGION PROPERTIES & STATISTICS
%--------------------------------------------------------

IMFO.stats = regionprops(BWMASK);

[IMFO.bi,IMFO.labs,IMFO.n,IMFO.a] = bwboundaries(BWMASK,'noholes');





% AGAIN PLOT BOUNDING COORDINATES AROUND ROIs
%--------------------------------------------------------

clc; close all
fh1 = figure('Units','normalized','OuterPosition',[.03 .06 .6 .8],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none',...
            'YDir','reverse','XColor','none','YColor','none'); hold on
ax2 = axes('Position',[.06 .06 .9 .9],'Color','none'); hold on

axis(ax1)
p=imagesc(BWMASK);p.CData=BWMASK; 
axis tight; colormap bone; 

axis(ax2)
for i = 1:numel(IMFO.stats)
    scatter(IMFO.bi{i}(:,2),IMFO.bi{i}(:,1),'.')
    hold on
end

pause(2)





clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS






%% SHOW ROI ACTIVITY AND GET MEAN ACTIVITY IN EACH ROI
%---------------------------------------------

I = rescale(IMG) .* BWMASK;

q = quantile(I(:),[.0001 .9999]);

close all; figure; a=axes; colormap(bone)

p = imagesc(I(:,:,1));   a.CLim=q;

for i = 1:size(I,3)
    p.CData = I(:,:,i);  pause(.04)
end

pause(1)



%% GET MEAN ACTIVITY IN EACH ROI
%--------------------------------------------------------

IM = rescale(IMG);
MUJ=[];
for i = 1:IMFO.n

    msk = IMFO.labs==i;

    for j = 1:size(IM,3)

        IMJ = IM(:,:,j);

        MUJ(i,j) = mean(IMJ(msk));
    end
end

ROIS = MUJ';

minROI = min(ROIS);

ROIS = ROIS - minROI;

ROIS = rescale(ROIS);



clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS




%% DETERMINE IF ACTIVITY IS SIMPLY RUNUP OR RUNDOWN
%--------------------------------------------------------


nbins = 9;

t = round(linspace(1,size(ROIS,2),nbins));

StartMu = mean(  ROIS( t(2):t(3)         ,:)  );
EndMu   = mean(  ROIS( t(end-3):t(end-2) ,:)  );

IRUN = StartMu - EndMu;
IRUNmu = mean(IRUN);
IRUNsd = std(IRUN);

RAN = (IRUN > (IRUNsd*2 + IRUNmu)) | (IRUN < (IRUNsd*-2 + IRUNmu));


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS RAN





%% AGAIN REMOVE ROIS THAT DID NOT PASS ALL THRESH TESTS
%--------------------------------------------------------

ROIS(:,RAN) = [];

F = find(RAN);

BW = IMFO.labs;
for i = 1:numel(F)
    BW(BW == F(i)) = 0;
end

BWMASK = BW > 0;

IMFO.stats = regionprops(BWMASK);
[IMFO.bi,IMFO.labs,IMFO.n,IMFO.a] = bwboundaries(BWMASK,'noholes');


% AGAIN PLOT BOUNDING COORDINATES AROUND ROIs
%--------------------------------------------------------
clc; close all
fh1 = figure('Units','normalized','OuterPosition',[.03 .06 .6 .8],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none',...
            'YDir','reverse','XColor','none','YColor','none'); hold on
ax2 = axes('Position',[.06 .06 .9 .9],'Color','none'); hold on
axis(ax1); p=imagesc(BWMASK);p.CData=BWMASK; 
axis tight; colormap bone; axis(ax2);
for i = 1:numel(IMFO.stats)
    scatter(IMFO.bi{i}(:,2),IMFO.bi{i}(:,1),'.');
    hold on
end
pause(2)
clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS






%% SHOW ROI ACTIVITY AND GET MEAN ACTIVITY IN EACH ROI
%---------------------------------------------

I = rescale(IMG) .* BWMASK;

q = quantile(I(:),[.0001 .9999]);

close all; figure; a=axes; colormap(bone)

p = imagesc(I(:,:,1));   a.CLim=q;

for i = 1:size(I,3)
    p.CData = I(:,:,i);  pause(.04)
end

pause(1)




































%%  PLOT MEAN ACTIVITY IN EACH ROI
%--------------------------------------------------------
clc; close all;
fh1 = figure('Units','pixels','Position',[10 35 1300 500],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
ax1.YLim = [0 1]; hold on

ph = plot(ROIS(:,1:3),'k','LineWidth',3);
pause(1)

% x = repmat(1:size(ROIS,1),3,1)';
% y = ROIS(:,1:3);
% for i=1:3
%     %ph1 = scatter(x(:,i),y(:,i),100,'.k'); hold on
%     ph2 = plot(x(:,i),y(:,i)); hold on
% end

for i = 4:size(ROIS,2)

    ph(1).YData = ph(2).YData;
    ph(2).YData = ph(3).YData;
    ph(3).YData = ROIS(:,i);
    ax1.YLim = [0 1];
    pause(.6)

end


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS



%########################################################################
%%  PLOT MEAN ACTIVITY FOR ALL ROIs
%########################################################################


clc; close all;
fh1 = figure('Units','normalized','Position',[.05 .08 .88 .85],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');

n = size(ROIS,2);
f = size(ROIS,1);

R = ROIS + repmat(1:n,f,1);

ph = plot(R,'LineWidth',3);
pause(1)

ROITABLE = table(ROIS);



fh1 = figure('Units','pixels','Position',[10 35 1300 500],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
ax1.YLim = [0 1]; hold on



ph = plot(ROIS,'LineWidth',5);
pause(1)


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS ROITABLE



%########################################################################
%%  EXPORT ROI ACTIVITY TRACES TO SPREADSHEET
%########################################################################
% d = char(datetime('now'));
% d = regexprep(d,':','_');
% d = regexprep(d,' ','_');
% d = regexprep(d,'-','_');

NIM.IMG = uint8(rescale(NIM.IMG).*255);
NIM.SMIM = uint8(rescale(NIM.SMIM).*255);
NIM.ABIM = uint8(rescale(NIM.ABIM).*255);



[path,name] = fileparts(PIX.info.Filename{1});

save(['ROI_ANALYSIS_' name '.mat'],'NIM','PIX','PC','MAGE','PIC','BWMASK',...
'IMFO','ROIS','ROITABLE');

writetable(ROITABLE,['ROI_ANALYSIS_' name '.csv'])






%###############################################################
%%      HAND-CLICK TO CHOOSE ROIS TO KEEP
%###############################################################

clc; close all


[AXE] = handpickroi(ROIS);


disp('CHOSEN ROIs:'); disp(AXE)


ROIL = ROIS(:,AXE);


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS ROITABLE ROIL

%########################################################################
%%  PLOT MEAN ACTIVITY FOR ALL ROIs
%########################################################################


clc; close all;
fh1 = figure('Units','normalized','Position',[.05 .08 .88 .85],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');

n = size(ROIL,2);
f = size(ROIL,1);

R = ROIL + repmat(1:n,f,1);

ph = plot(R,'LineWidth',3);
pause(1)

ROITABLE = table(ROIL);



fh1 = figure('Units','pixels','Position',[10 35 1300 500],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
ax1.YLim = [0 1]; hold on



ph = plot(ROIL,'LineWidth',5);
pause(1)


clearvars -except PIX IMG SMIM PC ABIM PCI IMAX IMV NIM MAGE PIC...
BWMASK IMFO ROIS ROITABLE ROIL


