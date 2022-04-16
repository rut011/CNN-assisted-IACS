
clc;clear all

addpath(pwd);
datadir = 'I:\VeriCells03282022\Monocytes\TDMSData'; % data directory 
storedir = 'I:\VeriCells03282022\Monocytes'; %store directory
foldername = 'ReconImageGallary';
cd(datadir)
data_listing = dir('*.tdms');
for k =1:length(data_listing)
% k = 1;
foldername = 'ReconImageGallary';
foldername = [foldername,num2str(k)];

cd (storedir);
    if ~(exist(foldername))
        mkdir(foldername);
    end
cd(datadir)
SF = 1; %sampling factor
SFDisplay = 10;%sampling factor for waveform plot
disp(['processing file No.',num2str(k),'/',num2str(length(data_listing))]);

tic
Data = convertTDMS(0,data_listing(k).name);

toc

%% parameters
VisualFlag = 0;
VisualStatFlag = 0;
SaveImgFlag = 1;
ImageIdx = 12;
FL1Thres = 0;
FL2Thres = 0;

%% Reorganize Image Data
ImageNum = length(Data.Data.MeasuredData(5).Data);
% ImageNum = 100;
ImageHeight = Data.Data.MeasuredData(6).Data(1);
ImageWidth = Data.Data.MeasuredData(7).Data(1);
TransSignal = -Data.Data.MeasuredData(8).Data/32.768;
FL1Signal = -Data.Data.MeasuredData(9).Data/32.768;
FL2Signal = -Data.Data.MeasuredData(10).Data/32.768;
ImageStack = [];
clear Data
for n = 1:ImageNum
    TempSignalTrans = TransSignal(1+(n-1)*ImageHeight*ImageWidth:n*ImageHeight*ImageWidth);
    TempSignalFL1 = FL1Signal(1+(n-1)*ImageHeight*ImageWidth:n*ImageHeight*ImageWidth);
    TempSignalFL2 = FL2Signal(1+(n-1)*ImageHeight*ImageWidth:n*ImageHeight*ImageWidth);
    ImageStack(n).TransImg = reshape(TempSignalTrans,[ImageWidth ImageHeight]);
    ImageStack(n).FL1Img = reshape(TempSignalFL1,[ImageWidth ImageHeight]);
    ImageStack(n).FL2Img = reshape(TempSignalFL2,[ImageWidth ImageHeight]);
    
end
clear TransSignal
clear FL1Signal
clear FL2Signal
for ImageIdx = 1:ImageNum
    disp(['processing No. ',num2str(ImageIdx),'/',num2str(ImageNum)]);


    %visualize
    if VisualFlag == 1
        figure(1);clf;
        imagesc(ImageStack(ImageIdx).TransImg);colormap(gray);colorbar;
        title(['Tranmission No.', num2str(ImageIdx)]);
        figure(2);clf;
        imagesc(ImageStack(ImageIdx).FL1Img);
        colorMap = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];
        colormap(colorMap);colorbar;
        title(['FL1 No.', num2str(ImageIdx)]);
        figure(3);clf;
        imagesc(ImageStack(ImageIdx).FL2Img);colormap(jet);colorbar;
        title(['Speed No.', num2str(ImageIdx)]);
    end

    %% process
%     TransScanRowstd = sum(abs(diff(ImageStack(ImageIdx).TransImg,[],1)));
%     [~,Transidx] = sort(TransScanRowstd);
%     TransScanRowBkg = mean(ImageStack(ImageIdx).TransImg(:,Transidx(1:10)),2);
%     TransImageBkgSub = ImageStack(ImageIdx).TransImg - repmat(TransScanRowBkg,1,size(ImageStack(ImageIdx).TransImg,2));
%     TransImageBkgSub = TransImageBkgSub - min(TransImageBkgSub,[],'all');
    TransImageBkgSub = ImageStack(ImageIdx).TransImg;
    TransImageBkgSub = TransImageBkgSub - min(TransImageBkgSub,[],'all');
    FL1ImageThres = ImageStack(ImageIdx).FL1Img;
    FL2ImageThres = ImageStack(ImageIdx).FL2Img;
%     FL1ImageThres (FL1ImageThres <=FL1Thres) = 0;

    %visualize
    if VisualFlag == 1
        figure(4);clf;
        imagesc(TransImageBkgSub);colormap(gray);colorbar;
        title(['Tranmission BkgSub No.', num2str(ImageIdx)]);

        figure(5);clf;
        imagesc(FL1ImageThres);
        colorMap = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];
        colormap(colorMap);colorbar;
        title(['FL1 Thres No.', num2str(ImageIdx)]);
    end

    %% speed detection 
    % PMT1Seg1Smooth = mean(ImageStack(n).FL2Img,1);
    % PMT1Seg1Smooth = imresize (PMT1Seg1Smooth, [1 100*size(PMT1Seg1Smooth,2)],'bilinear');
%     PMT1Seg1Smooth = reshape(ImageStack(ImageIdx).FL2Img,[],1);
% 
% %     SpeedThres = 410;
%     fs=1/(5e-6)*145;
%     Slide_width=0.00004;
%     sampling_n=8001;
%     sampling_point=1;
%     [ACF,lags,bounds] = autocorr(PMT1Seg1Smooth,sampling_n-1);
%     ACF_smooth = smooth(ACF,5);
%     [Peak,locsPeak]=findpeaks(ACF_smooth,'MinPeakDistance',500);
%     [peakloc]=find(locsPeak>2000&locsPeak<6000);
%     if isempty(peakloc)
%         SpeedCal = 0.17;
%     else 
%         [val, peakidx] = max(Peak(peakloc));
%         sec_max=locsPeak(peakloc(peakidx))-1;
%         SpeedCal = 32e-6/(sec_max/fs);
%     end

    %visualize
    if VisualFlag == 1
        figure(7),plot(lags(1:end),ACF_smooth(1:end));hold on;plot(sec_max,max(Peak(peakloc)),'r*');
        hold off;title('Autocorrelation');xlabel('sampling points');ylabel('ACF');
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(2)-4000,ylim(2)-0.1,['Detected speed: ', num2str(SpeedCal),'m/s'],'Color', 'r');
        xline(2000,'--','Color','g');
        xline(6000,'--','Color','g');
    end


    %% resize images
    SpeedCal = 0.15;
    xpixel = 0.26;
    ypixel = SpeedCal * 5;
    xFOV = xpixel*size(TransImageBkgSub,1);
    yFOV = ypixel*size(TransImageBkgSub,2);
    
    ImgresizeTrans = imresize(TransImageBkgSub, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
    ImgresizeFL1 = imresize(ImageStack(ImageIdx).FL1Img, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
    ImgresizeFL1Thres = imresize(FL1ImageThres, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
    MaskedFL1Pixel = nonzeros(ImgresizeFL1Thres);
    ImgresizeFL2 = imresize(ImageStack(ImageIdx).FL2Img, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
    ImgresizeFL2Thres = imresize(FL2ImageThres, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
    MaskedFL2Pixel = nonzeros(ImgresizeFL2Thres);

    %visualize
    if VisualFlag == 1
        Fig9=figure(9);clf;
        imagesc(ImgresizeTrans);colormap(gray);colorbar;ylabel('scanning [pixel]');xlabel('flow [pixel]');axis equal
        title(['Trans. speed resized No.', num2str(ImageIdx),' [pixel size: 0.25um]']);xlabel('flow pixel');ylabel('scanning pixel');
        Fig10=figure(10);clf;
        imagesc(ImgresizeFL1);colormap(colorMap);colorbar;ylabel('scanning [pixel]');xlabel('flow [pixel]');axis equal
        title(['FL1 speed resized No.', num2str(ImageIdx),' [pixel size: 0.25um]']);xlabel('flow pixel');ylabel('scanning pixel');
        Fig11=figure(11);clf;
        imagesc(ImgresizeFL1Thres);colormap(colorMap);colorbar;ylabel('scanning [pixel]');xlabel('flow [pixel]');axis equal
        title(['FL1Thres speed resized No.', num2str(ImageIdx),' [pixel size: 0.25um]']);xlabel('flow pixel');ylabel('scanning pixel');
    end

    %% process tranmission image
    % tic
    Iblur1 = imgaussfilt(ImgresizeTrans,1.5);
    [Gmag, Gdir] = imgradient(Iblur1);
    Gmagmean = median(Gmag,'all');
    Gmag(Gmag < Gmagmean) = 0;
    [~,threshold] = edge(Gmag,'Sobel');
    fudgeFactor = 0.16;
    BWs = edge(Iblur1,'Sobel',threshold * fudgeFactor);

    %dilate
    se90 = strel('line',4,90);
    se0 = strel('line',4,0);
    BWsdil = imdilate(BWs,[se90 se0]);
    %fill hole
%     BWsdil = imbinarize(Gmag);
    BWdfill = imfill(BWsdil,'holes');

    %smooth object
    seD = strel('diamond',2);
    BWsmooth = imerode(BWdfill ,seD);
    % BWfinal = imerode(BWfinal,seD);

    %remove small object
    se = strel('disk',5);
    afterOpening = imopen(BWsmooth,se);
    BWfinal = bwareaopen(afterOpening, 450);
    
    %remove object touching border
    BWfinal = imclearborder(BWfinal);

    MaskedTrans = ImgresizeTrans .* BWfinal;
    MaskedTransGmag = Gmag .* BWfinal;
    ImgresizeFL1Thres =  ImgresizeFL1Thres .* BWfinal;
    ImgresizeFL2Thres =  ImgresizeFL2Thres .* BWfinal;
%     ImgresizeFL1Thres (ImgresizeFL1Thres<0) = 0;
    % [xidx, yidx] = find(MaskedTrans);
    MaskedTransPixel = nonzeros(MaskedTrans);
    %gradient RMS calculation
    MaskedTransGmagPixel = nonzeros(MaskedTransGmag);
    Gradientmean = mean(MaskedTransGmagPixel,'all');
    % GradientVar = abs(MaskedTransPixel-Gradientmean);
    GradientRMS = rms(MaskedTransGmagPixel);

    %visualize
    if VisualFlag ==1
        figure(12);clf;
        subplot(241);
    %     imagesc(ImgresizeTrans);colormap(gray);colorbar;ylabel('scanning [pixel]');xlabel('flow [pixel]');axis equal
        imshow(mat2gray(ImgresizeTrans));
        title(['Trans. speed resized No.', num2str(ImageIdx),' [pixel size: 0.25um]']);xlabel('flow pixel');ylabel('scanning pixel');
        subplot(242);
        imshow(BWs)
        title('Binary Gradient Mask');
        subplot(243);
        imshow(BWsdil);
        title('Dilated Gradient Mask');
        subplot(244);
        imshow(BWdfill)
        title('Binary Image with Filled Holes')
        subplot(245);
        imshow(BWfinal);
        title('Segmented Image');
        subplot(246);
        imshow(labeloverlay(mat2gray(ImgresizeTrans),BWsmooth))
        title('smooth Mask Over Original Image');
        subplot(247);
        imshow(labeloverlay(mat2gray(ImgresizeTrans),BWfinal))
        title('final Mask Over Original Image');
        subplot(248);
        imshow(mat2gray(Gmag));
        title('Image Gradient Magnitude');
    end
    %% bounding box
    if (ImgresizeFL1Thres == 0)
        FL1AspectRatio = 0;
        FL1width = 0;
        FL1height = 0;
        FL1Area = 0;
    else
        [FL1MaskX, FL1Masky]= find(ImgresizeFL1Thres);
        FL1width = max(FL1MaskX,[],'all')-min(FL1MaskX,[],'all');
        FL1height = max(FL1Masky,[],'all')-min(FL1Masky,[],'all');
        FL1Area = length(FL1MaskX);
        if FL1width < FL1height
            FL1AspectRatio =  FL1width/FL1height;
        else
            FL1AspectRatio =  FL1height/FL1width;
        end
    end
    
    if (ImgresizeFL2Thres == 0)
        FL2AspectRatio = 0;
        FL2width = 0;
        FL2height = 0;
        FL2Area = 0;
    else
        [FL2MaskX, FL2Masky]= find(ImgresizeFL2Thres);
        FL2width = max(FL2MaskX,[],'all')-min(FL2MaskX,[],'all');
        FL2height = max(FL2Masky,[],'all')-min(FL2Masky,[],'all');
        FL2Area = length(FL2MaskX);
        if FL2width < FL2height
            FL2AspectRatio =  FL2width/FL2height;
        else
            FL2AspectRatio =  FL2height/FL2width;
        end
    end

    if (sum(MaskedTrans,'all') == 0)
        TransAspectRatio = 0;
        Transwidth = 0;
        Transheight = 0;
        TransArea = 0;
    else
        [TransMaskX, TransMasky]= find(MaskedTrans);
        Transwidth = max(TransMaskX,[],'all')-min(TransMaskX,[],'all');
        Transheight = max(TransMasky,[],'all')-min(TransMasky,[],'all');
        TransArea = length(TransMaskX);
        if Transwidth < Transheight
            TransAspectRatio =  Transwidth/Transheight;
        else
            TransAspectRatio =  Transheight/Transwidth;
        end
    end
    %% save information
    ImageStack(ImageIdx).TransImgSub = single(TransImageBkgSub);
    ImageStack(ImageIdx).FL1ImgThres = single(FL1ImageThres);
    ImageStack(ImageIdx).Speed = SpeedCal;
    ImageStack(ImageIdx).TransImgResize = single(ImgresizeTrans);
    ImageStack(ImageIdx).TransImgResizeMasked = uint8(mat2gray(MaskedTrans)*255) + uint8(100 * ~BWfinal);
    ImageStack(ImageIdx).FL1ImgResize = single(ImgresizeFL1);
    ImageStack(ImageIdx).FL1ThresImgResize = uint8(((ImgresizeFL1Thres-min(ImgresizeFL1Thres,[],'all'))./(max(ImgresizeFL1Thres,[],'all')-min(ImgresizeFL1Thres,[],'all'))).* BWfinal*255) + uint8(0 * ~BWfinal);
%     ImageStack(ImageIdx).FL1ThresImgResize = uint8(mat2gray(ImgresizeFL1Thres)*255) + uint8(0 * ~BWfinal);
    ImageStack(ImageIdx).FL2ImgResize = single(ImgresizeFL2);
    ImageStack(ImageIdx).FL2ThresImgResize = uint8(((ImgresizeFL2Thres-min(ImgresizeFL2Thres,[],'all'))./(max(ImgresizeFL2Thres,[],'all')-min(ImgresizeFL2Thres,[],'all'))).* BWfinal*255) + uint8(0 * ~BWfinal);
    
    ImageStack(ImageIdx).ID = ImageIdx;
    ImageStack(ImageIdx).Folder = k;
    ImageStack(ImageIdx).TransMaxInt = max(MaskedTrans,[],'all');
    ImageStack(ImageIdx).TransSumInt = sum(MaskedTrans,'all');
    ImageStack(ImageIdx).TransMeanInt = mean(MaskedTransPixel,'all');
    ImageStack(ImageIdx).TransGradientRMS = GradientRMS;
    ImageStack(ImageIdx).Contrast = max(MaskedTrans,[],'all') - min(MaskedTrans,[],'all');
    ImageStack(ImageIdx).TransArea = TransArea;
    ImageStack(ImageIdx).TransWidth = Transwidth;
    ImageStack(ImageIdx).TransHeight = Transheight;
    ImageStack(ImageIdx).TransAspectRatio = TransAspectRatio;
    ImageStack(ImageIdx).FL1MaxInt = max(ImgresizeFL1,[],'all');
    ImageStack(ImageIdx).FL1SumInt = sum(ImgresizeFL1,'all');
    ImageStack(ImageIdx).FL1MeanInt = mean(ImgresizeFL1,'all');
    ImageStack(ImageIdx).FL1ThresSumInt = sum(ImgresizeFL1Thres,'all');
    ImageStack(ImageIdx).FL1ThresMeanInt = mean(ImgresizeFL1Thres,'all');
    ImageStack(ImageIdx).FL1ThresMaxInt = max(ImgresizeFL1Thres,[],'all');
    ImageStack(ImageIdx).FL1Area = FL1Area;
    ImageStack(ImageIdx).FL1Width = FL1width;
    ImageStack(ImageIdx).FL1Height = FL1height;
    ImageStack(ImageIdx).FL1AspectRatio = FL1AspectRatio;
    
    ImageStack(ImageIdx).FL2MaxInt = max(ImgresizeFL2,[],'all');
    ImageStack(ImageIdx).FL2SumInt = sum(ImgresizeFL2,'all');
    ImageStack(ImageIdx).FL2MeanInt = mean(ImgresizeFL2,'all');
    ImageStack(ImageIdx).FL2ThresSumInt = sum(ImgresizeFL2Thres,'all');
    ImageStack(ImageIdx).FL2ThresMeanInt = mean(ImgresizeFL2Thres,'all');
    ImageStack(ImageIdx).FL2ThresMaxInt = max(ImgresizeFL2Thres,[],'all');
    ImageStack(ImageIdx).FL2Area = FL2Area;
    ImageStack(ImageIdx).FL2Width = FL2width;
    ImageStack(ImageIdx).FL2Height = FL2height;
    ImageStack(ImageIdx).FL2AspectRatio = FL2AspectRatio;
    
    % toc
end

%% Generate feature csv and mat
ID = extractfield(ImageStack,'ID');
Folder = extractfield(ImageStack,'Folder');
TransMaxInt = extractfield(ImageStack,'TransMaxInt'); 
TransSumInt = extractfield(ImageStack,'TransSumInt'); 
TransMeanInt = extractfield(ImageStack,'TransMeanInt'); 
TransGradientRMS = extractfield(ImageStack,'TransGradientRMS'); 
Contrast = extractfield(ImageStack,'Contrast'); 
TransArea = extractfield(ImageStack,'TransArea'); 
TransWidth = extractfield(ImageStack,'TransWidth'); 
TransHeight = extractfield(ImageStack,'TransHeight');
TransAspectRatio = extractfield(ImageStack,'TransAspectRatio');
FL1MaxInt = extractfield(ImageStack,'FL1MaxInt'); 
FL1SumInt = extractfield(ImageStack,'FL1SumInt'); 
FL1MeanInt = extractfield(ImageStack,'FL1MeanInt'); 
FL1ThresSumInt = extractfield(ImageStack,'FL1ThresSumInt'); 
FL1ThresMeanInt = extractfield(ImageStack,'FL1ThresMeanInt');
FL1ThresMaxInt = extractfield(ImageStack,'FL1ThresMaxInt');
FL1Area = extractfield(ImageStack,'FL1Area'); 
FL1Width = extractfield(ImageStack,'FL1Width'); 
FL1Height = extractfield(ImageStack,'FL1Height'); 
FL1AspectRatio = extractfield(ImageStack,'FL1AspectRatio');

FL2MaxInt = extractfield(ImageStack,'FL2MaxInt'); 
FL2SumInt = extractfield(ImageStack,'FL2SumInt'); 
FL2MeanInt = extractfield(ImageStack,'FL2MeanInt'); 
FL2ThresSumInt = extractfield(ImageStack,'FL2ThresSumInt'); 
FL2ThresMeanInt = extractfield(ImageStack,'FL2ThresMeanInt');
FL2ThresMaxInt = extractfield(ImageStack,'FL2ThresMaxInt');
FL2Area = extractfield(ImageStack,'FL2Area'); 
FL2Width = extractfield(ImageStack,'FL2Width'); 
FL2Height = extractfield(ImageStack,'FL2Height'); 
FL2AspectRatio = extractfield(ImageStack,'FL2AspectRatio'); 

CSV.data = [ID',Folder',TransMaxInt',TransSumInt',TransMeanInt',TransGradientRMS',Contrast',...
    TransArea',TransWidth',TransHeight',TransAspectRatio',FL1MaxInt',FL1SumInt',FL1MeanInt',...
    FL1ThresSumInt',FL1ThresMeanInt',FL1ThresMaxInt',FL1Area',FL1Width',FL1Height',FL1AspectRatio',...
    FL2MaxInt',FL2SumInt',FL2MeanInt',...
    FL2ThresSumInt',FL2ThresMeanInt',FL2ThresMaxInt',FL2Area',FL2Width',FL2Height',FL2AspectRatio'];

csv.colheaders = {'ID','Folder','TransMaxInt','TransSumInt','TransMeanInt','TransGradientRMS','Contrast',...
    'TransArea','TransWidth','TransHeight','TransAspectRatio','FL1MaxInt','FL1SumInt','FL1MeanInt',...
    'FL1ThresSumInt','FL1ThresMeanInt','FL1ThresMaxInt','FL1Area','FL1Width','FL1Height','FL1AspectRatio',...
    'FL2MaxInt','FL2SumInt','FL2MeanInt',...
    'FL2ThresSumInt','FL2ThresMeanInt','FL2ThresMaxInt','FL2Area','FL2Width','FL2Height','FL2AspectRatio'};

CSVTable = array2table(CSV.data,'VariableNames',{'ID','Folder','TransMaxInt','TransSumInt','TransMeanInt','TransGradientRMS','Contrast',...
    'TransArea','TransWidth','TransHeight','TransAspectRatio','FL1MaxInt','FL1SumInt','FL1MeanInt',...
    'FL1ThresSumInt','FL1ThresMeanInt','FL1ThresMaxInt','FL1Area','FL1Width','FL1Height','FL1AspectRatio',...
    'FL2MaxInt','FL2SumInt','FL2MeanInt',...
    'FL2ThresSumInt','FL2ThresMeanInt','FL2ThresMaxInt','FL2Area','FL2Width','FL2Height','FL2AspectRatio'});

cd(storedir);
fname = [data_listing(k).name(1:end-5),'features','.mat'];
save(fname,'CSV','-v7.3');

fnameCSV = [data_listing(k).name(1:end-5),'features','.csv'];
writetable(CSVTable,fnameCSV);




%%
if VisualStatFlag == 1
FL1ThresSumInt = extractfield(ImageStack,'FL1ThresSumInt');
figure(1);clf;
histogram(FL1ThresSumInt,1500);
xlim([30000, 3e+5]);
set(gca,'xscale','log')
title('FL1 Intensity Sum');

FL1MaxInt = extractfield(ImageStack,'FL1MaxInt');
figure(2);clf;
histogram(FL1MaxInt,500);

Aspect = extractfield(ImageStack,'FL1AspectRatio');
figure(3);clf;
histogram(Aspect,500);

AreaAll = extractfield(ImageStack,'FL1Area');

figure(4);clf;
scatter(FL1ThresSumInt,Aspect,'.');
xlabel('FL1 Intensity Sum'); ylabel('aspect ratio');title('aspect ratio vs. FL1 Intensity Sum');
[value AsIdx] =find(Aspect>0.5);

figure(5);clf;
scatter(AreaAll, Aspect,'.');
xlabel('area');ylabel('aspect ratio');title('are vs. aspect ratio');
end
%% save images
if SaveImgFlag == 1
    colorMapG = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];
    colorMapR = [linspace(0,1,256)', zeros(256,1), zeros(256,1)];
    store_dir = [storedir,'\',foldername];
    if ImageNum > 500
        ImageNum1 = ImageNum;
    else
        ImageNum1 = ImageNum;
    end
    parfor ImageIdx = 1:ImageNum1
        disp(['processing No. ',num2str(ImageIdx),'/',num2str(ImageNum1)]);
        cd(store_dir);
        TransName = [num2str(ImageIdx),'_TransRaw','.tif'];
        TransMaskedName = [num2str(ImageIdx),'_TransMasked','.tif'];
        FL1Name = [num2str(ImageIdx),'_FL1Raw','.tif'];
        FL2Name = [num2str(ImageIdx),'_FL2Raw','.tif'];
        FL1ThresName = [num2str(ImageIdx),'_FL1Masked','.tif'];
        FL2ThresName = [num2str(ImageIdx),'_FL2Masked','.tif'];
        TransProcess = ImageStack(ImageIdx).TransImgResizeMasked;
        Trans = mat2gray(ImageStack(ImageIdx).TransImgResize);
        FL1 = mat2gray(ImageStack(ImageIdx).FL1ImgResize);
        FL2 = mat2gray(ImageStack(ImageIdx).FL2ImgResize);
%         FL1Thres = mat2gray(ImageStack(ImageIdx).FL1ThresImgResize);
%         FL2Thres = mat2gray(ImageStack(ImageIdx).FL2ThresImgResize);
        imwrite(uint8(FL1*255),colorMapG, FL1Name);
        imwrite(uint8(FL2*255),colorMapR, FL2Name);
%         imwrite(uint8(FL1Thres*255),colorMap, FL1ThresName);
        imwrite(uint8(Trans*255),gray, TransName);
        imwrite(TransProcess,gray,TransMaskedName);

    end
end
end
%% save data
% cd(storedir);
% fname = [data_listing(k).name(1:end-5),'.mat'];
% save(fname,'ImageStack','-v7.3');