clc;clear
currentdir = pwd;
addpath(genpath(pwd));
datadir = 'I:\CHO-GFP\TDMS'; % data directory 
csvdir = 'I:\CHO-GFP'; %csv directory
storedir = 'I:\CHO-GFP\CHODatav2\HighGFP'; %store directory
CSVfile = 'HighGFP.csv';

foldername = CSVfile(1:end-4);
cd (storedir);

if ~(exist('ImageRaw','dir'))
    mkdir ImageRaw;
end

if ~(exist('ImageMask','dir'))
    mkdir ImageMask;
end

if ~(exist('ImageFile','dir'))
    mkdir ImageFile;
end

ImageMatStore = [storedir,'\ImageRaw'];
ImageMaskStore = [storedir,'\ImageMask'];
ImageFileStore = [storedir,'\ImageFile'];


tic
cd(csvdir);
CSVTable = readmatrix(CSVfile);

searchFolderList = CSVTable(:,2);
searchIDList = CSVTable(:,1);
clear ia
[C, ia] = unique(searchFolderList,'first');
C = [C;length(C)+1];
ia = [ia; length(searchFolderList)+1];

for k = 1
    cd(datadir);
    data_listing = dir('*.tdms');
    disp(['processing file No.',num2str(k),'/',num2str(length(data_listing))]);

    tic
    Data = convertTDMS(0,data_listing(k).name);

    toc

    %% parameters
    VisualFlag = 0;
    VisualStatFlag = 0;
    SaveImgFlag = 0;
    ImageIdx = 12;
    FL1Thres = 7;
    colorMap = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];
    ResizeWidth = 80;
    ResizeHeight = 112;

    %% Reorganize Image Data
    ImageNum = length(Data.Data.MeasuredData(5).Data);
    ImageHeight = Data.Data.MeasuredData(6).Data(1);
    ImageWidth = Data.Data.MeasuredData(7).Data(1);
    TransSignal = -Data.Data.MeasuredData(8).Data/32.768;
    FL1Signal = -Data.Data.MeasuredData(9).Data/32.768;
    FL2Signal = Data.Data.MeasuredData(10).Data/32.768;
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
    clear TempSignalTrans
    clear TempSignalFL1
    clear TempSignalFL2
    ImageList = int16(searchIDList(ia(C(k)):ia(C(k+1))-1));
    ImageStack = ImageStack(ImageList);
    for ImageIdx =  1:6000%length(ImageStack)
%     for ImageIdx =  1%:length(ImageStack)
        disp(['processing No. ',num2str(ImageIdx),'/',num2str(length(ImageStack))]);
        colorMap = [zeros(256,1),linspace(0,1,256)', zeros(256,1)];

        %visualize
        if VisualFlag == 1
            figure(1);clf;
            imagesc(ImageStack(ImageIdx).TransImg);colormap(gray);colorbar;
            title(['Tranmission No.', num2str(ImageIdx)]);
            figure(2);clf;
            imagesc(ImageStack(ImageIdx).FL1Img);
            
            colormap(colorMap);colorbar;
            title(['FL1 No.', num2str(ImageIdx)]);
            figure(3);clf;
            imagesc(ImageStack(ImageIdx).FL2Img);colormap(jet);colorbar;
            title(['Speed No.', num2str(ImageIdx)]);
        end

        %% process
        TransScanRowstd = sum(abs(diff(ImageStack(ImageIdx).TransImg,[],1)));
%         [~,Transidx] = sort(TransScanRowstd);
%         TransScanRowBkg = mean(ImageStack(ImageIdx).TransImg(:,Transidx(1:10)),2);
        
        if mean(TransScanRowstd(1:3))>=mean(TransScanRowstd(end-2:end))
            TransScanRowBkg = mean(ImageStack(ImageIdx).TransImg(:,1:3),2);
        else
            TransScanRowBkg = mean(ImageStack(ImageIdx).TransImg(:,end-2:end),2);
        end
        TransImageBkgSub = ImageStack(ImageIdx).TransImg - repmat(TransScanRowBkg,1,size(ImageStack(ImageIdx).TransImg,2));
        TransImageBkgSub = TransImageBkgSub - min(TransImageBkgSub,[],'all');
        
%         TransImageBkgSub = ImageStack(ImageIdx).TransImg - min(ImageStack(ImageIdx).TransImg,[],'all');

        FL1ImageThres = ImageStack(ImageIdx).FL1Img;
%         FL1ImageThres (FL1ImageThres <=FL1Thres) = 0;

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
%         PMT1Seg1Smooth = reshape(ImageStack(ImageIdx).FL2Img,[],1);
% 
%         SpeedThres = 410;
%         fs=1/(5e-6)*145;
%         Slide_width=0.00004;
%         sampling_n=8001;
%         sampling_point=1;
%         [ACF,lags,bounds] = autocorr(PMT1Seg1Smooth,sampling_n-1);
%         ACF_smooth = smooth(ACF,5);
%         [Peak,locsPeak]=findpeaks(ACF_smooth,'MinPeakDistance',500);
%         [peakloc]=find(locsPeak>2000&locsPeak<6000);
%         if isempty(peakloc)
%             SpeedCal = 0.17;
%         else 
%             [val, peakidx] = max(Peak(peakloc));
%             sec_max=locsPeak(peakloc(peakidx))-1;
%             SpeedCal = 32e-6/(sec_max/fs);
%         end

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
        
%         xpixel = 0.27;
%         ypixel = SpeedCal * 5;
%         xFOV = xpixel*size(TransImageBkgSub,1);
%         yFOV = ypixel*size(TransImageBkgSub,2);
% 
%         ImgresizeTrans = imresize(TransImageBkgSub, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
%         ImgresizeFL1 = imresize(ImageStack(ImageIdx).FL1Img, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
%         ImgresizeFL1Thres = imresize(FL1ImageThres, [round(xFOV/0.25) round(yFOV/0.25)],'bilinear');
%         MaskedFL1Pixel = nonzeros(ImgresizeFL1Thres);
%         
%         %image cropping
%         if size(ImgresizeTrans,2)>281
%             ImgresizeTrans = ImgresizeTrans(:,1:281);
%             ImgresizeFL1 = ImgresizeFL1(:,1:281);
%             ImgresizeFL1Thres = ImgresizeFL1Thres(:,1:281);
%             MaskedFL1Pixel = nonzeros(ImgresizeFL1Thres);
%         else
%             
%             ImgresizeTrans =[ImgresizeTrans,...
%                 repmat(ImgresizeTrans(:,end),1,281-size(ImgresizeTrans,2))];
%             ImgresizeFL1 =[ImgresizeFL1,...
%                 repmat(ImgresizeFL1(:,end),1,281-size(ImgresizeFL1,2))];
%             ImgresizeFL1Thres =[ImgresizeFL1Thres,...
%                 repmat(ImgresizeFL1Thres(:,end),1,281-size(ImgresizeFL1Thres,2))];
%             
%         end
        
                
        

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
%         Iblur1 = imgaussfilt(TransImageBkgSub,1.5);
%         [Gmag, Gdir] = imgradient(Iblur1);
%         Gmagmean = median(Gmag,'all');
%         Gmag(Gmag < Gmagmean) = 0;
%         [~,threshold] = edge(Gmag,'Sobel');
%         fudgeFactor = 0.163;
%         BWs = edge(Iblur1,'Sobel',threshold * fudgeFactor);
% 
%         %dilate
%         se90 = strel('line',4,90);
%         se0 = strel('line',4,0);
%         BWsdil = imdilate(BWs,[se90 se0]);
%         %fill hole
%     %     BWsdil = imbinarize(Gmag);
%         BWdfill = imfill(BWsdil,'holes');
% 
%         %smooth object
%         seD = strel('diamond',2);
%         BWsmooth = imerode(BWdfill ,seD);
%         % BWfinal = imerode(BWfinal,seD);
% 
%         %remove small object
%         se = strel('disk',5);
%         afterOpening = imopen(BWsmooth,se);
%         BWfinal = bwareaopen(afterOpening, 450);
% 
%         %remove object touching border
%         BWfinal = imclearborder(BWfinal);

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
        
        %Image Mask Generation
%         threshLevel = multithresh(mat2gray(ImgresizeFL1),2);
%         seg_I = imquantize(mat2gray(ImgresizeFL1),threshLevel);       
%         RGB = label2rgb(double(BWfinal)+seg_I-1,[0 0 1;0 1 0;1 0 0],'k');
        
        %resize images and masks for training
%         ImgresizeTrans = imresize(ImgresizeTrans,[ResizeWidth ResizeHeight],'bilinear');
%         ImgresizeFL1 = imresize(ImgresizeFL1,[ResizeWidth ResizeHeight],'bilinear');
%         ImgresizeFL1Thres = imresize(ImgresizeFL1Thres,[ResizeWidth ResizeHeight],'bilinear');
%         BWfinalresize = logical(imresize(BWfinal,[ResizeWidth ResizeHeight],'nearest'));
%         RGBresize = uint8(zeros(ResizeWidth, ResizeHeight,size(RGB,3)));
%         for layers = 1:size(RGB,3)
%             RGBresize(:,:,layers) = imresize(RGB(:,:,layers),[ResizeWidth ResizeHeight],'nearest');
%         end
%             
%         MaskedFL1Pixel = nonzeros(ImgresizeFL1Thres);
        
        %create masked transmission image
%         ImgresizeTransMasked = ImgresizeTrans .* BWfinalresize;
        
        %save RawImage .mat files
        cd(ImageMatStore);
        RGB = TransImageBkgSub;
        fname = [num2str(k),'-',num2str(ImageList(ImageIdx)),'_Images.mat'];
        parsave(fname,RGB);
        
        %save Mask .mat files
%         cd(ImageMaskStore);
%         fname = [num2str(k),'-',num2str(ImageList(ImageIdx)),'_Masks.mat'];
%         parsave(fname,RGB);
        
        %save image .jpg files
        cd(ImageFileStore);
        TransName = [num2str(k),'-',num2str(ImageList(ImageIdx)),'_TransRaw','.jpg'];
%         FL1Name = [num2str(k),'-',num2str(ImageList(ImageIdx)),'_FL1Raw','.jpg'];
        Trans = mat2gray(TransImageBkgSub);
%         FL1 = mat2gray(ImgresizeFL1);
%         imwrite(uint8(FL1*255),colorMap, FL1Name);
        imwrite(uint8(Trans*255),gray, TransName);
        
        
%         RGB = label2rgb(seg_I-1,[0 0 0;1 0 0],'k');
%         f=figure(13);clf;
% %         set(f,'visible','off');
%         imshow(uint8(floor(mat2gray(ImgresizeTrans)*255)));axis equal
% %         ylabel('scanning [pixel]');xlabel('flow [pixel]');
% %         title(['FL1Thres speed resized No.', num2str(ImageIdx),' [pixel size: 0.25um]']);xlabel('flow pixel');ylabel('scanning pixel');
%         hold on
%         h = imshow(RGBresize);
%         set(h, 'AlphaData', 0.3);
%         F = getframe(gca);
%         
%         TransMaskedName = [num2str(k),'-',num2str(ImageList(ImageIdx)),'_MaskOverlay','.jpg'];
%         imwrite(F.cdata, TransMaskedName);
%         axis off
        cd(storedir);
        % toc
    end
    cd(currentdir);
end






