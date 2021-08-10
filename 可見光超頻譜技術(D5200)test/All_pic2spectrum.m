
function [test]= funss()

%clear
%主路徑
main_folder_name = cd;

%% 超頻譜必要數據
% load('Output data\VIS Hyperspectral data.mat')
C = coder.load('Output data\C.txt');
M = coder.load('Output data\M.txt');
EV = coder.load('Output data\EV.txt');
light = coder.load('Output data\light.txt');

White_D65 = coder.load('Output data\Camera white point.txt');
White_light = coder.load('Output data\light white point.txt');
CMF = coder.load('spectrometer\CMF.txt');
%光譜儀=spectrometer
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
%% 超頻譜影像技術
%=========================選取欲處理的影像==================================

[FileName,PathName] = uigetfile('*.*','choseimg');
%選擇影像=choseimg
coder.varsize('pic0',[PathName,FileName]);
pic0 = imread([PathName,FileName]);
pic = imresize(pic0,0.2);
pic0 = imresize(pic0,0.2);

% 線性修正
%A_lin = rgb2lin(pic0);
%B_lin = chromadapt(A_lin,White_D65,'ColorSpace','linear-rgb');
%pic = lin2rgb(B_lin);
%figure('Name','線性亮度修正'), imshowpair( pic0 , pic ,'montage' ) ,title(['非線性響應轉換前  ','  線性響應轉換後'],'FontSize' ,14);

%%
%=====================超頻譜波段&影像大小調整================================
inputdata = inputdlg({'wavelengthmin(nm):','wavelengthmax(nm):','Spectrum resolution(nm):'},'Adjusting the wave length range',1,{'380','780','1'});
%波長最小值=wavelengthmin , 波長最大值=wavelengthmax, 頻譜解析度 = Spectrum resolution,
%調整波長範圍 = Adjusting the wavelength range
spectral_range =str2num (inputdata{1,1}):str2num (inputdata{3,1}):str2num(inputdata{2,1});

%pic_spectrum = ones(size(pic,1),size(pic,2),401);
R3 = menu('imgsize','allimg','select range');
%影像大小 = imgsize, 全圖 = allimg,  選取所需範圍 = select range

switch R3
    case 1
        pic = pic;
  
    case 2
        %擷取想的影像範圍
        figure('Name','Fetch org img'), imshow(pic);
        %擷取原始影像 = Fetch org img
    
        select_pic = imcrop(pic);
        pic=select_pic;
        close(figure(1));
        figure('Name','Partial image'), imshow(pic); 
        %部分影像 =  Partial image
   
end

%=============================超頻譜影像建置================================
XYZ_D65 = rgb2xyz(pic).*100;
XYZ = reshape(  (inv(Ma)*diag( (Ma*White_light)./(Ma*White_D65) )*Ma*reshape(XYZ_D65,[],3)')' ,size(pic,1), size(pic,2), 3); %size(pic,1), size(pic,2), 3 = [m n 3] = m,n,3 圖檔組成

extend = cat(3,ones(size(XYZ,1),size(XYZ,2),1),... %常數
    XYZ,... %一階
    XYZ(:,:,1).*XYZ(:,:,2), XYZ(:,:,2).*XYZ(:,:,3), XYZ(:,:,1).*XYZ(:,:,3), XYZ.^2,... %二階
    XYZ(:,:,1).*XYZ(:,:,2).*XYZ(:,:,3), XYZ.^3,... %三階
    XYZ(:,:,1).*  XYZ(:,:,2).^2 , XYZ(:,:,1).*  XYZ(:,:,3).^2 , XYZ(:,:,1).^2.*  XYZ(:,:,2),... 
    XYZ(:,:,2).*  XYZ(:,:,3).^2 , XYZ(:,:,1).^2.*  XYZ(:,:,3) , XYZ(:,:,2).^2.*  XYZ(:,:,3));


CorrectXYZ = reshape( (C*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), 3);

extend = cat(3,CorrectXYZ, CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2), CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3), CorrectXYZ(:,:,1).*CorrectXYZ(:,:,3), CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3));
%{
extend = cat(3,ones(size(CorrectXYZ,1),size(CorrectXYZ,2),1),... %常數
    CorrectXYZ,... %一階
    CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2), CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3), CorrectXYZ(:,:,1).*CorrectXYZ(:,:,3), CorrectXYZ.^2,... %二階
    CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3), CorrectXYZ.^3,... %三階
    CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2).^2 , CorrectXYZ(:,:,1).*CorrectXYZ(:,:,3).^2 , CorrectXYZ(:,:,1).^2.*CorrectXYZ(:,:,2),... 
    CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3).^2 , CorrectXYZ(:,:,1).^2.*CorrectXYZ(:,:,3) , CorrectXYZ(:,:,2).^2.*CorrectXYZ(:,:,3));
%}

%[ 二維模擬頻譜 = 三維擴展矩陣轉二維* M' * EV(所選特定頻率的範圍,:)' ] 再用reshape轉成三維
pic_spectrum = reshape(  (EV(spectral_range-379,:)*M*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), size(spectral_range,2));
pic_spectrum2 = reshape(reshape(pic_spectrum,[],401).*(ones(size(pic,1)*size(pic,2),1)*light'),size(pic,1),size(pic,2),401);

%save('pic_spectrum.mat','pic_spectrum','-v7.3')
%% 影像再現程式
R1 = menu('Img reproduction chose','change light ','no changelight(org light)','no Img reproduction');
%影像再現選擇 = Img reproduction chose 換置光源=change light 
%不換置光源(原始光源) = no changelight(org light) 不進行影像再現=no Img reproduction
%============================光源選擇=======================================
if (R1 == 0) || (R1 == 3)
    return
    
elseif (R1 == 1) %換置光源
    [uselight_FileName,uselight_PathName] = uigetfile('*.*','chose light data');
    %選擇欲換置之光源數據檔 = chose light data
    uselight = load([uselight_PathName,uselight_FileName]);
    
elseif (R1 == 2) %不換置光源
    uselight = light;
    
end
%============================影像呈現=======================================
new_pic = zeros(100000,1);
if (R1 == 1) || (R1 == 2)
    R2 = menu('Img reproduction','all wave/sub wave img reproduction','Individual narrow image reproduction');
    %影像再現選擇 = Img reproduction 全波段/特定波段之影像再現 = all wave/sub wave img reproduction
    %各別窄波段影像再現 = Individual narrow image reproduction
    
    %用於再現影像的光源參數
    uselight_k = 100/sum(CMF(:,2).*uselight(:,1));
    White_uselight =uselight_k.*(uselight'*CMF)'./100;
    
    switch R2
      case 1 %全波段影像再現________________________________________________
        new_picXYZ = reshape(uselight_k.*(reshape(pic_spectrum,[],size(spectral_range,2))*(repmat(uselight(spectral_range-379,:),1,3).*CMF(spectral_range-379,:))),size(pic,1), size(pic,2), 3);
        new_pic = xyz2rgb(reshape( (inv(Ma)*diag( (Ma*White_D65)./(Ma*White_uselight) )*Ma*reshape(new_picXYZ,[],3)')'  ,size(pic,1), size(pic,2), size(pic,3))/100).*255;
        
        figure('Name','Img reproduction result'), imshowpair( pic0 , uint8(new_pic) ,'montage' )
        % 影像再現結果 = Img reproduction result
        title('org img(left)   |   reproduction img(right)','FontSize' ,18)
        %原始影像(左) = org img(left) 再現影像(右) = reproduction img(right)
        
      case 2 %各別窄波段影像再現 (影像將儲存於資料夾中)______________________
        cd(uigetdir(main_folder_name,'save path'))
        %儲存路徑 = save path
        for i = spectral_range
            uselight_spectral_range_k = 100./max(max(CMF(spectral_range-379,:))).*uselight(i-379,1);
            new_picXYZ = reshape(uselight_spectral_range_k.*reshape(pic_spectrum(:,:,spectral_range == i),[],1).*uselight(i-379,:)*CMF(i-379,:),size(pic,1), size(pic,2), 3);
            new_pic = xyz2rgb(reshape( (inv(Ma)*diag( (Ma*White_D65)./(Ma*White_uselight) )*Ma*reshape(new_picXYZ,[],3)')',size(pic,1), size(pic,2), 3)/100).*255;
            
            imwrite(uint8(new_pic),[num2str(i),'nm.tiff'])
        end
        cd(main_folder_name)
   
    end
    return
end
%% 多點頻譜可視畫
clear pic_select_spectrum
clc

fig = figure('Name','chose pixel take Spectrum data','position', [100 50 1700 700]);
%選取 pixel 取得頻譜data = chose pixel take Spectrum data
subplot(1,2,1), imshow( uint8(new_pic));
title(' chose more pixel ','FontSize' ,14)
%可選取多點單pixel = chose more pixel
  
set(gcf,'keypressfcn','key = get(gcf,''CurrentCharacter'');') %鍵盤事件
  figure_state = fig;
  key = 0;
  select_times = 0; %選取次數
  position = zeros(1,4);
  while ishandle(figure_state) == 1 %figure被關閉時中止迴圈
      if (key == 13) %按下Enter，執行1次圈選
        select_times = select_times +1;
        h = imrect;
        position(select_times,:) = fix(getPosition(h));
        key = 0;
      elseif (key == 27) %按下Esc，結束圈選
        break
      else
        pause(0.2) %delay time
      end
  end  
im=1:size(position,1);
pic_select_spectrum=[pic_spectrum(position(im,2),position(im,1),:),[],1,size(position,1)]
for i = 1:size(position,1)
    pic_select_spectrum(:,i) = reshape(pic_spectrum(position(i,2),position(i,1),:),[],1);
end
subplot(1,2,2),plot(380:780,pic_select_spectrum)

title('Spectrum result','FontSize' ,18)
%頻譜結果 = Spectrum result
axis([380 780 0 1])
xlabel({'wavelength(nm)'}),ylabel('Reflectivity(a.u.)');
%axis auto

%% 儲存反射頻譜
choice = questdlg('save Spectrum img data?', 'Menu', 'Yes','No' , 'Yes');
%儲存頻譜影像資料 = save Spectrum img data
switch choice
      case 'Yes'
           disp('save Spectrum data please wait....');
           %存取頻譜檔案中請稍後 = save Spectrum data please wait....
           save('pic_VIS_spectrum.mat','pic_spectrum','-v7.3')
           clc
      case 'No'
           disp('no save');
           %不儲存 = no save
           clc
end      
end