%==========================================================================
% 任意圈選影像範圍之可見光超頻譜程式 step.1
% 日期:2017.12.13
% 『內容』: 資料庫標示影像目標物
% 『提醒』: 可至ROI_SelectObject_save進行目標物存檔(請勿clear變數)
%==========================================================================
%feature('memstats')%檢視可用記憶體指令
clear
%主路徑
main_folder_name = cd;
%current_path = pwd;
%% 超頻譜必要數據
% load('Output data\VIS Hyperspectral data.mat')
C = load('Output data\C.txt');
M = load('Output data\M.txt');
EV = load('Output data\EV.txt');
light = load('light\SUN.txt');

White_D65 = load('Output data\Camera white point.txt');
White_light = load('Output data\light white point.txt');
CMF = load('光譜儀\CMF.txt');
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
%% 超頻譜影像技術
%=========================選取欲處理的影像==================================
[FileName,PathName] = uigetfile('*.*','選擇影像');
pic = imread([PathName,FileName]);

% 線性修正
A_lin = rgb2lin(pic);
B_lin = chromadapt(A_lin,White_D65,'ColorSpace','linear-rgb');
pic = lin2rgb(B_lin);

%=====================超頻譜波段&影像大小調整================================
inputdata = inputdlg({'波長最小值(nm):','波長最大值(nm):','頻譜解析度(nm):'},'調整波長範圍',1,{'380','780','1'});
spectral_range =str2num (inputdata{1,1}):str2num (inputdata{3,1}):str2num(inputdata{2,1});

%pic_spectrum = ones(size(pic,1),size(pic,2),401);
pic = imresize(pic,0.2);
 
%=========================超頻譜影像建置====================================
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
pic_spectrum = reshape(  (EV(spectral_range-379,:)*M*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), size(spectral_range,2));
%pic_spectrum2 = reshape(reshape(pic_spectrum,[],401).*(ones(size(pic,1)*size(pic,2),1)*light'),size(pic,1),size(pic,2),size(spectral_range,2));
%%  彩色影像再現

uselight_k = 100/sum(CMF(:,2).*light(:,1));
White_uselight =uselight_k.*(light'*CMF)'./100;

new_picXYZ = reshape(uselight_k.*(reshape(pic_spectrum,[],size(spectral_range,2))*(repmat(light(spectral_range-379,:),1,3).*CMF(spectral_range-379,:))),size(pic,1), size(pic,2), 3);
new_pic = xyz2rgb(reshape( (inv(Ma)*diag( (Ma*White_D65)./(Ma*White_uselight) )*Ma*reshape(new_picXYZ,[],3)')'  ,size(pic,1), size(pic,2), size(pic,3))/100).*255;
new_pic = uint8(new_pic);
figure(1),imshow(new_pic)

clear uselight_k new_picXYZ XYZ XYZ_D65 extend CorrectXYZ  Ma CMF%清除無用變數
%% 圈選任意範圍
clear position mask select_pic%重製(方便此處按run section)

disp('可圈擇目標物種類:tree、clouds、roof、tube、lamp、electric_tower、ship、ocean、airplane、car ');
inputdata = inputdlg({'想圈選的數量 : '},'目標物種類',[1  40],{'2'});
number = str2num (inputdata{1,1});
%===========================圈選mask程式====================================

%多圈選範圍之位置暫存_______________________________________________________
for z =1:number 
clear  col row BW BW2 position0 %清除變數

figure(1),imshow(new_pic)
title('請選擇: tree 、 clouds 、 roof 、 tube 、 lamp 、electric_tower 、ship 、ocean、airplane、car 其中一種','Interpreter','none');
items = inputdlg({ ['第',num2str(z),'個想圈選目標物名稱 : '] , ' 是否裁減邊框(y/n) :' },'目標物種類',[1  40]);

if char(items{2,1}) == 'y'
    title('請裁減範圍邊框');
    h = imrect;
    position0 = wait(h);
    position0 = fix(position0);
    new_pic1 = new_pic( position0(1,2) : position0(1,2)+position0(1,4) , position0(1,1) : position0(1,1)+position0(1,3),:); 
else
    position0 = zeros(1,2);
    new_pic1 = new_pic; 
end        
disp(['請選取 ',char(items{1,1}),' 任意範圍'])
[BW, col , row] = roipoly(new_pic1);

%讀取所圈選目標物描框的顏色..................................................
object = items{1,1};
color = struct( 'tree' ,[0,1,0] ,'clouds',[0,0,1],'roof',[1,0,0] ,'tube',[0.8,0.8,0.8] ,'lamp',[0.2 0.2 0.2],'electric_tower',[0,0,0],'ship',[1 1 0],'ocean',[0 1 1],'airplane',[1 0 1],'car',[1 1 1]);
color_object(z,:) = color.(object);

%圈選目標物position暫存.....................................................
position(z) =struct('x',[col]+position0(1,1),'y',[row]+position0(1,2));%每次圈選的矩陣位置

%圈選目標物mask分類.........................................................
BW2 = uint8(BW);
    switch object
        case 'tree'
             BW2(find(BW2>0))=25;
        case 'clouds'   
             BW2(find(BW2>0))=50;
        case 'roof'  
             BW2(find(BW2>0))=75;
        case 'tube' 
             BW2(find(BW2>0))=100;
        case 'lamp' 
             BW2(find(BW2>0))=125; 
        case 'electric_tower'
             BW2(find(BW2>0))=150; 
        case 'ship'
             BW2(find(BW2>0))=175; 
        case 'ocean'
             BW2(find(BW2>0))=200; 
        case 'airplane'
             BW2(find(BW2>0))=225; 
        case 'car'
             BW2(find(BW2>0))=250; 
    end
if char(items{2,1}) == 'y'    
    BW3 = uint8(zeros(size(pic,1),size(pic,2)));
    BW3(position0(1,2) : position0(1,2)+position0(1,4) , position0(1,1) : position0(1,1)+position0(1,3))=BW3(position0(1,2) : position0(1,2)+position0(1,4), position0(1,1) : position0(1,1)+position0(1,3))+BW2;
else
    BW3 = BW2;
end
single_BW(:,:,z) = BW3;
select_object{z,1} = object;
end
%__________________________________________________________________________
allBW =  uint8(sum(single_BW,3));%全部的圈選範圍

select_pic = zeros(size(pic,1),size(pic,2),3);%放置記憶體空間
%任意範圍內擷取_____________________________________________________________
tic %開始計算運算時間
for i= 1:size(allBW,1)
    for j=1:size(allBW,2)
        
         if allBW(i,j)~=0
             
             for k=1:3
             select_pic(i,j,k)= new_pic(i,j,k);
             end      
             
         else 
             select_pic(i,j,:)=0;
         end
    
    end
end
time = toc; %結束計算運算時間
%__________________________________________________________________________
close(figure(1))
clear  col row BW BW1 BW2 BW3 i j inputdata new_pic1 color  object %清除無用變數
%% 圈選處理方式
select_pic = uint8(select_pic);
R1 = menu('處理所圈選方式:','描框','只顯示圈選範圍');
figure('Name','圈選範圍內影像'),

switch  R1
  case 1
     %畫框
     imshow(new_pic),title('所擷取的範圍影像'); 
     hold on
     for k=1:number 
     plot(getfield(position,{k},'x'),getfield(position,{k},'y'),'color',color_object(k,:) ,'LineWidth',2)
    text( position(k).x(2),position(k).y(2),char(select_object{k,1}),'Color','k','HorizontalAlignment','right','FontSize',14,'Interpreter','none'); %標示名稱
     end
     
  case 2
    %圈選外的遮蔽
    imshow(select_pic,[]),title(['擷取影像 運算time : ',num2str(time),'  seconds.']);
end
%% 任意範圍內頻譜擷取
disp('請等待讀取時間....')  
tic %開始計算運算時間
[all_row , all_col] = find(select_pic(:,:,1) > 0 | select_pic(:,:,2) > 0 | select_pic(:,:,3) > 0);

pic_spectrum2 = zeros(size(spectral_range,2),size(all_row,1));%放置記憶體空間
for i=1:size(all_row,1)
pic_spectrum2(:,i) = reshape( pic_spectrum(all_row(i,1),all_col(i,1),:),[],size(spectral_range,2))';
end


figure('Name','圈選範圍內頻譜'), plot( spectral_range(1:end) , pic_spectrum2 )
axis([ spectral_range(1,1)  spectral_range(1,end) 0 1])
time = toc; %結束計算運算時間
title(['圈選pixel數量 : ',num2str(size(all_row,1)),'  個點  ' ,'   運算time : ',num2str(time),' s.']);

clc
clear k i z
%% 存檔mask(不用)
%{
cd(uigetdir(main_folder_name,'mask儲存路徑'))
string = 'JPG';
change_string = 'tif';
FileName1 = strrep(FileName ,string ,change_string);
imwrite(allBW,['mask_',FileName1])
cd(main_folder_name)
%}