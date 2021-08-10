%==========================================================================
% 可見光超頻譜影像技術模組建置
% 日期:2017.12.13
%『內容』: 建模組- 藉由影像RGB channels建立380~780nm頻帶
%『提醒』: 透過模擬&量測頻譜比較，目的找出C,M,EV
%==========================================================================
clear
%主路徑
main_folder_name = cd;   % ~~cd:顯示當前資料夾位置
%% 選擇圖檔或紀錄檔
msgbox(' "單一色塊影像" 、 "單一不均勻色塊影像" 選項要選個別24張色塊 ','告知說明','Help');
R1 = menu('選取圖檔按鈕選項','單一色塊影像','單一不均勻色塊影像','color checker色板影像','txt數據');

switch  R1
    
    case 0
    msgbox('無選擇按鈕', 'Wanring msg' ,'warn')
    clc
    return
    
    case 1
    %單一色塊找RGB
    [FileName,PathName] = uigetfile('*.*','選擇24種色塊其中之影像','MultiSelect','on');
    
        for i = 1:size(FileName,2)
        
            color_picture = imread([PathName,char(FileName(:,i))]);
            camera_RGB(:,i) = mean(reshape(color_picture,[],3),1); %(每張影像必須只包含單一色塊)
        end

    case 2
    %讀取多色塊圖檔找RGB (om顯微鏡)
    [FileName,PathName] = uigetfile('*.*','選擇圖檔中有多種色塊','MultiSelect','on');
    
        for i = 1:size(FileName,2)
        
            color_picture = imread([PathName,char(FileName(:,i))]);
            figure(i), imshow(color_picture);
            
            h = imrect;
            position = wait(h);
            position = fix(position);

            camera_RGB(:,i) =  median(reshape(color_picture(position(1,2):position(1,2)+position(1,4),position(1,1):position(1,1)+position(1,3),:),[],3));
            %~~position(1,2)表示 第一列第二行的值 就是y ,同理position(1,2)+position(1,4)就是y+d
            
        end
        
         
    case 3
    %24色塊板 取24個RGB
    [FileName,PathName] = uigetfile('*.JPG','讀取拍攝的color checker色板');  
    color_picture = imread([PathName,FileName]);
    
    figure, imshow(color_picture);
    camera_RGB = [];
        
        for i = 1:24  %圈24次
            h = imrect;
            position = wait(h);
            position = fix(position);
            camera_RGB(:,i) = fix(mean(reshape(color_picture(position(1,2):position(1,2)+position(1,4),position(1,1):position(1,1)+position(1,3),:),[],3)));
        end
        
    case 4
    %讀取已匯集完成的cameraRGB
    [FileName,PathName] = uigetfile('*.txt','讀取camera RGB txt');     %uigetfile(文件類型設定,對話框的標示)
    camera_RGB = load([PathName,FileName]);
end

%% 線性色適應修正

D65_XYZ = whitepoint('d65');

camera_RGB_3D = reshape(uint8(camera_RGB)',24,1,3);
A_lin = rgb2lin(camera_RGB_3D); % ~~ sRGB-> 線性RGB ; B=rgb2lin(A):撤銷A影像sRGB值中的gamma修正，使B影像具線性RGB值
B_lin = chromadapt(A_lin,D65_XYZ,'ColorSpace','linear-rgb'); %~~ B = chromadapt(A,illuminant) :根據場景光源調整A的色彩平衡。光源必須與輸入圖像處於相同的色彩空間。
camera_RGB = double(reshape( lin2rgb(B_lin) ,[] ,3)');

%(1)Camera RGB 轉 XYZ
camera_XYZ = rgb2xyz(camera_RGB'./255)'.*100;
%(2)Camera linRGB 轉 XYZ
camera_linRGB = double(reshape( B_lin ,[],3))';
camera_linXYZ = rgb2xyz(camera_linRGB'./255 ,'colorSpace','linear-rgb')'.*100;

clear camera_RGB_3D A_lin B_lin
%% 選擇反射頻譜與光源檔
R2 = menu('選取反射頻譜選項','新量測反射頻譜','標準color_Rspectrum.txt數據');

switch  R2
    case 1
    %內容:24色塊光譜儀數據轉換成反射頻譜
    Data=[0,380:1:780]';
    
         for i = 1:24
         Read= load ([ '光譜儀\spectrum_data\', num2str(i) ,'.txt']); %讀取1~24.原始量測檔案
         SPData = [i,interp1(Read(:,1),Read(:,2),380:1:780)]; %內插
         Data(:,i+1)=[SPData]';
         end
         
    spectrum_data = Data(2:402 , 2:25);
    spectrum_data = (spectrum_data-min(spectrum_data(:)))/(max(spectrum_data(:))-min(spectrum_data(:)));
    color_Rspectrum = spectrum_data ./ repmat(spectrum_data(:,19) , 1 , 24) ;

    case 2
    %內容:讀取標準24色塊反射頻譜
    [FileName,PathName] = uigetfile('*.txt','讀取24色塊測量頻譜 txt');  
    color_Rspectrum = load([PathName,FileName]);
end
%%
%==========================================================================
%反射頻譜 轉 XYZ
%CMF→配色函數
%color_Rspectrum→放置24色塊物體反射頻譜(不含光源)
%light→放置光源頻譜
%==========================================================================
CMF = load('光譜儀\CMF.txt');

%讀取光源檔
[FileName,PathName] = uigetfile('*.txt','讀取量測白點光源');     
light = load([PathName,FileName]);

k =100/sum(CMF(:,2).*light(:,1));
spectrum_XYZ = k.*((color_Rspectrum.*(repmat(light,1,24)))'*CMF)';
%% 相機色適應轉換 (D65 → 量測白點光源)
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
%D65頻譜轉XYZ求法 or 標準直接用[0.95047; 1.00000; 1.08883]
White_D65 = load('light\D65.txt');

k_D65 = 100/sum(CMF(:,2).*White_D65(:,1));
D65_XYZ = k_D65.*(White_D65'*CMF)'./100;    

%'量測白點光源'頻譜轉XYZ求法
light_XYZ = k.*(light'*CMF)'./100;
camera_XYZ = inv(Ma) * diag(( Ma*light_XYZ)./(Ma*D65_XYZ)) *Ma *camera_XYZ;

%% 校正矩陣C
extend_camera_XYZ = [  ones(1,24);... %~~常數
                       camera_XYZ;... %~~一階 X Y Z
                       camera_XYZ(1,:).*camera_XYZ(2,:);  camera_XYZ(2,:).* camera_XYZ(3,:);camera_XYZ(1,:).* camera_XYZ(3,:);  camera_XYZ.^2;...　%~~二階 XY YZ XZ
                       camera_XYZ(1,:).*camera_XYZ(2,:).*camera_XYZ(3,:); camera_XYZ.^3;... 　%~~三階　XYZ　X^3 Y^3 Z^3　　
                       camera_XYZ(1,:).*camera_XYZ(2,:).^2; camera_XYZ(1,:).*camera_XYZ(3,:).^2;  camera_XYZ(1,:).^2.*camera_XYZ(2,:);... %~~三階 XY^2 XZ^2 X^2Y
                       camera_XYZ(2,:).*camera_XYZ(3,:).^2; camera_XYZ(1,:).^2.*camera_XYZ(3,:); camera_XYZ(2,:).^2.*camera_XYZ(3,:);... %~~三階 YZ^2 X^2Z Y^2Z
                    ];

C = (spectrum_XYZ)*pinv(extend_camera_XYZ); %求C校正矩陣
correction_XYZ = C * extend_camera_XYZ; %得camera在sun校正後的XYZ

%均方根誤差
i=1:24 ;
RMSE_XYZ= sqrt ( sum( (correction_XYZ(:,i) - spectrum_XYZ(:,i)) .^2) / 3 )' ;   
RMSE_allXYZ = sum ( RMSE_XYZ(:)/24);

% Coefficient 相關係數評估(可用cftool線性回歸工具)
%(1)非線性時XYZ(線性色適應校正過)
R=corrcoef(spectrum_XYZ(2,19:24),camera_XYZ(2,19:24));
R2_original = R(2)^2;
%(2)線性時XYZ
R=corrcoef(spectrum_XYZ(2,19:24),camera_linXYZ(2,19:24));
R2_lin = R(2)^2;

%% PCA分析&模擬頻譜

%PCA分析:得EV
[COEFF, SCORE, LATENT, TSQUARED, explained, mu]  = pca(color_Rspectrum','Algorithm','eig','Centered',false);
EV = COEFF(:,1:12); %特徵向量
alpha = SCORE(:,1:12)'; %特徵值

%主成分EXPLAINED
figure('Name',['主成分比重 : ',num2str(sum(explained(1:12))),'%'],'position', [500 400 1000 500]), %~~ 'position',[left bottom width height]
for j=1:12
   subplot(3,4,j),
   plot( EV(:,j) )
   title(['第',num2str(j),'主成分','(' ,num2str(explained(j,:) ) ,'%)' ])    
end

%擴展矩陣:得M
extend = [spectrum_XYZ;spectrum_XYZ(1,:).*spectrum_XYZ(2,:); spectrum_XYZ(2,:).*spectrum_XYZ(3,:);spectrum_XYZ(1,:).*spectrum_XYZ(3,:); 
          spectrum_XYZ(1,:).*spectrum_XYZ(2,:).*spectrum_XYZ(3,:) 
         ];
%{
extend = [  ones(1,24);... 
            spectrum_XYZ;... 
            spectrum_XYZ(1,:).*spectrum_XYZ(2,:); spectrum_XYZ(2,:).* spectrum_XYZ(3,:);spectrum_XYZ(1,:).* spectrum_XYZ(3,:); spectrum_XYZ.^2;...
            spectrum_XYZ(1,:).*spectrum_XYZ(2,:).*spectrum_XYZ(3,:); spectrum_XYZ.^3;...
            spectrum_XYZ(1,:).*spectrum_XYZ(2,:).^2; spectrum_XYZ(1,:).*spectrum_XYZ(3,:).^2;  spectrum_XYZ(1,:).^2.*spectrum_XYZ(2,:);... 
            spectrum_XYZ(2,:).*spectrum_XYZ(3,:).^2; spectrum_XYZ(1,:).^2.*spectrum_XYZ(3,:);spectrum_XYZ(2,:).^2.*spectrum_XYZ(3,:);...
         ];
%}
M = alpha*pinv(extend);

%模擬頻譜
extend_correction_XYZ = [correction_XYZ;correction_XYZ(1,:).*correction_XYZ(2,:); correction_XYZ(2,:).*correction_XYZ(3,:); correction_XYZ(1,:).*correction_XYZ(3,:);
                         correction_XYZ(1,:).*correction_XYZ(2,:).*correction_XYZ(3,:) 
                        ];
%{
extend_correction_XYZ  = [  ones(1,24);... 
            correction_XYZ;... 
            correction_XYZ(1,:).*correction_XYZ(2,:); correction_XYZ(2,:).* correction_XYZ(3,:);correction_XYZ(1,:).* correction_XYZ(3,:);  correction_XYZ.^2;...
            correction_XYZ(1,:).*correction_XYZ(2,:).*correction_XYZ(3,:); correction_XYZ.^3;...
            correction_XYZ(1,:).*correction_XYZ(2,:).^2; correction_XYZ(1,:).*correction_XYZ(3,:).^2;  correction_XYZ(1,:).^2.*correction_XYZ(2,:);... 
            correction_XYZ(2,:).*correction_XYZ(3,:).^2; correction_XYZ(1,:).^2.*correction_XYZ(3,:); correction_XYZ(2,:).^2.*correction_XYZ(3,:);...
         ];
 %}              
simulate_spectrum = EV*M*extend_correction_XYZ; 

simulate_spectrum_XYZ= k.*((simulate_spectrum.*(repmat(light,1,24)))'*CMF)';

%均方根誤差
i=1:24 ;
RMSE_spectrum = sqrt ( sum( (simulate_spectrum(:,i) - color_Rspectrum(:,i)) .^2) / 401 )' ;   
RMSE_allspectrum= sum ( RMSE_spectrum(:)/24);
%% 頻譜色差計算
Lab1=xyz2lab(spectrum_XYZ'/100);
Lab2=xyz2lab(simulate_spectrum_XYZ'/100);

LabRGB1 = uint8(xyz2rgb( (inv(Ma)*diag( (Ma*D65_XYZ)./(Ma*light_XYZ) )*Ma *spectrum_XYZ)' /100).*255);
LabRGB2 = uint8(xyz2rgb( (inv(Ma)*diag( (Ma*D65_XYZ)./(Ma*light_XYZ) )*Ma *simulate_spectrum_XYZ)'/100).*255);


Total_CIE76 = sum((Lab1-Lab2).^2,2).^0.5;
Total_AvgCIE76 = mean(Total_CIE76); %CIE76色差
Total_CIE2000 = deltaE00(Lab1, Lab2)';
Total_AvgCIE2000 = mean(Total_CIE2000);%CIE2000平均色差

%% 模擬誤差結果
R3 = menu('頻譜選取','各色塊頻譜比較','全頻譜比較','色差圖示比較');

switch R3
    
   case 1
   %各色塊頻譜比較
   for i=1:24
   figure(i),plot(380:780,simulate_spectrum(:,i),380:780,color_Rspectrum(:,i));
   axis([380 780 0 max(color_Rspectrum(:,19))])
   legend('模擬頻譜','量測頻譜')
   end


   case 2
   %全色塊頻譜比較
   figure('position', [100 300 1600 600]),
   
   subplot(1,2,1),plot(380:780,color_Rspectrum);
   axis([380 780 0 max(color_Rspectrum(:))])dsq
   title('標準頻譜')
   
   subplot(1,2,2),plot(380:780,simulate_spectrum);
   axis([380 780 0 max(simulate_spectrum(:))])
   title('模擬頻譜')

   
   
   case 3
   %色差顏色比較(放大看)
   L1=repmat(reshape(LabRGB1,24,1,3),1,5);
   L2=repmat(uint8(ones(24,1,3).*255),1,1);
   L3=repmat(reshape(LabRGB2,24,1,3),1,5);
   hFigure = imtool([L1 L2 L3]); %overview可放大縮小/pix vaule可看RGB
   set(hFigure,'NumberTitle','Off','Name','色差比較(實際/模擬)');

end
   
clear R R1 R2 R3 i j k
%% 存檔

cd('Output data')

save('C.txt','C','-ascii')                              %校正矩陣
save('M.txt','M','-ascii')                              %轉換矩陣
save('EV.txt','EV','-ascii')                            %特徵向量
save('light.txt','light','-ascii')                      %光源
%save('CameraRGB.txt','CameraRGB','-ascii')              %Camera RGB
save('Camera white point.txt','D65_XYZ','-ascii')     %相機白點
save('light white point.txt','light_XYZ','-ascii')    %光源白點
save('VIS Hyperspectral data.mat','C','M','EV','-v7.3')

cd(main_folder_name)
