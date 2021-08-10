%==========================================================================
% �i�����W�W�мv���޳N�Ҳիظm
% ���:2017.12.13
%�y���e�z: �ؼҲ�- �ǥѼv��RGB channels�إ�380~780nm�W�a
%�y�����z: �z�L����&�q���W�Ф���A�ت���XC,M,EV
%==========================================================================
clear
%�D���|
main_folder_name = cd;   % ~~cd:��ܷ�e��Ƨ���m
%% ��ܹ��ɩά�����
msgbox(' "��@����v��" �B "��@�����æ���v��" �ﶵ�n��ӧO24�i��� ','�i������','Help');
R1 = menu('������ɫ��s�ﶵ','��@����v��','��@�����æ���v��','color checker��O�v��','txt�ƾ�');

switch  R1
    
    case 0
    msgbox('�L��ܫ��s', 'Wanring msg' ,'warn')
    clc
    return
    
    case 1
    %��@�����RGB
    [FileName,PathName] = uigetfile('*.*','���24�ئ���䤤���v��','MultiSelect','on');
    
        for i = 1:size(FileName,2)
        
            color_picture = imread([PathName,char(FileName(:,i))]);
            camera_RGB(:,i) = mean(reshape(color_picture,[],3),1); %(�C�i�v�������u�]�t��@���)
        end

    case 2
    %Ū���h������ɧ�RGB (om��L��)
    [FileName,PathName] = uigetfile('*.*','��ܹ��ɤ����h�ئ��','MultiSelect','on');
    
        for i = 1:size(FileName,2)
        
            color_picture = imread([PathName,char(FileName(:,i))]);
            figure(i), imshow(color_picture);
            
            h = imrect;
            position = wait(h);
            position = fix(position);

            camera_RGB(:,i) =  median(reshape(color_picture(position(1,2):position(1,2)+position(1,4),position(1,1):position(1,1)+position(1,3),:),[],3));
            %~~position(1,2)��� �Ĥ@�C�ĤG�檺�� �N�Oy ,�P�zposition(1,2)+position(1,4)�N�Oy+d
            
        end
        
         
    case 3
    %24����O ��24��RGB
    [FileName,PathName] = uigetfile('*.JPG','Ū�����᪺color checker��O');  
    color_picture = imread([PathName,FileName]);
    
    figure, imshow(color_picture);
    camera_RGB = [];
        
        for i = 1:24  %��24��
            h = imrect;
            position = wait(h);
            position = fix(position);
            camera_RGB(:,i) = fix(mean(reshape(color_picture(position(1,2):position(1,2)+position(1,4),position(1,1):position(1,1)+position(1,3),:),[],3)));
        end
        
    case 4
    %Ū���w�׶�������cameraRGB
    [FileName,PathName] = uigetfile('*.txt','Ū��camera RGB txt');     %uigetfile(��������]�w,��ܮت��Х�)
    camera_RGB = load([PathName,FileName]);
end

%% �u�ʦ�A���ץ�

D65_XYZ = whitepoint('d65');

camera_RGB_3D = reshape(uint8(camera_RGB)',24,1,3);
A_lin = rgb2lin(camera_RGB_3D); % ~~ sRGB-> �u��RGB ; B=rgb2lin(A):�M�PA�v��sRGB�Ȥ���gamma�ץ��A��B�v����u��RGB��
B_lin = chromadapt(A_lin,D65_XYZ,'ColorSpace','linear-rgb'); %~~ B = chromadapt(A,illuminant) :�ھڳ��������վ�A����m���šC���������P��J�Ϲ��B��ۦP����m�Ŷ��C
camera_RGB = double(reshape( lin2rgb(B_lin) ,[] ,3)');

%(1)Camera RGB �� XYZ
camera_XYZ = rgb2xyz(camera_RGB'./255)'.*100;
%(2)Camera linRGB �� XYZ
camera_linRGB = double(reshape( B_lin ,[],3))';
camera_linXYZ = rgb2xyz(camera_linRGB'./255 ,'colorSpace','linear-rgb')'.*100;

clear camera_RGB_3D A_lin B_lin
%% ��ܤϮg�W�лP������
R2 = menu('����Ϯg�W�пﶵ','�s�q���Ϯg�W��','�з�color_Rspectrum.txt�ƾ�');

switch  R2
    case 1
    %���e:24������л��ƾ��ഫ���Ϯg�W��
    Data=[0,380:1:780]';
    
         for i = 1:24
         Read= load ([ '���л�\spectrum_data\', num2str(i) ,'.txt']); %Ū��1~24.��l�q���ɮ�
         SPData = [i,interp1(Read(:,1),Read(:,2),380:1:780)]; %����
         Data(:,i+1)=[SPData]';
         end
         
    spectrum_data = Data(2:402 , 2:25);
    spectrum_data = (spectrum_data-min(spectrum_data(:)))/(max(spectrum_data(:))-min(spectrum_data(:)));
    color_Rspectrum = spectrum_data ./ repmat(spectrum_data(:,19) , 1 , 24) ;

    case 2
    %���e:Ū���з�24����Ϯg�W��
    [FileName,PathName] = uigetfile('*.txt','Ū��24������q�W�� txt');  
    color_Rspectrum = load([PathName,FileName]);
end
%%
%==========================================================================
%�Ϯg�W�� �� XYZ
%CMF���t����
%color_Rspectrum����m24�������Ϯg�W��(���t����)
%light����m�����W��
%==========================================================================
CMF = load('���л�\CMF.txt');

%Ū��������
[FileName,PathName] = uigetfile('*.txt','Ū���q�����I����');     
light = load([PathName,FileName]);

k =100/sum(CMF(:,2).*light(:,1));
spectrum_XYZ = k.*((color_Rspectrum.*(repmat(light,1,24)))'*CMF)';
%% �۾���A���ഫ (D65 �� �q�����I����)
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
%D65�W����XYZ�D�k or �зǪ�����[0.95047; 1.00000; 1.08883]
White_D65 = load('light\D65.txt');

k_D65 = 100/sum(CMF(:,2).*White_D65(:,1));
D65_XYZ = k_D65.*(White_D65'*CMF)'./100;    

%'�q�����I����'�W����XYZ�D�k
light_XYZ = k.*(light'*CMF)'./100;
camera_XYZ = inv(Ma) * diag(( Ma*light_XYZ)./(Ma*D65_XYZ)) *Ma *camera_XYZ;

%% �ե��x�}C
extend_camera_XYZ = [  ones(1,24);... %~~�`��
                       camera_XYZ;... %~~�@�� X Y Z
                       camera_XYZ(1,:).*camera_XYZ(2,:);  camera_XYZ(2,:).* camera_XYZ(3,:);camera_XYZ(1,:).* camera_XYZ(3,:);  camera_XYZ.^2;...�@%~~�G�� XY YZ XZ
                       camera_XYZ(1,:).*camera_XYZ(2,:).*camera_XYZ(3,:); camera_XYZ.^3;... �@%~~�T���@XYZ�@X^3 Y^3 Z^3�@�@
                       camera_XYZ(1,:).*camera_XYZ(2,:).^2; camera_XYZ(1,:).*camera_XYZ(3,:).^2;  camera_XYZ(1,:).^2.*camera_XYZ(2,:);... %~~�T�� XY^2 XZ^2 X^2Y
                       camera_XYZ(2,:).*camera_XYZ(3,:).^2; camera_XYZ(1,:).^2.*camera_XYZ(3,:); camera_XYZ(2,:).^2.*camera_XYZ(3,:);... %~~�T�� YZ^2 X^2Z Y^2Z
                    ];

C = (spectrum_XYZ)*pinv(extend_camera_XYZ); %�DC�ե��x�}
correction_XYZ = C * extend_camera_XYZ; %�ocamera�bsun�ե��᪺XYZ

%����ڻ~�t
i=1:24 ;
RMSE_XYZ= sqrt ( sum( (correction_XYZ(:,i) - spectrum_XYZ(:,i)) .^2) / 3 )' ;   
RMSE_allXYZ = sum ( RMSE_XYZ(:)/24);

% Coefficient �����Y�Ƶ���(�i��cftool�u�ʦ^�k�u��)
%(1)�D�u�ʮ�XYZ(�u�ʦ�A���ե��L)
R=corrcoef(spectrum_XYZ(2,19:24),camera_XYZ(2,19:24));
R2_original = R(2)^2;
%(2)�u�ʮ�XYZ
R=corrcoef(spectrum_XYZ(2,19:24),camera_linXYZ(2,19:24));
R2_lin = R(2)^2;

%% PCA���R&�����W��

%PCA���R:�oEV
[COEFF, SCORE, LATENT, TSQUARED, explained, mu]  = pca(color_Rspectrum','Algorithm','eig','Centered',false);
EV = COEFF(:,1:12); %�S�x�V�q
alpha = SCORE(:,1:12)'; %�S�x��

%�D����EXPLAINED
figure('Name',['�D������ : ',num2str(sum(explained(1:12))),'%'],'position', [500 400 1000 500]), %~~ 'position',[left bottom width height]
for j=1:12
   subplot(3,4,j),
   plot( EV(:,j) )
   title(['��',num2str(j),'�D����','(' ,num2str(explained(j,:) ) ,'%)' ])    
end

%�X�i�x�}:�oM
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

%�����W��
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

%����ڻ~�t
i=1:24 ;
RMSE_spectrum = sqrt ( sum( (simulate_spectrum(:,i) - color_Rspectrum(:,i)) .^2) / 401 )' ;   
RMSE_allspectrum= sum ( RMSE_spectrum(:)/24);
%% �W�Ц�t�p��
Lab1=xyz2lab(spectrum_XYZ'/100);
Lab2=xyz2lab(simulate_spectrum_XYZ'/100);

LabRGB1 = uint8(xyz2rgb( (inv(Ma)*diag( (Ma*D65_XYZ)./(Ma*light_XYZ) )*Ma *spectrum_XYZ)' /100).*255);
LabRGB2 = uint8(xyz2rgb( (inv(Ma)*diag( (Ma*D65_XYZ)./(Ma*light_XYZ) )*Ma *simulate_spectrum_XYZ)'/100).*255);


Total_CIE76 = sum((Lab1-Lab2).^2,2).^0.5;
Total_AvgCIE76 = mean(Total_CIE76); %CIE76��t
Total_CIE2000 = deltaE00(Lab1, Lab2)';
Total_AvgCIE2000 = mean(Total_CIE2000);%CIE2000������t

%% �����~�t���G
R3 = menu('�W�п��','�U����W�Ф��','���W�Ф��','��t�ϥܤ��');

switch R3
    
   case 1
   %�U����W�Ф��
   for i=1:24
   figure(i),plot(380:780,simulate_spectrum(:,i),380:780,color_Rspectrum(:,i));
   axis([380 780 0 max(color_Rspectrum(:,19))])
   legend('�����W��','�q���W��')
   end


   case 2
   %������W�Ф��
   figure('position', [100 300 1600 600]),
   
   subplot(1,2,1),plot(380:780,color_Rspectrum);
   axis([380 780 0 max(color_Rspectrum(:))])dsq
   title('�з��W��')
   
   subplot(1,2,2),plot(380:780,simulate_spectrum);
   axis([380 780 0 max(simulate_spectrum(:))])
   title('�����W��')

   
   
   case 3
   %��t�C����(��j��)
   L1=repmat(reshape(LabRGB1,24,1,3),1,5);
   L2=repmat(uint8(ones(24,1,3).*255),1,1);
   L3=repmat(reshape(LabRGB2,24,1,3),1,5);
   hFigure = imtool([L1 L2 L3]); %overview�i��j�Y�p/pix vaule�i��RGB
   set(hFigure,'NumberTitle','Off','Name','��t���(���/����)');

end
   
clear R R1 R2 R3 i j k
%% �s��

cd('Output data')

save('C.txt','C','-ascii')                              %�ե��x�}
save('M.txt','M','-ascii')                              %�ഫ�x�}
save('EV.txt','EV','-ascii')                            %�S�x�V�q
save('light.txt','light','-ascii')                      %����
%save('CameraRGB.txt','CameraRGB','-ascii')              %Camera RGB
save('Camera white point.txt','D65_XYZ','-ascii')     %�۾����I
save('light white point.txt','light_XYZ','-ascii')    %�������I
save('VIS Hyperspectral data.mat','C','M','EV','-v7.3')

cd(main_folder_name)
