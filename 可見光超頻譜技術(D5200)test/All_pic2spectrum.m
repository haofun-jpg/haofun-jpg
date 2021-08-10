
function [test]= funss()

%clear
%�D���|
main_folder_name = cd;

%% �W�W�Х��n�ƾ�
% load('Output data\VIS Hyperspectral data.mat')
C = coder.load('Output data\C.txt');
M = coder.load('Output data\M.txt');
EV = coder.load('Output data\EV.txt');
light = coder.load('Output data\light.txt');

White_D65 = coder.load('Output data\Camera white point.txt');
White_light = coder.load('Output data\light white point.txt');
CMF = coder.load('spectrometer\CMF.txt');
%���л�=spectrometer
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
%% �W�W�мv���޳N
%=========================������B�z���v��==================================

[FileName,PathName] = uigetfile('*.*','choseimg');
%��ܼv��=choseimg
coder.varsize('pic0',[PathName,FileName]);
pic0 = imread([PathName,FileName]);
pic = imresize(pic0,0.2);
pic0 = imresize(pic0,0.2);

% �u�ʭץ�
%A_lin = rgb2lin(pic0);
%B_lin = chromadapt(A_lin,White_D65,'ColorSpace','linear-rgb');
%pic = lin2rgb(B_lin);
%figure('Name','�u�ʫG�׭ץ�'), imshowpair( pic0 , pic ,'montage' ) ,title(['�D�u���T���ഫ�e  ','  �u���T���ഫ��'],'FontSize' ,14);

%%
%=====================�W�W�Ъi�q&�v���j�p�վ�================================
inputdata = inputdlg({'wavelengthmin(nm):','wavelengthmax(nm):','Spectrum resolution(nm):'},'Adjusting the wave length range',1,{'380','780','1'});
%�i���̤p��=wavelengthmin , �i���̤j��=wavelengthmax, �W�иѪR�� = Spectrum resolution,
%�վ�i���d�� = Adjusting the wavelength range
spectral_range =str2num (inputdata{1,1}):str2num (inputdata{3,1}):str2num(inputdata{2,1});

%pic_spectrum = ones(size(pic,1),size(pic,2),401);
R3 = menu('imgsize','allimg','select range');
%�v���j�p = imgsize, ���� = allimg,  ����һݽd�� = select range

switch R3
    case 1
        pic = pic;
  
    case 2
        %�^���Q���v���d��
        figure('Name','Fetch org img'), imshow(pic);
        %�^����l�v�� = Fetch org img
    
        select_pic = imcrop(pic);
        pic=select_pic;
        close(figure(1));
        figure('Name','Partial image'), imshow(pic); 
        %�����v�� =  Partial image
   
end

%=============================�W�W�мv���ظm================================
XYZ_D65 = rgb2xyz(pic).*100;
XYZ = reshape(  (inv(Ma)*diag( (Ma*White_light)./(Ma*White_D65) )*Ma*reshape(XYZ_D65,[],3)')' ,size(pic,1), size(pic,2), 3); %size(pic,1), size(pic,2), 3 = [m n 3] = m,n,3 ���ɲզ�

extend = cat(3,ones(size(XYZ,1),size(XYZ,2),1),... %�`��
    XYZ,... %�@��
    XYZ(:,:,1).*XYZ(:,:,2), XYZ(:,:,2).*XYZ(:,:,3), XYZ(:,:,1).*XYZ(:,:,3), XYZ.^2,... %�G��
    XYZ(:,:,1).*XYZ(:,:,2).*XYZ(:,:,3), XYZ.^3,... %�T��
    XYZ(:,:,1).*  XYZ(:,:,2).^2 , XYZ(:,:,1).*  XYZ(:,:,3).^2 , XYZ(:,:,1).^2.*  XYZ(:,:,2),... 
    XYZ(:,:,2).*  XYZ(:,:,3).^2 , XYZ(:,:,1).^2.*  XYZ(:,:,3) , XYZ(:,:,2).^2.*  XYZ(:,:,3));


CorrectXYZ = reshape( (C*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), 3);

extend = cat(3,CorrectXYZ, CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2), CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3), CorrectXYZ(:,:,1).*CorrectXYZ(:,:,3), CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3));
%{
extend = cat(3,ones(size(CorrectXYZ,1),size(CorrectXYZ,2),1),... %�`��
    CorrectXYZ,... %�@��
    CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2), CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3), CorrectXYZ(:,:,1).*CorrectXYZ(:,:,3), CorrectXYZ.^2,... %�G��
    CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3), CorrectXYZ.^3,... %�T��
    CorrectXYZ(:,:,1).*CorrectXYZ(:,:,2).^2 , CorrectXYZ(:,:,1).*CorrectXYZ(:,:,3).^2 , CorrectXYZ(:,:,1).^2.*CorrectXYZ(:,:,2),... 
    CorrectXYZ(:,:,2).*CorrectXYZ(:,:,3).^2 , CorrectXYZ(:,:,1).^2.*CorrectXYZ(:,:,3) , CorrectXYZ(:,:,2).^2.*CorrectXYZ(:,:,3));
%}

%[ �G�������W�� = �T���X�i�x�}��G��* M' * EV(�ҿ�S�w�W�v���d��,:)' ] �A��reshape�ন�T��
pic_spectrum = reshape(  (EV(spectral_range-379,:)*M*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), size(spectral_range,2));
pic_spectrum2 = reshape(reshape(pic_spectrum,[],401).*(ones(size(pic,1)*size(pic,2),1)*light'),size(pic,1),size(pic,2),401);

%save('pic_spectrum.mat','pic_spectrum','-v7.3')
%% �v���A�{�{��
R1 = menu('Img reproduction chose','change light ','no changelight(org light)','no Img reproduction');
%�v���A�{��� = Img reproduction chose ���m����=change light 
%�����m����(��l����) = no changelight(org light) ���i��v���A�{=no Img reproduction
%============================�������=======================================
if (R1 == 0) || (R1 == 3)
    return
    
elseif (R1 == 1) %���m����
    [uselight_FileName,uselight_PathName] = uigetfile('*.*','chose light data');
    %��ܱ����m�������ƾ��� = chose light data
    uselight = load([uselight_PathName,uselight_FileName]);
    
elseif (R1 == 2) %�����m����
    uselight = light;
    
end
%============================�v���e�{=======================================
new_pic = zeros(100000,1);
if (R1 == 1) || (R1 == 2)
    R2 = menu('Img reproduction','all wave/sub wave img reproduction','Individual narrow image reproduction');
    %�v���A�{��� = Img reproduction ���i�q/�S�w�i�q���v���A�{ = all wave/sub wave img reproduction
    %�U�O���i�q�v���A�{ = Individual narrow image reproduction
    
    %�Ω�A�{�v���������Ѽ�
    uselight_k = 100/sum(CMF(:,2).*uselight(:,1));
    White_uselight =uselight_k.*(uselight'*CMF)'./100;
    
    switch R2
      case 1 %���i�q�v���A�{________________________________________________
        new_picXYZ = reshape(uselight_k.*(reshape(pic_spectrum,[],size(spectral_range,2))*(repmat(uselight(spectral_range-379,:),1,3).*CMF(spectral_range-379,:))),size(pic,1), size(pic,2), 3);
        new_pic = xyz2rgb(reshape( (inv(Ma)*diag( (Ma*White_D65)./(Ma*White_uselight) )*Ma*reshape(new_picXYZ,[],3)')'  ,size(pic,1), size(pic,2), size(pic,3))/100).*255;
        
        figure('Name','Img reproduction result'), imshowpair( pic0 , uint8(new_pic) ,'montage' )
        % �v���A�{���G = Img reproduction result
        title('org img(left)   |   reproduction img(right)','FontSize' ,18)
        %��l�v��(��) = org img(left) �A�{�v��(�k) = reproduction img(right)
        
      case 2 %�U�O���i�q�v���A�{ (�v���N�x�s���Ƨ���)______________________
        cd(uigetdir(main_folder_name,'save path'))
        %�x�s���| = save path
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
%% �h�I�W�Хi���e
clear pic_select_spectrum
clc

fig = figure('Name','chose pixel take Spectrum data','position', [100 50 1700 700]);
%��� pixel ���o�W��data = chose pixel take Spectrum data
subplot(1,2,1), imshow( uint8(new_pic));
title(' chose more pixel ','FontSize' ,14)
%�i����h�I��pixel = chose more pixel
  
set(gcf,'keypressfcn','key = get(gcf,''CurrentCharacter'');') %��L�ƥ�
  figure_state = fig;
  key = 0;
  select_times = 0; %�������
  position = zeros(1,4);
  while ishandle(figure_state) == 1 %figure�Q�����ɤ���j��
      if (key == 13) %���UEnter�A����1�����
        select_times = select_times +1;
        h = imrect;
        position(select_times,:) = fix(getPosition(h));
        key = 0;
      elseif (key == 27) %���UEsc�A�������
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
%�W�е��G = Spectrum result
axis([380 780 0 1])
xlabel({'wavelength(nm)'}),ylabel('Reflectivity(a.u.)');
%axis auto

%% �x�s�Ϯg�W��
choice = questdlg('save Spectrum img data?', 'Menu', 'Yes','No' , 'Yes');
%�x�s�W�мv����� = save Spectrum img data
switch choice
      case 'Yes'
           disp('save Spectrum data please wait....');
           %�s���W���ɮפ��еy�� = save Spectrum data please wait....
           save('pic_VIS_spectrum.mat','pic_spectrum','-v7.3')
           clc
      case 'No'
           disp('no save');
           %���x�s = no save
           clc
end      
end