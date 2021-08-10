%==========================================================================
% ���N���v���d�򤧥i�����W�W�е{�� step.1
% ���:2017.12.13
% �y���e�z: ��Ʈw�Хܼv���ؼЪ�
% �y�����z: �i��ROI_SelectObject_save�i��ؼЪ��s��(�Ф�clear�ܼ�)
%==========================================================================
%feature('memstats')%�˵��i�ΰO������O
clear
%�D���|
main_folder_name = cd;
%current_path = pwd;
%% �W�W�Х��n�ƾ�
% load('Output data\VIS Hyperspectral data.mat')
C = load('Output data\C.txt');
M = load('Output data\M.txt');
EV = load('Output data\EV.txt');
light = load('light\SUN.txt');

White_D65 = load('Output data\Camera white point.txt');
White_light = load('Output data\light white point.txt');
CMF = load('���л�\CMF.txt');
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
%% �W�W�мv���޳N
%=========================������B�z���v��==================================
[FileName,PathName] = uigetfile('*.*','��ܼv��');
pic = imread([PathName,FileName]);

% �u�ʭץ�
A_lin = rgb2lin(pic);
B_lin = chromadapt(A_lin,White_D65,'ColorSpace','linear-rgb');
pic = lin2rgb(B_lin);

%=====================�W�W�Ъi�q&�v���j�p�վ�================================
inputdata = inputdlg({'�i���̤p��(nm):','�i���̤j��(nm):','�W�иѪR��(nm):'},'�վ�i���d��',1,{'380','780','1'});
spectral_range =str2num (inputdata{1,1}):str2num (inputdata{3,1}):str2num(inputdata{2,1});

%pic_spectrum = ones(size(pic,1),size(pic,2),401);
pic = imresize(pic,0.2);
 
%=========================�W�W�мv���ظm====================================
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
pic_spectrum = reshape(  (EV(spectral_range-379,:)*M*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), size(spectral_range,2));
%pic_spectrum2 = reshape(reshape(pic_spectrum,[],401).*(ones(size(pic,1)*size(pic,2),1)*light'),size(pic,1),size(pic,2),size(spectral_range,2));
%%  �m��v���A�{

uselight_k = 100/sum(CMF(:,2).*light(:,1));
White_uselight =uselight_k.*(light'*CMF)'./100;

new_picXYZ = reshape(uselight_k.*(reshape(pic_spectrum,[],size(spectral_range,2))*(repmat(light(spectral_range-379,:),1,3).*CMF(spectral_range-379,:))),size(pic,1), size(pic,2), 3);
new_pic = xyz2rgb(reshape( (inv(Ma)*diag( (Ma*White_D65)./(Ma*White_uselight) )*Ma*reshape(new_picXYZ,[],3)')'  ,size(pic,1), size(pic,2), size(pic,3))/100).*255;
new_pic = uint8(new_pic);
figure(1),imshow(new_pic)

clear uselight_k new_picXYZ XYZ XYZ_D65 extend CorrectXYZ  Ma CMF%�M���L���ܼ�
%% �����N�d��
clear position mask select_pic%���s(��K���B��run section)

disp('�i��ܥؼЪ�����:tree�Bclouds�Broof�Btube�Blamp�Belectric_tower�Bship�Bocean�Bairplane�Bcar ');
inputdata = inputdlg({'�Q��諸�ƶq : '},'�ؼЪ�����',[1  40],{'2'});
number = str2num (inputdata{1,1});
%===========================���mask�{��====================================

%�h���d�򤧦�m�Ȧs_______________________________________________________
for z =1:number 
clear  col row BW BW2 position0 %�M���ܼ�

figure(1),imshow(new_pic)
title('�п��: tree �B clouds �B roof �B tube �B lamp �Belectric_tower �Bship �Bocean�Bairplane�Bcar �䤤�@��','Interpreter','none');
items = inputdlg({ ['��',num2str(z),'�ӷQ���ؼЪ��W�� : '] , ' �O�_�������(y/n) :' },'�ؼЪ�����',[1  40]);

if char(items{2,1}) == 'y'
    title('�е���d�����');
    h = imrect;
    position0 = wait(h);
    position0 = fix(position0);
    new_pic1 = new_pic( position0(1,2) : position0(1,2)+position0(1,4) , position0(1,1) : position0(1,1)+position0(1,3),:); 
else
    position0 = zeros(1,2);
    new_pic1 = new_pic; 
end        
disp(['�п�� ',char(items{1,1}),' ���N�d��'])
[BW, col , row] = roipoly(new_pic1);

%Ū���Ұ��ؼЪ��y�ت��C��..................................................
object = items{1,1};
color = struct( 'tree' ,[0,1,0] ,'clouds',[0,0,1],'roof',[1,0,0] ,'tube',[0.8,0.8,0.8] ,'lamp',[0.2 0.2 0.2],'electric_tower',[0,0,0],'ship',[1 1 0],'ocean',[0 1 1],'airplane',[1 0 1],'car',[1 1 1]);
color_object(z,:) = color.(object);

%���ؼЪ�position�Ȧs.....................................................
position(z) =struct('x',[col]+position0(1,1),'y',[row]+position0(1,2));%�C����諸�x�}��m

%���ؼЪ�mask����.........................................................
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
allBW =  uint8(sum(single_BW,3));%���������d��

select_pic = zeros(size(pic,1),size(pic,2),3);%��m�O����Ŷ�
%���N�d���^��_____________________________________________________________
tic %�}�l�p��B��ɶ�
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
time = toc; %�����p��B��ɶ�
%__________________________________________________________________________
close(figure(1))
clear  col row BW BW1 BW2 BW3 i j inputdata new_pic1 color  object %�M���L���ܼ�
%% ���B�z�覡
select_pic = uint8(select_pic);
R1 = menu('�B�z�Ұ��覡:','�y��','�u��ܰ��d��');
figure('Name','���d�򤺼v��'),

switch  R1
  case 1
     %�e��
     imshow(new_pic),title('���^�����d��v��'); 
     hold on
     for k=1:number 
     plot(getfield(position,{k},'x'),getfield(position,{k},'y'),'color',color_object(k,:) ,'LineWidth',2)
    text( position(k).x(2),position(k).y(2),char(select_object{k,1}),'Color','k','HorizontalAlignment','right','FontSize',14,'Interpreter','none'); %�ХܦW��
     end
     
  case 2
    %���~���B��
    imshow(select_pic,[]),title(['�^���v�� �B��time : ',num2str(time),'  seconds.']);
end
%% ���N�d���W���^��
disp('�е���Ū���ɶ�....')  
tic %�}�l�p��B��ɶ�
[all_row , all_col] = find(select_pic(:,:,1) > 0 | select_pic(:,:,2) > 0 | select_pic(:,:,3) > 0);

pic_spectrum2 = zeros(size(spectral_range,2),size(all_row,1));%��m�O����Ŷ�
for i=1:size(all_row,1)
pic_spectrum2(:,i) = reshape( pic_spectrum(all_row(i,1),all_col(i,1),:),[],size(spectral_range,2))';
end


figure('Name','���d���W��'), plot( spectral_range(1:end) , pic_spectrum2 )
axis([ spectral_range(1,1)  spectral_range(1,end) 0 1])
time = toc; %�����p��B��ɶ�
title(['���pixel�ƶq : ',num2str(size(all_row,1)),'  ���I  ' ,'   �B��time : ',num2str(time),' s.']);

clc
clear k i z
%% �s��mask(����)
%{
cd(uigetdir(main_folder_name,'mask�x�s���|'))
string = 'JPG';
change_string = 'tif';
FileName1 = strrep(FileName ,string ,change_string);
imwrite(allBW,['mask_',FileName1])
cd(main_folder_name)
%}