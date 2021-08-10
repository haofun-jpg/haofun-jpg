%==========================================================================
% �Ŧ��@�פ��R���i�����W�W�е{�� v1.70821 beta
% ���:2018.6.20
%�y���e�z: 
% 1.�W�ЫG���k�@��
% 2.���RPM2.5 PM10
%==========================================================================
clear
%�D���|
main_folder_name = cd;

%% �W�W�Х��n�ƾ�
C = load('../Output data\C.txt');
M = load('../Output data\M.txt');
EV = load('../Output data\EV.txt');
light = load('../Output data\light.txt');
CMF = load('../���л�\CMF.txt');
Ma = [0.40024 0.70760 -0.08081 ;-0.22603 1.16532 0.04570 ; 0  0  0.91822];
White_D65 = load('../Output data\Camera white point.txt');
White_light = load('../Output data\light white point.txt');

day = 1;
sourcePM25 = load('../Air pollution PM2.5 PM10\PM252.txt');
sourcePM25 = sourcePM25(:,day);
sourcePM10 = load('../Air pollution PM2.5 PM10\PM102.txt');
sourcePM10 = sourcePM10(:,day);
PM25 = interp1(1:9,sourcePM25,1:0.18:9,'spline')';
PM10 = interp1(1:9,sourcePM10,1:0.18:9,'spline')';

%Ū���Ϥ�
[FileName,PathName] = uigetfile('*.*','���45�i�ɶ��I�v��','MultiSelect','on');

picture = imread([PathName,char(FileName(:,1))]);
figure(1), imshow(picture);
set(gcf,'keypressfcn','key = get(gcf,''CurrentCharacter'');') %��L�ƥ�
figure_state = figure(1);
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

%% �W�W�мv���޳N
spectral_range = 380-379:780-379;

data = [];
for i = 1:size(FileName,2)
    picture = imread([PathName,char(FileName(:,i))]);

RGB = ones(1,3);
for j = 1:size(position,1)
RGB(j,:) = median(reshape(picture(position(j,2):position(j,2)+position(j,4),position(j,1):position(j,1)+position(j,3),:),[],3),1);
end
% RGB = double(reshape(imresize(picture,0.01),[],3));

XYZ = rgb2xyz(RGB./255)'.*100;
XYZ = inv(Ma) * diag( (Ma*White_light)./(Ma*White_D65) ) *Ma *XYZ;
extend = [  ones(1,size(XYZ,2));... 
            XYZ;... 
            XYZ(1,:).*XYZ(2,:);  XYZ(2,:).* XYZ(3,:);XYZ(1,:).*XYZ(3,:); XYZ.^2;...
            XYZ(1,:).*XYZ(2,:).*XYZ(3,:);XYZ.^3;...
            XYZ(1,:).*XYZ(2,:).^2; XYZ(1,:).*XYZ(3,:).^2;  XYZ(1,:).^2.*XYZ(2,:);... 
            XYZ(2,:).*XYZ(3,:).^2; XYZ(1,:).^2.*XYZ(3,:); XYZ(2,:).^2.*XYZ(3,:)
         ];
CorrectXYZ = (C*extend);
extend = [CorrectXYZ;CorrectXYZ(1,:).*CorrectXYZ(2,:); CorrectXYZ(2,:).*CorrectXYZ(3,:); CorrectXYZ(1,:).*CorrectXYZ(3,:); CorrectXYZ(1,:).*CorrectXYZ(2,:).*CorrectXYZ(3,:)];
data(:,:,i) = EV(spectral_range,:)*M*extend;

end

%% �G���k�@��
data2 = reshape(max(data,[],2),[],45);
[coeff,score,latent,tsquared,explained,mu] = pca(data2','Algorithm','eig','Centered',false); %���R�Ӷ����G�ׯS�x
data3 = data ./ repmat(reshape(score(:,1)./max(score(:,1)),1,1,[]),size(data,1),size(data,2));
data3avg = reshape(mean(data3,1),[],45);
data3_550nm = reshape(data3(550-379,:,:),[],45);
% dataX = data./repmat(coeff(:,1),1,size(position,1));
% dataX2 = mean(dataX./repmat(dataX(1,:),size(FileName,2),1),2);

%% 45�Ӯɶ��I�G���k�@�Ƶ��G
figure('Name','Source'),
for k = 1:45
subplot(5,9,k)
plot(380:780,data(:,:,k))
title(num2str(k))
axis([380 780 0 1])
end

figure('Name','Normalize'),
for k = 1:45
subplot(5,9,k)
plot(380:780,data3(:,:,k))
title(num2str(k))
axis([380 780 0 1])
end

%% �䤤9 �Ӯɶ��I�G���k�@�Ƶ��G(��ܮɶ�)
figure('Name','Source'),
s=[1 8 12 17 23 27 30 43 36];
for k = 1:9
subplot(3,3,k)
plot(380:780,data(:,:,s(1,k)))
title(['Time: ' num2str(k+8) ':00'])
axis([380 780 0 1])
xlabel('Wavelength (nm)')
ylabel('Reflectance (a.u.)')
end

figure('Name','Normalize'),
for k = 1:9
subplot(3,3,k)
plot(380:780,data3(:,:,s(1,k)))
title(['Time: ' num2str(k+8) ':00'])
axis([380 780 0 1])
xlabel('Wavelength (nm)')
ylabel('Reflectance (a.u.)')
end
%% ����PM2.5�@��
M_avg = max(data3avg,[],2);
M_550 = max(data3_550nm,[],2);
a = 0.005;
b = 0.005;
simPM25 = (repmat(log(M_avg+b),1,45)-log(data3avg))./a';

%%
newsimPM25 = [];
for k = 1:size(simPM25,1)
    for j = 1:9
        for i = 1:5
            d = abs(sourcePM25(j,1) - simPM25(k,(j-1)*5+1));
            if abs(sourcePM25(j,1) - simPM25(k,(j-1)*5+i))<= d
            newsimPM25(k,j) = simPM25(k,(j-1)*5+i);
            end
        end
    end
end

figure
%plot(newsimPM25','p')
boxplot([newsimPM25])
hold on
plot(sourcePM25,'o')

%% ����PM10�@��
M_avg = max(data3avg,[],2);
M_550 = max(data3_550nm,[],2);
a = 0.001837;
b = 0.005;

simPM10 = (repmat(log(M_avg+b),1,45)-log(data3avg))./a';

%% 
newsimPM10 = [];
for k = 1:size(simPM10,1)
    for j = 1:9
        for i = 1:5
            d = abs(sourcePM10(j,1) - simPM10(k,(j-1)*5+1));
            if abs(sourcePM10(j,1) - simPM10(k,(j-1)*5+i))<= d
            newsimPM10(k,j) = simPM10(k,(j-1)*5+i);
            end
        end
    end
end
figure
%plot(newsimPM10','p')
boxplot([newsimPM10])
hold on
plot(sourcePM10,'o')