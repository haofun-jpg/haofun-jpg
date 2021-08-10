%==========================================================================
% ���N���v���d��s�ɤ��i�����W�W�е{�� step.2
% ���:2018.2.22
%�y���e�z: �i��each object specturm�s�ɦܸ�Ʈw�������w��Ƨ�
%==========================================================================
inputdata = inputdlg({'�Q�s�ɪ�object specturm : ','����v������ɶ�'},'�ؼЪ��x�s',[1  40],{'�Ҧp:tree�Broof...','�Ҧp:20171214_1214'});
object = char(inputdata{1,1});
time = char(inputdata{2,1});
switch object
        case 'tree'
             label=25;
        case 'clouds'   
             label=50;
        case 'roof'  
             label=75;
        case 'tube' 
             label=100;
        case 'lamp' 
             label=125; 
        case 'electric_tower'
             label=150; 
        case 'ship'
             label=175; 
        case 'ocean'
             label=200; 
        case 'airplane'
             label=225; 
        case 'car'
             label=250; 
    end
%%  �S�wobject mask�^��
for i= 1:size(allBW,1)
    for j=1:size(allBW,2)
        
         if allBW(i,j) == label
             
             for k=1:3
             select_pic(i,j,k)= new_pic(i,j,k);
             end      
             
         else 
             select_pic(i,j,:)=0;
         end
    
    end
end
%% �S�wobject�W���^��

[all_row , all_col] = find(select_pic(:,:,1) > 0 | select_pic(:,:,2) > 0 | select_pic(:,:,3) > 0);

pic_spectrum2 = zeros(size(spectral_range,2),size(all_row,1));%��m�O����Ŷ�
for i=1:size(all_row,1)
pic_spectrum2(:,i) = reshape( pic_spectrum(all_row(i,1),all_col(i,1),:),[],size(spectral_range,2))';
end


figure('Name','���d���W��'), plot( spectral_range(1:end) , pic_spectrum2 )
axis([ spectral_range(1,1)  spectral_range(1,end) 0 1])
title(['   �ؼЪ����� : ',object , '   pixels�ƶq : ',num2str(size(all_row,1)),'  ���I  ' ]);

%% �s��
cd(uigetdir(['..��X����\�v����Ʈwv.2\VIR\',object] ,'object specturm��Ʈw�x�s���|'))
save([time,'_',object,'.mat'],'pic_spectrum2','-v7.3')
cd(main_folder_name)