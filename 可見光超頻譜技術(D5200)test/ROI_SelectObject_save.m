%==========================================================================
% 任意圈選影像範圍存檔之可見光超頻譜程式 step.2
% 日期:2018.2.22
%『內容』: 進行each object specturm存檔至資料庫內的指定資料夾
%==========================================================================
inputdata = inputdlg({'想存檔的object specturm : ','拍攝影像日期時間'},'目標物儲存',[1  40],{'例如:tree、roof...','例如:20171214_1214'});
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
%%  特定object mask擷取
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
%% 特定object頻譜擷取

[all_row , all_col] = find(select_pic(:,:,1) > 0 | select_pic(:,:,2) > 0 | select_pic(:,:,3) > 0);

pic_spectrum2 = zeros(size(spectral_range,2),size(all_row,1));%放置記憶體空間
for i=1:size(all_row,1)
pic_spectrum2(:,i) = reshape( pic_spectrum(all_row(i,1),all_col(i,1),:),[],size(spectral_range,2))';
end


figure('Name','圈選範圍內頻譜'), plot( spectral_range(1:end) , pic_spectrum2 )
axis([ spectral_range(1,1)  spectral_range(1,end) 0 1])
title(['   目標物種類 : ',object , '   pixels數量 : ',num2str(size(all_row,1)),'  個點  ' ]);

%% 存檔
cd(uigetdir(['..整合測試\影像資料庫v.2\VIR\',object] ,'object specturm資料庫儲存路徑'))
save([time,'_',object,'.mat'],'pic_spectrum2','-v7.3')
cd(main_folder_name)