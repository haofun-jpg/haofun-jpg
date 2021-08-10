%clear
[FileName,PathName] = uigetfile('*.*','選取影像','MultiSelect','on');
input = inputdlg('Delay Time (s)','Delay Time (s)');
Delay = str2num(input{:});

I = imread([PathName,char(FileName(:,1))]);
[X,Map] = gray2ind(I,2^8);
imwrite(X, Map, [PathName,'GIF.gif'], 'GIF', 'WriteMode', 'overwrite', 'DelayTime', Delay, 'LoopCount', Inf );
for i = 2:size(FileName,2)
I = imread([PathName,char(FileName(:,i))]);
[X,Map] = gray2ind(I,2^8);
imwrite(X, Map, [PathName,'GIF.gif'], 'GIF', 'WriteMode', 'append', 'DelayTime', Delay);
end