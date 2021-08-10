
cd(uigetdir(main_folder_name,'Àx¦s¸ô®|'))
for i = 1:401
    plot([380+(i-1) 380+(i-1) 380+i 380+i],[0 1 1 0])
    hold on
    plot(380:780,reshape(pic_spectrum(403,244,:),[],1))
    axis([380 780 0 max(pic_spectrum(403,244,:))])
    xlabel('wavelength (nm)')
    ylabel('Intensity (a.u.)')
    saveas(gcf,[num2str(380+i-1),'nm.png'])
    hold off
end