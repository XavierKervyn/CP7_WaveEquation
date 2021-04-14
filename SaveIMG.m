function SaveIMG(filename)
    saveas(gcf,filename,'fig')
    saveas(gcf,filename,'jpg')
    %saveas(gcf,filename,'eps')
    print(gcf,filename,'-depsc')
end