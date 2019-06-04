%%%%%%% Nested function to make basic and PBoC-style heatmap plots %%%%%%%
function makeHeatmapPlots(Data, Title, YLabel, FileName,xtick_string,ytick_string,Colormap_heat,write_string,figPath,paperFigPath,PixelSize,lb,ub)
    snippet_fig = figure;
    colormap(Colormap_heat)
    snippet_im = imagesc(Data);
    title(Title)
    snip_size = size(Data,1);
    a = PixelSize;
%     if contains(FileName, 'rel')
%         caxis([relEnrich_lb relEnrich_ub])
%     else
%         caxis([lb ub])
%     end
    caxis([lb ub])
    h = colorbar;
    ylabel('\mum','FontSize',12)
    xlabel('\mum','FontSize',12)
    ylabel(h,YLabel,'FontSize',12)
    eval(xtick_string)
    eval(ytick_string)
    saveas(snippet_fig,[figPath write_string FileName '.png']);

    % make paper fig in PBoC style
    snippet_ax = gca;
    StandardFigurePBoC(snippet_im,snippet_ax);
    saveas(snippet_fig, [paperFigPath write_string FileName '.pdf'])
end