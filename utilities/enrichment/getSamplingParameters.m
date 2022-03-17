function [min_nucleus_area, max_nucleus_area, snippet_size, minSampleSep,...
          minEdgeSep, roi_rad_spot_pix] = getSamplingParameters(PixelSize,proteinSamplingInfo)
  
  % set min and max acceptable area for nucleus segmentation (NL: need to make this more dynamic to account for different nuclear cycles)   
  min_nucleus_area = round(pi*(2 ./ PixelSize).^2);
  max_nucleus_area = round(pi*(4 ./ PixelSize).^2);

  % set snippet to be 3um in size
  snippet_size = round(proteinSamplingInfo.snippet_size_um ./ PixelSize);

  % set min separation between control and locus 
  minSampleSep = round(proteinSamplingInfo.minSampleSepUm ./ PixelSize);
  minEdgeSep = round(proteinSamplingInfo.minEdgeSepUm ./ PixelSize);

  % calculate ROI size in pixels for spot and control
  roi_rad_spot_pix = round(proteinSamplingInfo.ROIRadiusSpot ./ PixelSize);