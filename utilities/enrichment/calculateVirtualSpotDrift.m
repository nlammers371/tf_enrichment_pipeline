function driftTol = calculateVirtualSpotDrift(refVecStruct,PixelSize)

  % calculate average frame-over-frame particle drift from data
  lin_diff_vec = diff(refVecStruct.particle_order_ref);
  x_diff_vec = diff(refVecStruct.spot_x_ref);
  y_diff_vec = diff(refVecStruct.spot_y_ref);
  dr_vec = sqrt(x_diff_vec.^2+y_diff_vec.^2);
  dr_vec = dr_vec(lin_diff_vec==0);

  % sets sigma of movement for virtual spot
  driftTol = nanmedian(dr_vec);%*PixelSize;