function driftTol = calculateVirtualSpotDrift(RefStruct,SetID)

  % calculate average frame-over-frame particle drift from data
  set_filter_vec = RefStruct.setID_ref==SetID;
  lin_diff_vec = diff(RefStruct.particle_index_ref(set_filter_vec));  
  
  x_diff_vec = diff(RefStruct.spot_x_ref(set_filter_vec));
  y_diff_vec = diff(RefStruct.spot_y_ref(set_filter_vec));
  dr_vec = sqrt(x_diff_vec.^2+y_diff_vec.^2);
  
  dr_vec = dr_vec(lin_diff_vec==0);

  % sets sigma of movement for virtual spot
  driftTol = nanmedian(dr_vec);%*PixelSize;