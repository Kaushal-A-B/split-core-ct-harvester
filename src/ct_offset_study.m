% Flux linkage with the line conductor moved off-centre in the window, open secondary and
% near short circuit (i2 = -0.9*i1/N2), relative to the centred case.

clear; clc;
ct_setup;
ct_params;
I1pk = I1_design*sqrt(2);
wx = bw/2;  ox = wx + c;

offsets = [0 0; 3 0; -3 0; 0 4; 0 -4; 3 4; -3 -4];
cases = [0, -0.9*I1pk/N2];

openfemm(1);
res = zeros(size(offsets, 1), 4);
for k = 1:size(offsets, 1)
  ct_draw_model(offsets(k,1), offsets(k,2), 0, []);
  mi_saveas(fullfile(modeldir, 'ct_offset.fem'));
  mi_modifycircprop('Primary', 1, I1pk);
  for q = 1:2
    mi_modifycircprop('Secondary', 1, cases(q));
    mi_analyze(1);  mi_loadsolution;
    cp = mo_getcircuitproperties('Secondary');
    pa = mo_getpointvalues(wx, 0);  pb = mo_getpointvalues(ox, 0);
    res(k, 2*q-1) = cp(3);
    res(k, 2*q)   = (pa(1) - pb(1))/(c*1e-3);
    mo_close;
  end
  mi_close;
  fprintf('offset (%+.0f, %+.0f) mm: lambda open %.5f Wb, loaded %.5f Wb\n', offsets(k,:), res(k,1), res(k,3));
end
closefemm;

fprintf('\noffset x, y [mm]   flux linkage open [%%]   flux linkage loaded [%%]   B open [%%]\n');
for k = 1:size(offsets, 1)
  fprintf('   %+3.0f, %+3.0f        %+8.3f             %+8.3f               %+8.3f\n', offsets(k,:), ...
          100*(res(k,1)/res(1,1) - 1), 100*(res(k,3)/res(1,3) - 1), 100*(res(k,2)/res(1,2) - 1));
end
save('-text', 'ct_offset.dat', 'offsets', 'res');
