% Flux density in the steel at the current peak (open circuit), from ct_oc_<I>A.dat.

ct_setup;
ct_params;
graphics_toolkit('gnuplot');
for I = [300 500]
  load(sprintf('ct_oc_%dA.dat', I));
  Bs = Bmap(:,3)/stack_f;
  fig = figure('visible', 'off', 'position', [100 100 700 560]);
  graphics_toolkit(fig, 'gnuplot');
  hold on;
  edges = 0:0.2:2.2;
  cm = jet(numel(edges));
  for q = 1:numel(edges)-1
    k = Bs >= edges(q) & Bs < edges(q+1);
    if any(k), plot(Bmap(k,1), Bmap(k,2), 's', 'markersize', 4, 'markerfacecolor', cm(q,:), 'color', cm(q,:)); end
  end
  k = Bs >= edges(end);
  if any(k), plot(Bmap(k,1), Bmap(k,2), 's', 'markersize', 4, 'markerfacecolor', [1 0 1], 'color', [1 0 1]); end
  axis equal;  grid on;  xlabel('x [mm]');  ylabel('y [mm]');
  title(sprintf('%g A rms: B in the steel at the current peak, max %.2f T (colour steps of 0.2 T)', I, max(Bs)));
  print(fig, sprintf('fig_bmap_%dA.png', I), '-dpng', '-r100');
end
