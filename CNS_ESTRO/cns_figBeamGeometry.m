function cns_figBeamGeometry(ct, cst, beamSets, base)
% cns_figBeamGeometry  Beam directions over the anatomy, to check angles BEFORE dose calculation.
%
%   cns_figBeamGeometry(ct, cst, study.beamCandidates.sets, 'D:/.../fig_beamGeometry')
%
% One row per beam set: axial slice and sagittal slice through the PTV centre
% (CT grey), with the outlines of eyes, lenses, optic nerves, chiasm,
% brainstem, cochleae and hippocampi PROJECTED onto the plane (so eyes and
% ears show even when they are not on that slice), and one arrow per beam
% pointing from the source to the target centre.
% Directions use matRad's own rotation (matRad_getRotationMatrix, LPS):
% the axial panel is in radiological view (patient left on the right,
% anterior at the top). The csv lists each beam's unit source direction.
% Patient anatomy: keep the output on hospital storage.

sets = cns_toCell(beamSets);
iT = find(strcmp(cst(:,2), 'PTV'), 1);
if isempty(iT), iT = find(strcmp(cst(:,2), 'CTV'), 1); end
[r, c, s] = ind2sub(ct.cubeDim, cst{iT,4}{1});
cen = [mean(r) mean(c) mean(s)];                 % row (y), col (x), slice (z)
zc = round(cen(3)); xc = round(cen(2));
HU = double(ct.cubeHU{1});
oars = {'Eye_L', 'Eye_R', 'Lens_L', 'Lens_R', 'OpticNerve_L', 'OpticNerve_R', 'Chiasm', 'Brainstem', ...
    'Cochlea_L', 'Cochlea_R', 'Hippocampus_L', 'Hippocampus_R'};
P = cns_okabeIto();
f = cns_journalFigure(17.5, 6.2 * numel(sets) + 1);
csv = {};
for k = 1:numel(sets)
    S = sets{k};
    g = S.gantryAngles(:)'; cpl = S.couchAngles(:)';
    for v = 1:2                                    % 1 axial, 2 sagittal
        ax = subplot(numel(sets), 2, 2 * (k - 1) + v); hold(ax, 'on');
        if v == 1
            img = HU(:, :, zc);
        else
            img = squeeze(HU(:, xc, :))';          % rows = z, cols = y
            img = flipud(img);                     % superior at the top
        end
        imagesc(ax, img, [-200 300]); colormap(ax, gray); axis(ax, 'image', 'off');
        set(ax, 'YDir', 'reverse');                % row 1 at the top (imagesc does not set it under hold on)
        drawOutline(ax, ct, cst, 'PTV', v, zc, xc, [0.85 0 0], '-', 1.0);
        for o = 1:numel(oars)
            drawOutline(ax, ct, cst, oars{o}, v, [], [], [1 1 1] * 0.95, ':', 0.6);
        end
        for b = 1:numel(g)
            R = matRad_getRotationMatrix(g(b), cpl(b));
            d = (R * [0; -1; 0])';                 % unit vector iso -> source, LPS (x, y, z)
            col = P(1 + mod(b, 7) + (b >= 7), :);
            if v == 1
                u = [d(1) / ct.resolution.x, d(2) / ct.resolution.y];   % (col, row) per mm
                p0 = [cen(2) cen(1)];
            else
                u = [d(2) / ct.resolution.y, -d(3) / ct.resolution.z];  % (col = y, row = -z)
                p0 = [cen(1), ct.cubeDim(3) - cen(3) + 1];
            end
            if norm(u) < 1e-6
                plot(ax, p0(1), p0(2), 'o', 'Color', col, 'MarkerSize', 6, 'LineWidth', 1.2);  % beam along the view axis
            else
                L = 90;                             % mm of arrow
                q = p0 + L * u;
                plot(ax, [q(1) p0(1)], [q(2) p0(2)], '-', 'Color', col, 'LineWidth', 1.6);
                plot(ax, p0(1), p0(2), '.', 'Color', col, 'MarkerSize', 10);
                text(ax, q(1), q(2), sprintf(' G%g C%g', g(b), cpl(b)), 'Color', col, 'FontSize', 7, ...
                    'FontWeight', 'bold');
            end
            if v == 1
                csv(end+1, :) = {S.id, g(b), cpl(b), d(1), d(2), d(3)}; %#ok<AGROW>
            end
        end
        if v == 1
            title(ax, sprintf('%s: axial (patient left on the right)', strrep(S.id, '_', ' ')), 'FontWeight', 'normal');
        else
            title(ax, 'sagittal (anterior on the left, superior up)', 'FontWeight', 'normal');
        end
    end
end
cns_saveJournalFigure(f, base, {'set', 'gantry', 'couch', 'srcDirX_left', 'srcDirY_posterior', 'srcDirZ_superior'}, csv);
end

function drawOutline(ax, ct, cst, name, view, zc, xc, col, ls, lw)
i = find(strcmp(cst(:,2), name), 1);
if isempty(i) || isempty(cst{i,4}{1}), return; end
m = false(ct.cubeDim); m(cst{i,4}{1}) = true;
if view == 1
    if isempty(zc), s = any(m, 3); else, s = m(:, :, zc); end
else
    if isempty(xc), s = squeeze(any(m, 2)); else, s = squeeze(m(:, xc, :)); end
    s = flipud(s');
end
if any(s(:))
    contour(ax, double(s), [0.5 0.5], 'LineColor', col, 'LineStyle', ls, 'LineWidth', lw);
end
end
