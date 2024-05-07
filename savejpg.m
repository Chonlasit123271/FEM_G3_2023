function savejpg(structure,Quads,State)
%%This function was created to plot in many view point and save jpg
% structure is structural coordinate
% Coord is structural face

% Plot the hexahedral element
views = [0 0; 0 90; 90 0; 45 45]; % Define the view angles for X-Z, Y-Z, and Z-Y planes
view_axis = ["X-Z axis","Y-Z axis","Z-Y axis","X-Y-Z axis"];
for i = 1:4
    figure; % Create a new figure for each plot
    hold on;
    for j = 1:size(Quads, 1)
        face_nodes = structure(Quads(j, :), :);
        patch(face_nodes(:, 1), face_nodes(:, 2), face_nodes(:, 3), 'b');
    end
    xlabel('X-axis: WIDTH(meter)');
    ylabel('Y-axis: DEPTH(meter)');
    zlabel('Z-axis: HEIGHT(meter)');
    title([char(State) 'Hexahedral Element - ' view_axis(i)]);
    grid on;
    view(views(i, :)); % Set the view angle
    hold off;
    exportgraphics(gcf, [char(State) 'Hexahedral_Element_' char(view_axis(i)) '.jpg']); % Save the figure as a JPEG image
end
end