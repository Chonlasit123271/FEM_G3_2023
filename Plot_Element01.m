function Plot_Element01(DefE, NodeCoord)
% This function plots mesh with index of nodes and elements
% -----------------------------------------------------------------------
% Meaning of variables
% DefE: matrix defines elements
% NodeCoord: matrix contains coordinates of nodes
% -----------------------------------------------------------------------

% Calculate the minimum distance between nodes
min_distance = min(pdist(NodeCoord));

% Set the limit for plotting elements dynamically
limit = min_distance*1.1;

figure(1);
Plot_Node(NodeCoord);
views = [45 45]; % Define the view angles for X-Z, Y-Z, and Z-Y planes
for vp = 1:1
% Plot horizontal and vertical members
    for i = 1:size(DefE, 1)
        for j = 1:8
            node1 = DefE(i, j);
            if node1 ~= 0 % If the node index is valid
                for k = (j + 1):8
                    node2 = DefE(i, k);
                    if node2 ~= 0 % If the node index is valid
                        % Check if the nodes are connected horizontally or vertically
                        if (NodeCoord(node1, 1) == NodeCoord(node2, 1) && NodeCoord(node1, 2) == NodeCoord(node2, 2)) || ...
                           (NodeCoord(node1, 1) == NodeCoord(node2, 1) && NodeCoord(node1, 3) == NodeCoord(node2, 3)) || ...
                           (NodeCoord(node1, 2) == NodeCoord(node2, 2) && NodeCoord(node1, 3) == NodeCoord(node2, 3))
                            % Ensure that inclined members are not plotted
                            if NodeCoord(node1, 1) == NodeCoord(node2, 1) || NodeCoord(node1, 2) == NodeCoord(node2, 2) || NodeCoord(node1, 3) == NodeCoord(node2, 3)
                                % Calculate the length of the element
                                length_element = norm(NodeCoord(node1,:) - NodeCoord(node2,:));
                                % Check if the length is less than or equal to the limit
                                if length_element <= limit
                                    plot3([NodeCoord(node1, 1), NodeCoord(node2, 1)], [NodeCoord(node1, 2), NodeCoord(node2, 2)], [NodeCoord(node1, 3), NodeCoord(node2, 3)], 'b');
                                    hold on;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Connect missing horizontal nodes
    for i = 1:size(NodeCoord, 1)
        for j = i+1:size(NodeCoord, 1)
            % Check if the nodes are on the same horizontal plane
            if NodeCoord(i, 2) == NodeCoord(j, 2)
                % Calculate the length of the segment
                length_segment = norm(NodeCoord(i,:) - NodeCoord(j,:));
                % Plot the segment only if its length is less than or equal to the limit
                if length_segment <= limit
                    plot3([NodeCoord(i, 1), NodeCoord(j, 1)], [NodeCoord(i, 2), NodeCoord(j, 2)], [NodeCoord(i, 3), NodeCoord(j, 3)], 'b');
                    hold on;
                end
            end
        end
    end
xlabel('X-axis: WIDTH(meter)');
ylabel('Y-axis: DEPTH(meter)');
zlabel('Z-axis: HEIGHT(meter)');
title(['Node and Connectivity']);
view(views(vp, :)); % Set the view angle
exportgraphics(gcf, ['Node_and_Connectivity_1.jpg']); % Save the figure as a JPEG image
axis equal
grid on
hold off;

end
end