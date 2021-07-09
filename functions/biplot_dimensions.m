function biplot_dimensions(coefs,score,dims,oblabels)
%Generate colored biplot For specific data

    figure
    hbi = biplot(coefs(:,dims),'Scores',score(:,dims),'Marker','o','ObsLabels',oblabels);
    for ii = 1:length(hbi)-length(oblabels)
    set(hbi(ii), 'Color', [0.5 0.5 0.5]);
    end
    
    
    
    for ii = length(hbi)-length(oblabels):length(hbi)
        userdata = get(hbi(ii), 'UserData');
        if ~isempty(userdata)
            if oblabels(userdata) == "Healthy"
                color = [0 0.75 0];
            elseif oblabels(userdata) == "EES"
                color = [0 0 0.75];
            end
            if oblabels(userdata) == "NoEES"
                color = [0.75 0 0];
            end
            if oblabels(userdata) == "Human" 
                color = [0 0.6 0.6];
            end
            if oblabels(userdata) == "NHP" 
                color =  [0.6 0 0.6];
            end
            if oblabels(userdata) == "Rodent" 
                color =  [0.6 0.6 0];
            end

            set(hbi(ii), 'Color', color);
            set(hbi(ii), 'MarkerFaceColor', color);

        else
            set(hbi(ii), 'Color', 'k');
        end
    end
    title(strcat("Plot along the eigenvectors ",mat2str(dims)));
    xlabel(strcat("Component ",mat2str(dims(1))));
    ylabel(strcat("Component ",mat2str(dims(2))));
    
    
    set(gcf,'color','w');

end

