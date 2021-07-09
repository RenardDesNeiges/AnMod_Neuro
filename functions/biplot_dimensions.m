function biplot_dimensions(coefs,score,dims,oblabels)
%Generate colored biplot For specific data

    figure
    hbi = biplot(coefs(:,dims),'Scores',score(:,dims),'ObsLabels',oblabels);
    for ii = 1:length(hbi)-length(oblabels)
    set(hbi(ii), 'Color', [0.5 0.5 0.5]);
    end
    
    for ii = length(hbi)-length(oblabels):length(hbi)
    userdata = get(hbi(ii), 'UserData');
    if ~isempty(userdata)
        if oblabels(userdata) == "Healthy"
            set(hbi(ii), 'Color', [0 0.75 0]);
        elseif oblabels(userdata) == "EES"
            set(hbi(ii), 'Color', [0 0 0.75]);
        end
        if oblabels(userdata) == "NoEES"
            set(hbi(ii), 'Color', [0.75 0 0]);
        end
        if oblabels(userdata) == "Human" 
            set(hbi(ii), 'Color', [0 0.6 0.2]);
        end
        if oblabels(userdata) == "NHP" 
            set(hbi(ii), 'Color', [0.6 0 0.2]);
        end
    else
            set(hbi(ii), 'Color', 'k');
        end
    end
    
    set(gcf,'color','w');

end

