function bar_eigenvector(coefs,dim,field_names)
%bar_eigenvector Plots the dim-th eigenvector of coefs as barh

       set(groot,'defaultAxesTickLabelInterpreter','none');  
       eig = coefs(:,dim);
       figure 
       barh(eig)
       set(gca,'YTick', 1:size(coefs,1),'yticklabel',field_names)
       set(gcf,'color','w');
       title(strcat("Content of the component ", mat2str(dim)))
       
end

