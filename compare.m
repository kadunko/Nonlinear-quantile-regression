function compare
close all

for m = 1:2
    figure(m)
    if m == 1 %these numbers come from 'tables.pdf', attached
        mat1 = [1.2232	1.2212	1.4413	1.2551	2.1088;
            1.0061	1.0046	1.1774	1.0303	1.6827;
            0.5614	0.5608	0.6439	0.5734	0.8669;
            0.0972	0.0971	0.1032	0.0981	0.1245];
        
        mat2 = [1.2486	1.4271	1.2431	1.3290	1.4409;
            1.0230	1.1595	1.0188	1.0832	1.1694;
            0.5678	0.6285	0.5638	0.5905	0.6300;
            0.0946	0.0986	0.0945	0.0958	0.1016];
        
        mat3 = [1.2742	1.3005	1.3353	1.2722	1.7511;
            1.0424	1.0643	1.0927	1.0417	1.4088;
            0.5756	0.5871	0.6013	0.5754	0.7417;
            0.0967	0.0977	0.0987	0.0967	0.1126];
        
        mat4 = [1.0127	1.6829	1.0927	1.4566	0.9406;
            0.7875	1.3493	0.8969	1.1751	0.7798;
            0.4593	0.7069	0.4995	0.6218	0.4441;
            0.0820	0.0976	0.0858	0.0921	0.0816];
        
    elseif m == 2
        mat1 = [1.4888	1.2219	1.4552	1.2481	2.4613;
            1.1115	1.0060	1.1927	1.0279	1.8367;
            0.6174	0.5653	0.6565	0.5776	0.9273;
            0.1195	0.1171	0.1246	0.1180	0.1458];
        
        mat2 = [1.3728	1.4344	1.2434	1.3381	1.5969;
            1.0980	1.1666	1.0218	1.0906	1.2342;
            0.6054	0.6357	0.5677	0.5983	0.6619;
            0.1170	0.1189	0.1153	0.1157	0.12406];
        
        mat3 = [1.4197	1.3048	1.3417	1.2738	1.9919;
            1.1358	1.0689	1.1013	1.0446	1.5121;
            0.6110	0.5931	0.6091	0.5815	0.7864;
            0.1203	0.1178	0.1198	0.1166	0.1344];
        
        mat4 = [1.2186	1.6944	1.0840	1.4775	0.9719;
            0.9008	1.3595	0.8917	1.1884	0.7873;
            0.4887	0.7151	0.4974	0.6316	0.4543;
            0.1018	0.1188	0.1054	0.1121	0.1058];
    end
    
    y = linspace(.5, 4.5, 100);
    y = y';
    ymin = zeros(4,1);
    ymax = ymin;
    
    mat1 = mat1./mat1(:,1);
    mat1 = mat1(:,2:5);  
        ymin(1) = min(min(mat1)) - .1;
        ymax(1) = max(max(mat1)) + .1;  
    mat2 = mat2./mat2(:,1);
    mat2 = mat2(:,2:5);  
        ymin(2) = min(min(mat2)) - .1;
        ymax(2) = max(max(mat2)) + .1;      
    mat3 = mat3./mat3(:,1);
    mat3 = mat3(:,2:5);  
        ymin(3) = min(min(mat3)) - .1;
        ymax(3) = max(max(mat3)) + .1;      
    mat4 = mat4./mat4(:,1);
    mat4 = mat4(:,2:5);  
        ymin(4) = min(min(mat4)) - .1;
        ymax(4) = max(max(mat4)) + .1; 
     ylim1 = min(ymin);
     ylim2 = max(ymax);
    for k = 1:4
        subplot(1,4,k)
        if k==1
            labelstring = '$\sigma(x) \propto {1}$';
            mat = mat1;
        elseif k==2
            labelstring = '$\sigma(x) \propto {x+x_{N}}$';
            mat = mat2;
        elseif k==3
            labelstring = '$\sigma(x) \propto \sqrt{x+x_{N}}$';
            mat = mat3;
        elseif k==4
            labelstring = '$\sigma(x) \propto {x+.25*x_{N}}$';
            mat = mat4;
        end
        markers = ['o'; '+'; '*'; 'x'];
        c = ['k'; 'm'; 'b'; 'r'];
     hold on
        for j = 1:4
            for i = 1:4
                scatter(j, mat(i,j), c(i), markers(i));
            end
        end
        xlim([.5 4.5]);
        ylim([ylim1 ylim2]);
        plot(y, ones(length(y)));
        xticks([1 2 3 4]);
        xticklabels({'SL1', 'SL2', 'SL3', 'SL4'});
        xlabel(labelstring, 'Interpreter', 'latex');
        if k==2
            strA = {'$\nu = .20$'; '$\nu = .35$'; '$\nu = .65$'; '$\nu = .95$'};
            legend(strA,'Interpreter', 'latex', 'Location', 'northwest');
        end
    hold off
    end
end