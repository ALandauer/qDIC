% close all
clear
clc

folders = {'Star1NoNoise_results_ss_mtf2x','Star2Noise_results_ss_mtf2x','Star3NoNoiseStrain_results_ss_mtf2x',...
    'Star4NoiseStrain_results_ss_mtf2x','Star5LargeNoisy_results_ss_mtf2x','Star6StrainNoisy_results_ss_mtf2x'};
ss_sw_pairs  = {[1,1],[2,3],[3,4],[4,5],[5,6],[6,7]};
SS = [8,16,24,32,40,48,56,64];
step = 1;
SW = [5:2:19];

cut_line = 257;
for jj = [3]
    cnt = 0;
    for ii = 1:6%[8,16,24,32,40,48,56,64]%[8,16,32,40,48,56]%[9,19,29,39,49,59]
        
        fprintf('SS%i, SW %i\n',SS(ss_sw_pairs{ii}(1)), SW(ss_sw_pairs{ii}(2)))
        
        file_in = ['.',filesep,folders{jj},filesep,strcat('results_qDIC_ss',num2str(SS(ss_sw_pairs{ii}(1))),'.mat')];
        file_in_img = ['.',filesep,folders{jj},filesep,strcat('results_qDIC_ss',num2str(SS(ss_sw_pairs{ii}(1))),'_img.mat')];
        load(file_in)
        img_mat = load(file_in_img);
        
        [newGridX,newGridY] = meshgrid([gridPoints{1,2}(1,1):gridPoints{1,2}(1,end)],[gridPoints{1,1}(1,1):gridPoints{1,1}(end,1)]);
        for kk = 1:1
            u_y{kk}{1} = interp2(gridPoints{1,2},gridPoints{1,1},u{kk}{1},newGridX,newGridY,'cubic');
            u_y{kk}{2} = interp2(gridPoints{1,2},gridPoints{1,1},u{kk}{2},newGridX,newGridY,'cubic');
        end
        
        m2px = [1,1];
        [Fij, J] = calculateFij_2d(u,dm,m2px,['optimal',num2str(SW(ss_sw_pairs{ii}(2)))]);
        [Eij_cell, eij_cell] = calculateEij_2d(Fij);
        
        for kk = 1:1
            cnt = cnt + 1;
            e_22{kk} = 100*interp2(gridPoints{1,2},gridPoints{1,1},Eij_cell{kk}{2,2},newGridX,newGridY,'cubic');
            exp_data_strain{jj}(cnt,:) = e_22{kk}(cut_line,1:end);
            
            mask = ~isnan(e_22{kk}(7:end-6,9:end-9));
            mask = 0.7*mask;
            
            figure
            ax1 = axes;
            image(ax1,img_mat.cellIMG{kk+1})
            colormap(ax1,gray(256))
            axis image
            
            ax2 = axes;
            colormap(ax2,jet(256))
            imagesc(ax2,e_22{kk}(7:end-6,9:end-9),'AlphaData',mask)
            colormap(ax2,'jet')
            axis image
            
            if kk ~= 1
                caxis(ax2, [-0.1,0.1])
                flag = 'noise';
            else
                caxis(ax2, [-5,5])
                flag = 'deformed';
            end
            
            ax1.Visible = 'off';
            ax2.Visible = 'off';
            linkprop([ax1 ax2],'Position');
            colorbar
            
            set(gcf,'Units','Inches');
            set(gcf,'Position',[3 3 12 9]);
            
            if jj == 3
                type = 'NoNoise';
            elseif jj == 4
                type = 'small';
            else
                type = 'large';
            end
            saveas(gcf,sprintf('strain_eyy_image%s_%s_SS%i_SW%i.png',type,flag,SS(ss_sw_pairs{ii}(1)), SW(ss_sw_pairs{ii}(2))))
        end
        
    end
end
% %
% % SS = [8,16,24,32,40,48,56,64];
% % step = 1;
% % SW = [5:2:19];
% %
% % all = 0;
% % for ii = 1:length(SS)
% %
% %    for jj = 1:length(SW)
% %         all = all+1;
% %         VSG(all) = ((SW(jj)-1)*step) + SS(ii);
% %     end
% % end

%%
figure
plot(exp_data_strain{jj})
xlabel('pixel pos')
ylabel('disp ampl')


%%
figure
imagesc(e_22{1}),colorbar,axis image,caxis([-0.5,0.5])


