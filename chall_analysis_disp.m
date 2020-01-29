% close all
clear
clc

folders = {'Star1NoNoise_results_ss_mtf2x','Star2Noise_results_ss_mtf2x','Star3NoNoiseStrain_results_ss_mtf2x',...
    'Star4NoiseStrain_results_ss_mtf2x','Star5LargeNoisy_results_ss_mtf2x','Star6StrainNoisy_results_ss_mtf2x'};

cut_line = 257;
for jj = 5
    cnt = 0;
    for ii = [8,16,24,32,40,48,56,64]%[8,16,32,40,48,56]%[9,19,29,39,49,59]
                    
            file_in = ['.',filesep,folders{jj},filesep,strcat('results_qDIC_ss',num2str(ii),'.mat')];
            file_in_img = ['.',filesep,folders{jj},filesep,strcat('results_qDIC_ss',num2str(ii),'_img.mat')];
            load(file_in)
            img_mat = load(file_in_img);
            
            [newGridX,newGridY] = meshgrid([gridPoints{1,2}(1,1):gridPoints{1,2}(1,end)],[gridPoints{1,1}(1,1):gridPoints{1,1}(end,1)]);

        for kk = 1:2
            cnt = cnt + 1;
            
            u_y = interp2(gridPoints{1,2},gridPoints{1,1},u{kk}{2},newGridX,newGridY,'cubic');
            
            exp_data{jj}(cnt,:) = u_y(cut_line,1:end);
            
            mask = ~isnan(u_y(7:end-6,2:end));
            mask = 0.7*mask;
            
            figure
            ax1 = axes;
            image(ax1,img_mat.cellIMG{kk+1})
            colormap(ax1,gray(256))
            axis image
            
            ax2 = axes;
            colormap(ax2,jet(256))
            imagesc(ax2,u_y(7:end-6,2:end),'AlphaData',mask)
            colormap(ax2,'jet')
            axis image
            
            if kk == 1
                caxis(ax2, [-0.03,0.03])
                flag = 'noise';
            else
                caxis(ax2, [-0.5,0.5])
                flag = 'deformed';
            end
            
            ax1.Visible = 'off';
            ax2.Visible = 'off';
            linkprop([ax1 ax2],'Position');
            colorbar
            
            set(gcf,'Units','Inches');
            set(gcf,'Position',[3 3 12 9]);
            
            if jj == 1
                type = 'NoNoise';
            elseif jj == 2
                type = 'Noise';
            else
                type = 'Large';
            end
            
            saveas(gcf,sprintf('displacement_v_image%s_%s_SS%i.png',type,flag,ii))
        end
    end
end



%%
figure
plot(exp_data{jj}')
xlabel('pixel pos')
ylabel('disp ampl')


%%
% figure
% imagesc(u{1}{2}(12:end-11,12:end-11)),colorbar,axis image,caxis([-0.5,0.5])


