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
for jj = 5
    cnt = 0;
    for ii = 1:6%[8,16,24,32,40,48,56,64]%[8,16,32,40,48,56]%[9,19,29,39,49,59]
        
        fprintf('SS%i, SW %i',SS(ss_sw_pairs{ii}(1)), SW(ss_sw_pairs{ii}(2)))
        
        file_in = ['.',filesep,folders{jj},filesep,strcat('results_qDIC_ss',num2str(SS(ss_sw_pairs{ii}(1))),'.mat')];
        load(file_in)
        
        [newGridX,newGridY] = meshgrid([gridPoints{1,2}(1,1):gridPoints{1,2}(1,end)],[gridPoints{1,1}(1,1):gridPoints{1,1}(end,1)]);
        for kk = 1:2
            u_y{kk}{1} = interp2(gridPoints{1,2},gridPoints{1,1},u{kk}{1},newGridX,newGridY,'cubic');
            u_y{kk}{2} = interp2(gridPoints{1,2},gridPoints{1,1},u{kk}{2},newGridX,newGridY,'cubic');
        end
        
        m2px = [1,1];
        [Fij, J] = calculateFij_2d(u,dm,m2px,['optimal',num2str(SW(ss_sw_pairs{ii}(2)))]);
        [Eij, eij_cell] = calculateEij_2d(Fij);
        
        for kk = 1:2
            cnt = cnt + 1;
            exp_data{jj}(cnt,:) = u_y(cut_line,1:end);
        end
        
    end
end

SS = [8,16,24,32,40,48,56,64];
step = 1;
SW = [5:2:19];

all = 0;
for ii = 1:length(SS)
 
   for jj = 1:length(SW)
        all = all+1;
        VSG(all) = ((SW(jj)-1)*step) + SS(ii);
    end
end

%%
figure
plot(exp_data{jj}')
xlabel('pixel pos')
ylabel('disp ampl')


%%
figure
imagesc(u_y),colorbar,axis image,caxis([-0.5,0.5])


