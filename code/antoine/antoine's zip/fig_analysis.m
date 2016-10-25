im_pattern (:,:,1) = (imread('Patterns\Tpic.jpg'));
im_pattern (:,:,2) = (imread('Patterns\Oblikpic.jpg'));
im_pattern (:,:,3) = (imread('Patterns\squares2.jpg'));
im_pattern (:,:,4) = (imread('Patterns\2_horiz2.jpg'));
im_pattern (:,:,5) = (imread('Patterns\Triangles_aligned.jpg'));
im_pattern (:,:,6) = (imread('Patterns\Triangles_same_COM.jpg'));
im_pattern (:,:,7) = (imread('Patterns\2bars_aligned.jpg'));
im_pattern (:,:,8) = (imread('Patterns\2bars_same_COM.jpg'));

nb_im = size (im_pattern,3);
for i = 1: nb_im;
[vals_l_r2(:,:,i),vals_r_r2(:,:,i),ths]=vf_panoconv_average_centre(im_pattern(:,:,i),'r2');
[vals_l_r4(:,:,i),vals_r_r4(:,:,i),ths]=vf_panoconv_average_centre(im_pattern(:,:,i),'r4');
end
vals_r2 = [vals_l_r2;vals_r_r2];
vals_r4 = [vals_l_r4;vals_r_r4];

   
   %Get the individuals activity shifted along the figure take only the left ones
    figure(); 
 nb_im = size (im_pattern,3);
 for i = 1:nb_im; %for each image
  subplot (3,nb_im,i); imshow (im_pattern(:,:,i));
  shifted_r4 = vals_l_r4(:,:,i);
  for j = 1:size (vals_l_r4,1);
       shifted_r4 (j,:) = shifted_r4(j,:)+(j*0.4);
  end
  shifted_r2 = vals_l_r2(:,:,i);
  for j = 1:size (vals_l_r2,1);
       shifted_r2 (j,:) = shifted_r2(j,:)+(j*0.4);
  end
  
  subplot (3,nb_im,i+nb_im); plot(ths',shifted_r4 (:,:)');hold on;
  axis ([-180 180 -1 6]); title('r4');
  subplot (3,nb_im,i+(nb_im*2)); plot(ths',shifted_r2 (:,:)');
  axis ([-180 180 -1 6]); title('r2');
 end
 
 
 %Get the difference when facing one pattern (36) vs the other (18);
  
%  nb_im = size (im_pattern,3);
%  for i = 1:nb_im;
%      diff_r2(:,i) = abs(vals_r2(:,36,i)-vals_r2(:,18,i));
%      diff_r4(:,i) = abs(vals_r4(:,36,i)-vals_r4(:,18,i));
%  end
%  figure()
%  barwitherr([std(diff_r2);std(diff_r4)]', [mean(diff_r2);mean(diff_r4)]');
%  
 
 Q1=[27:45]; Q2=[9:27];
  for i = 1:nb_im;
     dif_r2 = abs(vals_r2(:,Q1,i)-vals_r2(:,Q2,i));
     dif_r4 = abs(vals_r4(:,Q1,i)-vals_r4(:,Q2,i));
     
     diff_r2(:,i)= mean(dif_r2')';
     diff_r4(:,i)= mean(dif_r4')';
  end
figure()
 subplot (2,1,1) ; barwitherr([std(diff_r2);std(diff_r4)]', [mean(diff_r2);mean(diff_r4)]');
 title ('r2 then r4');
 subplot (2,1,2) ; barwitherr([std(diff_r2)]', [mean(diff_r2)]');
  title ('r2 only');

  
%     nb_im = size (im_pattern,3); 
%     for i = 1: nb_im;
%     diff_quadrant_v6(:,i) = compare_quadrant (vals_r2(:,:,i));
%     diff_quadrant_v4(:,i) = compare_quadrant (vals_r4(:,:,i));
%   end %Carefull as image now along dim 2 (kernel along dim 1)
%   subplot (2,1,1); boxplot (diff_quadrant_v6); title ('Activity difference between opposite and adjacent quadrant r2')
%   subplot (2,1,2); boxplot (diff_quadrant_v4); title ('r4')
%   
%     barwitherr([std(diff_quadrant_v6);std(diff_quadrant_v4)]', [mean(diff_quadrant_v6);mean(diff_quadrant_v4)]');
%   

%FIG1 2 std deviation for Tpics
im_pattern (:,:,1) = (imread('Patterns\Tpic.jpg'));
i = 1 ; %nb_im; %for each image
  k=1

    mins =[]; maxs =[];
  for ker = 1:size(vals_r6_pattern,1); %for each kernel find min and man of curve
    f = find (vals_r6_pattern (ker,:,i)> mean (vals_r6_pattern(ker,:,i)) + (std(vals_r6_pattern (ker,:,i))));%spot one std above average
    maxs = [maxs f];
    f = find (vals_r6_pattern (ker,:,i)< mean (vals_r6_pattern(ker,:,i)) - (std(vals_r6_pattern (ker,:,i))));%spot one std below average
    mins = [mins f];
  end
          n_max_r6 = histc (maxs,[1:length(vals_r6_pattern)]); 
          n_min_r6 = histc (mins,[1:length(vals_r6_pattern)]); 
       
           mins =[]; maxs =[];
  for ker = 1:size(vals_r4_pattern,1); %for each kernel find min and man of curve
    f = find (vals_r4_pattern (ker,:,i)> mean (vals_r4_pattern(ker,:,i)) + (std(vals_r4_pattern (ker,:,i))));%spot one std above average
    maxs = [maxs f];
    f = find (vals_r4_pattern (ker,:,i)< mean (vals_r4_pattern(ker,:,i)) - (std(vals_r4_pattern (ker,:,i))));%spot one std below average
    mins = [mins f];
  end
                  n_max_r4 = histc (maxs,[1:length(vals_r4_pattern)]); 
                  n_min_r4 = histc (mins,[1:length(vals_r4_pattern)]); 
                  
  subplot (3,1,1); 
  imshow (im_pattern(:,:,i)); hold on;

  subplot (3,1,2);
 plot (medfilt1(n_max_r6,10),'r');xlim([1 length(vals_r6_pattern)]);ylim([0 10]);hold on;
 plot (n_min_r6,'b');xlim([1 length(vals_r6_pattern)]);ylim([0 10]);hold on;
 title 'r2 '
 
  subplot (3,1,3);
 plot (n_max_r4,'r');xlim([1 length(vals_r4_pattern)]);ylim([0 5]);hold on;
 plot (n_min_r4,'b');xlim([1 length(vals_r4_pattern)]);ylim([0 5]);hold on;
 title 'r4'
 
 