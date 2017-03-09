im_pattern (:,:,1) = (imread('Patterns\Tpic.jpg'));
im_pattern (:,:,2) = (imread('Patterns\Bpic.jpg'));
im_pattern (:,:,3) = (imread('Patterns\Oblikpic.jpg'));
im_pattern (:,:,4) = (imread('Patterns\Oblikpic_20deg_up.jpg'));
im_pattern (:,:,5) = (imread('Patterns\squares2.jpg'));
im_pattern (:,:,6) = (imread('Patterns\Triangles_aligned.jpg'));
im_pattern (:,:,7) = (imread('Patterns\Triangles_same_COM.jpg'));
im_pattern (:,:,8) = (imread('Patterns\2bars_aligned.jpg'));
im_pattern (:,:,9) = (imread('Patterns\2bars_same_COM.jpg'));
im_pattern (:,:,10) = (imread('Patterns\2_horiz2.jpg'));
im_pattern (:,:,11) = (imread('Patterns\2_verti2.jpg'));

%im_pattern=im_nat;
%im_pattern=im_bars;
 nb_im = size (im_pattern,3);
for i = 1: nb_im;
[vals_l_r2(:,:,i),vals_r_r2(:,:,i),ths]=vf_panoconv_average_centre(im_pattern(:,:,i),'r2');
[vals_l_r4(:,:,i),vals_r_r4(:,:,i),ths]=vf_panoconv_average_centre(im_pattern(:,:,i),'r4');
end
vals_r2 = [vals_l_r2;vals_r_r2];
vals_r4 = [vals_l_r4;vals_r_r4];


% 
%   subplot (3,nb_im,i); imshow (im_pattern(:,:,i));
%   subplot (3,nb_im,i+nb_im); plot(ths',vals_l (:,:,i)');hold on;
%   title('l');
%   subplot (3,nb_im,i+(nb_im*2)); plot(ths',vals_r (:,:,i)');hold on;
%   title('r');

  %----------------rIDF----------------
  nb_im = size (im_pattern,3);
  for i = 1:nb_im;
  ref_r2 = vals_r2 (:,19,i); %Reference while facing middle of the pics.
  ref_r4 = vals_r4 (:,19,i);
    
  diff_r2=[]; diff_r4=[];
    for j = 1:size(vals_r2,2); %For each azimutal pixel
        diff_r2(j) = mean(abs(vals_r2(:,j,i)-ref_r2));
        diff_r4(j) = mean(abs(vals_r4(:,j,i)-ref_r4));
    end
    rIDF_r2 (i,:)=diff_r2;
    rIDF_r4 (i,:)=diff_r4;
  end
  figure();
  for i = 1:nb_im;
  subplot (3,nb_im,i);imshow (im_pattern(:,:,i));
  subplot (3,nb_im,i+nb_im); plot (ths,rIDF_r2(i,:)); xlim([-180 180]);title ('r2'); 
  subplot (3,nb_im,i+(nb_im*2)); plot (ths,rIDF_r4(i,:));xlim([-180 180]);title ('r4'); 
  end
      
  %---------rIDF only left eyes-------
  
    nb_im = size (im_pattern,3);
  for i = 1:nb_im;
  ref_r2 = vals_l_r2 (:,37,i); %Reference while facing middle of the pics.
  ref_r4 = vals_l_r4 (:,37,i);
    
  diff_r2=[]; diff_r4=[];
    for j = 1:size(vals_r2,2); %For each azimutal pixel
        diff_r2(j) = mean(abs(vals_l_r2(:,j,i)-ref_r2));
        diff_r4(j) = mean(abs(vals_l_r4(:,j,i)-ref_r4));
    end
    rIDF_l_r2 (i,:)=diff_r2;
    rIDF_l_r4 (i,:)=diff_r4;
  end
%   figure();
%   for i = 1:nb_im;
%   subplot (3,nb_im,i);imshow (im_pattern(:,:,i));
%   subplot (3,nb_im,i+nb_im);     plot (ths,rIDF_l_r2(i,:)); xlim([-180 180]);title ('r2'); 
%   subplot (3,nb_im,i+(nb_im*2)); plot (ths,rIDF_l_r4(i,:));xlim([-180 180]);title ('r4'); 
%   end
  
  %Only Right eye
  nb_im = size (im_pattern,3);
  for i = 1:nb_im;
  ref_r2 = vals_r_r2 (:,37,i); %Reference while facing middle of the pics.
  ref_r4 = vals_r_r4 (:,37,i);
    
  diff_r2=[]; diff_r4=[];
    for j = 1:size(vals_r2,2); %For each azimutal pixel
        diff_r2(j) = mean(abs(vals_r_r2(:,j,i)-ref_r2));
        diff_r4(j) = mean(abs(vals_r_r4(:,j,i)-ref_r4));
    end
    rIDF_r_r2 (i,:)=diff_r2;
    rIDF_r_r4 (i,:)=diff_r4;
  end

  figure(); %Both eyes
  for i = 1:nb_im;
  
  subplot (3,nb_im,i);imshow (im_pattern(:,:,i));
  
  subplot (3,nb_im,i+nb_im); 
  plot (ths,rIDF_l_r2(i,:)); hold on;
  plot (ths,rIDF_r_r2(i,:),'r'); xlim([-180 180]);title ('r2'); 
  
  subplot (3,nb_im,i+(nb_im*2)); 
  plot (ths,rIDF_l_r4(i,:)); hold on;
  plot (ths,rIDF_r_r4(i,:),'r');xlim([-180 180]);title ('r4'); 
  
  end
  
  
  %------Test for retinal position invariance-------
  
  im_oblik = [3 4];
  
  ref_r2 = vals_r2 (:,37,3); %Reference image oblik centred (3).
  ref_r4 = vals_r4 (:,37,3);
  
  k=0;
  for i = im_oblik; %Test correct image (3) then the elvated one (4)
  k=k+1;
      diff_r2=[]; diff_r4=[];
    for j = 1:size(vals_r2,2); %For each azimutal pixel
        diff_r2(j) = mean(abs(vals_r2(:,j,i)-ref_r2));
        diff_r4(j) = mean(abs(vals_r4(:,j,i)-ref_r4));
    end
    rIDF_oblik_r2 (k,:)=diff_r2;
    rIDF_oblik_r4 (k,:)=diff_r4;

  subplot (3,2,k);imshow (im_pattern(:,:,i));
  subplot (3,2,k+2); plot (ths,rIDF_oblik_r2(k,:)); xlim([-180 180]);title ('r2'); 
  subplot (3,2,k+(2*2)); plot (ths,rIDF_oblik_r4(k,:));xlim([-180 180]);title ('r4'); 
  end

