function [select, coord] = getroi(array)
% [select, coord] = getroi(array)  Select and return region of interest
%
% Select and return image region of interest (ROI) via a GUI window and
% 'right-button-mouse' operation. If the mouse button is clicked and
%  released without movement, the entire displayed image will be selected.
%   array  (uint8)  - input image array(nlin, npix [, ncolor])
%   select (double) - output ROI as an array(newlin, newpix[, ncolor])
%   coord is list of coordinates of ROI (upperleft(x,y),lowerright(x,y))
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

dim = size(array);
nlin = dim(1);
npix = dim(2);
if size(dim)==[1 2];
  ncol =1;
 else;
  ncol = dim(3);
 end;
 screen = 0.95*(get(0, 'ScreenSize')); 
% 0.95 is to allow a tolerance which I need so very large narrow images stay
% visible

% screen = get(0, 'ScreenSize');
  
% Set aspect ratio approx to that of array
rat = npix/nlin;
if rat<0.25     % following lines make ROI selection of narrow images easier
 rat=0.25;
else if rat>4
 rat =4
end;
end;
if nlin>=npix;
   if nlin > 0.5*screen(4);
    ht = min(nlin, 0.8*screen(4));    % This change helps with large images
   else;
    ht = 0.5*screen(4);
   end;
   wid = round(ht*rat);
 else;
   if npix > 0.5*screen(3);
    wid = min(npix, 0.8*screen(3));   % This change helps with large images
   else;
    wid = 0.5*screen(3);
   end;
   ht = round(wid/rat);
end;
  pos = round([screen(3)/10 screen(4)/10 wid ht]); 
figure(1), set(gcf,'Position', pos);

disp(' ');
disp('Select ROI with right mouse button, no move = all');
disp(' ');

temp = class(array);
if temp(1:5) ~= 'uint8' 
  imagesc( double(array)/double(max(max(max(array))))), 
	colormap('gray'),
	title('Select ROI');
  else

if ncol == 1;
    imagesc(array), 
	colormap('gray'),
	title('Select ROI');
  else
    imagesc(array), 
	colormap('gray'),
    title('Select ROI');  
  end
end
%axis off
 
junk=waitforbuttonpress;
ul=get(gca,'CurrentPoint');
final_rect=rbbox;
lr=get(gca,'CurrentPoint');
ul=round(ul(1,1:2));
lr=round(lr(1,1:2));

if ul(1,1) > lr(1,1);             % sort x coordinates
   mtemp = ul(1,1);
   ul(1,1) = lr(1,1);
   lr(1,1) = mtemp;
 end;
 if ul(1,2) > lr(1,2);             % sort y coordinates
   mtemp = ul(1,2);
   ul(1,2) = lr(1,2);
   lr(1,2) = mtemp;
 end; 

roi = [lr(2)-ul(2)  lr(1)-ul(1)];  % if del x,y <10 pixels, select whole array
if roi(1)<10;
    ul(2) =1;
    lr(2) =nlin;
 end;
 if roi(2)<10;
    ul(1) =1;
    lr(1) =npix;
end;
select=double( array(ul(2):lr(2), ul(1):lr(1), :) );
coord = [ul(:,:), lr(:,:)];
close;
