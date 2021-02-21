function [fval,scaled_image,scaled_cnvs] = cost_fcn_v3(cnvs,im,grid_num)
    persistent scaled_im mask grid_pxl_size
    if nargin == 3
        grid_pxl_size = round(size(im,1)/grid_num);
        scaled_im = blkproc(im, [grid_pxl_size grid_pxl_size], 'mean2');       
        R = round(size(scaled_im,1)/2);
        [X,Y] = meshgrid(1:size(scaled_im,1));
        mask = logical(((X-R).^2+(Y-R).^2)>R^2);
        scaled_im(mask) = 0;
        scaled_im = scaled_im./sum(scaled_im(:));
        scaled_image = scaled_im;
        fval = 0;
        return;
    end
    
    scaled_cnvs = blkproc(cnvs, [grid_pxl_size grid_pxl_size], 'mean2');
    scaled_cnvs(mask) = 0;
    scaled_cnvs = scaled_cnvs./sum(scaled_cnvs(:));
	fval = sum((scaled_im(:)-scaled_cnvs(:)).^2);
    
    scaled_image = scaled_im;
    
%     figure(666); clf;
%     subplot(1,2,1);
%     imshow (scaled_im);
%     subplot(1,2,2);
%     imshow (scaled_cnvs);
%     title(['f=' num2str(fval)]);
%     pause(0.5);
end