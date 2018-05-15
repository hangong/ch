function fusion = imShowPair(im1,im2)

sz = size(im1);
fusion = zeros([sz,3]);

fusion(:,:,2) = im1/max(im1(:));
fusion(:,:,1) = im2/max(im2(:));
fusion(:,:,3) = fusion(:,:,1);

imshow(fusion);

end
