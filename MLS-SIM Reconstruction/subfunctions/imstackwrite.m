function imstackwrite(filename,imstack)
filename = char(filename);
imwrite(imstack(:,:,1),filename, 'WriteMode', 'overwrite');
for ii=2:size(imstack,3)
    imwrite(imstack(:,:,ii),filename,'WriteMode','append');
end
end