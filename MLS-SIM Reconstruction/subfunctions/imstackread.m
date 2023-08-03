function imstack=imstackread(filename, frame_idc)
filename = char(filename);
info=imfinfo(filename);
H = info(1).Height;
W = info(1).Width;
D = size(info, 1);
if info(1).BitDepth==8
    img_type = 'uint8';
else
    img_type = 'uint16';
end

if nargin < 2
    frame_idc = 1:D;
else
    frame_idc = frame_idc(:);
end
assert(max(frame_idc) <= D);
D = numel(frame_idc);

imstack=zeros(H, W, D, img_type);
img = Tiff(filename, 'r');
for nD = 1:D
    img.setDirectory(frame_idc(nD));
    imstack(:, :, nD) = img.read();
end
disp(['load: ' filename ': finished']);
end