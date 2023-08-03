function SSIM = satu(SIM, factor)
SSIM = saturate(SIM./max(SIM(:)).*inv_saturate(factor));
SSIM = SSIM ./ max(SSIM(:));
end