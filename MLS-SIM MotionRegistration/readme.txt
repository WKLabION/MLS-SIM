Prerequisites:
PC with Windows 10 operating system
At least 16GB RAM
MATLAB version R2021a or higher

Installation guide:
Install all the software in the prerequisites in order. The typical install time on a "normal" desktop computer is about 30 minutes.

Expected output and runtime:
Motion-corrected image stacks from raw image stack sequences. The typical runtime varies for different data size.

Instructions for use:
(1) run <template_gen_final.m> to generate a template image <Template.tif>
(2) run <line_registration_final.m> to generate image registration files <XXX_reg_13_lines.mat>
(3) run <frame_crt_final.m> to generate motion corrected images
 
