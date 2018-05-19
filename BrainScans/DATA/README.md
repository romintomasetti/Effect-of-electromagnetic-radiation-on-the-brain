Source:
	http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?alias=subject20_crisp&download=1
	[ACCESSED ON May 4 2018]


####################################################################################
# What is inside this directory ?

	1) subject20_crisp_v.rawb
		This is the file downloaded from the link above.
		It consists in a MRI of a brain.

	2) subject20_crisp_v_rawb.vti
		This is the previous raw bytes file, converted in .vti format for visualization in VTI format.

	3) materials_100by119by100.vti
		This is the result of a resampling (downsampling) of the file subject20_crisp_v_rawb.vti with Paraview.

	4) materials_100by119by100.bin
		This is the file materials_100by119by100.vti in binary format.
	
	5) materials_104by123by104.bin
		This is the file resulting from cleaning of materials_100by119by100.bin in MATLAB.
		This is the file we use for simulations.