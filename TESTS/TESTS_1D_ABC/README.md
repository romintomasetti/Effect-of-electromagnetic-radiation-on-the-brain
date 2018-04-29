This directory contains input files that are usefull for running the 12 possible 1D cases.

The convention for normal to surface is 'unit *outward* normal'.

The convention for file names is the following:

	1) FACE_EX or FACE_EY or FACE_EZ:
		The emitting face has a normal e_x, e_y or e_z, respectively.

	2) FACE_Minus_EX or FACE_Minus_EY or FACE_Minus_EZ:
		The emitting face has a normal -e_x, -e_y or -e_z, respectively.

	3) Electric_along_ + a direction (either X, Y or Z):
		Specifies the direction of the electric field.

The 12 cases are:

	1) Propagation in the +x direction, emitting face with normal -e_x:

		1.a) Electric field along y and magnetic field along z:
			IMPOSED=FACE_Minus_EX_Electric_along_Z

		1.b) Electric field along z and magnetic field along y:
			IMPOSED=FACE_Minus_EX_Electric_along_Y

	2) Propagation in the -x direction, emitting face with normal +e_x:

		2.a) Electric field along y and magnetic field along z:
			IMPOSED=FACE_EX_Electric_along_Z

		2.b) Electric field along z and magnetic field along y.
			IMPOSED=FACE_EX_Electric_along_Y

	3) Propagation in the +y direction, emtting face with normal -e_y:
	
		3.a) Electric field along z, magnetic field along x:
			IMPOSED=FACE_Minus_EY_Electric_along_Z

		3.b) Electric field along x, magnetic field along z:
			IMPOSED=FACE_Minus_EY_Electric_along_X

	4) Propagation in the -y direction, emitting face with normal +e_y:

		4.a) Electric field along z, magnetic field along x:
			IMPOSED=FACE_EY_Electric_along_Z

		4.b) Electric field along x, magnetic field along z:
			IMPOSED=FACE_EY_Electric_along_X

	5) Propagation in the +z direction, emitting face normal -e_z:
	
		5.a) Electric field along y, magnetic field along x:
			IMPOSED=FACE_Minus_EZ_Electric_along_Y

		5.b) Electric field along x, magnetic field along y:
			IMPOSED=FACE_Minus_EZ_Electric_along_X

	6) Propagation in the -z direction, emitting face normal +e_z:
	
		5.a) Electric field along y, magnetic field along x:
			IMPOSED=FACE_EZ_Electric_along_Y

		5.b) Electric field along x, magnetic field along y:
			IMPOSED=FACE_EZ_Electric_along_X


